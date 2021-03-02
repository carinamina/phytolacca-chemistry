###Data wrangling for all herbivory datasets

###########
#setup
###########
setwd("/Users/carina/Documents/R_working_directory")
lats = read.csv("lats_long_names.csv", header=TRUE)
library(afex)
library(pbkrtest)
library(MuMIn)
library(arm)
library(effects)
library(lsmeans)
library(car)
library(nlme)
library(bbmle)
library(multcomp)
library(plyr)
library(ggplot2)
library(googlesheets)

lme_results <- function(lme_model)
{
  a<- hist((resid(lme_model) - mean(resid(lme_model), na.rm=T)) / sd(resid(lme_model), na.rm=T), freq=F); curve(dnorm, add = TRUE)
  b<-qqnorm(lme_model)
  c<-plot(lme_model)
  d<-summary(lme_model)
  #e<-intervals(lme_model)[1]
  f<-anova(lme_model, type = "marginal", test = "F")
  g<-nobs(lme_model)
  
  list(a,b,c,d,f,g)
}

lme_slopes <- function(lme_model)
{
  b_mature = summary(lme_model)$coef$fixed[1]
  m_mature = summary(lme_model)$coef$fixed[2]
  b_young = summary(lme_model)$coef$fixed[1]+summary(lme_model)$coef$fixed[3]
  m_young = summary(lme_model)$coef$fixed[2]+summary(lme_model)$coef$fixed[4]
  coef <- as.vector(c(b_mature,m_mature,b_young,m_young))
  names(coef) <- c("b_mature","m_mature","b_young","m_young")
  slopes <- as.data.frame(t(coef))
  return(slopes)
}

#############
#Palatability
#############
palat <- read.csv("palatability_PHRI_20180327.csv", header = TRUE)
sum(is.na(palat))
palat$X <- NULL
str(palat)
# cupID = cup identifier (experimental unit)
# line = maternal line of diet eaten in that cup
# lat = latitude of origin of plant diet
# pop = population of origin of plant diet
# age = leaf age
# cons_init_day = ln(leaf area consumed+1) per initial caterpillar per day
# mass_init_day = ln(final larval biomass+1) per initial caterpillar per day

popcheck <- unique(palat[c("pop","lat")])
#goofy but I'm having trouble figuring out sorting the df. Open this and sort by lat: the three southern-most and northern-most populations are the same as the populations used in other analyses: With, Hwy 27, MacArthur; Caesar, KBS, McPhee
rm(popcheck)

regions <- c("tropical", "tropical","tropical", "subtropical", "subtropical", "subtropical",NA,NA,NA,NA,NA,NA,NA, "temperate", "temperate", "temperate")
regions <- cbind(regions, sort(unique(palat$lat)))
colnames(regions) <- c("region", "lat")
regional <- droplevels(na.omit(merge(regions, palat, by = "lat")))

rm(regions)

#############
#Field herbivory
#############
wholedoc_trop <- gs_title("Tropical herbivory 2016_census 1-3_RAW_checked")
trop <- as.data.frame(gs_read_csv(wholedoc_trop, ws = "Sheet1", as.is = TRUE))
trop$season <- as.factor(trop$season)
trop$season <- revalue(trop$season, c("1" = "dry","2" = "dry-wet", "3" = "wet"))
#how many leaves were marked at each population at each timepoint
nrow(trop)/3/3
wholedoc_MI <- gs_title("herbivory 2017_census 1-2_checked")
mi <- as.data.frame(gs_read_csv(wholedoc_MI, ws = "Sheet1", as.is = TRUE))
nrow(mi)/3/2
nrow(trop)+nrow(mi)
mi$season <- as.factor(mi$season)
mi$season <- revalue(mi$season, c("1" = "mid","2" = "late"))
names(trop)
names(mi)

field <- rbind(trop, mi)
#str(field)
rm(trop, mi, wholedoc_MI, wholedoc_trop)

field$indiv <- as.factor(paste(field$pop, field$ID, sep = "_"))
field$ID <- NULL
field$pop <- as.factor(field$pop)
field <- droplevels(subset(field, field$pop != "BEN"))
#total marked leaves
dim(field)
#tropical marked
dim(subset(field, field$pop == "Tiri" | field$pop == "Gav" | field$pop == "Bella"))

field$age <- revalue(field$age, c("juv" = "young"))
field$age <- as.factor(field$age)
#remove those with NA date measured, which were dead or couldn't be found
field <- field[which(is.na(field$date_meas) == FALSE),]
field$date_mark <- as.Date(field$date_mark, "%m/%d/%Y")
field$date_meas <- as.Date(field$date_meas, "%m/%d/%Y")

field$interval <- as.numeric(field$date_meas - field$date_mark)

field$final <- rep(NA, nrow(field))
for(i in 1:nrow(field)){ 
  #subtract "damaged" areas (due to disease or whatever) from total area before calculating %herbiv
  if (is.na(field$meas_dam[i]) == FALSE)
    field$final[i] = field$meas_herbiv_area[i]/(field$meas_total_area[i]-field$meas_dam[i])*100
  #calculate %herbiv for partially eaten
  else if (is.na(field$meas_herbiv_area[i]) == FALSE)
    field$final[i] = field$meas_herbiv_area[i]/field$meas_total_area[i]*100
  #herbiv is 0 if uneaten or "not mature"
  else if (is.na(field$meas_no_herbiv[i]) == FALSE)
    field$final[i] = 0
  else if (is.na(field$nm[i]) == FALSE)
    field$final[i] = 0
  #herbiv is 100 if completely eaten or "meristem damage"
  else if (is.na(field$meas_100_herbiv[i]) == FALSE)
    field$final[i] = 100
  else if (is.na(field$md[i]) == FALSE)
    field$final[i] = 100
  #otherwise no data (e.g. if dead or "can't find leaf")
  else
    field$final[i] = NA
}

field$initial <- rep(NA, nrow(field))
for(i in 1:nrow(field)){
  if (is.na(field$mark_herbiv_area[i]) == FALSE)
    field$initial[i] = field$mark_herbiv_area[i]/field$mark_TOTAL_area[i]*100
  else if (is.na(field$mark_no_herbiv[i]) == FALSE)
    field$initial[i] = 0
  else
    field$initial[i] = NA
}

field$standing <- ifelse(field$age == "young", field$final, field$final - field$initial)
#some are negative bc of measurement error, but since we are taking mean per plant I think it's fine to leave them
field$rate <- field$standing/field$interval

field_indiv <- ddply(field, c("pop","indiv","season","age"), summarise,
                date_meas = first(date_meas),
                standing = mean(standing,na.rm=T),
                rate = mean(rate,na.rm=T),
                interval = mean(interval, na.rm=T)
)
field_indiv <- droplevels(na.omit(field_indiv))

hist(field_indiv$standing)
hist(field_indiv$rate)
#still some negatives but does it matter?

field_indiv$region <- as.factor(ifelse(field_indiv$pop == "Bella" | field_indiv$pop == "Gav" | field_indiv$pop == "Tiri", "tropical", "temperate"))
field_indiv$pop_season <- as.factor(paste(field_indiv$pop, field_indiv$season, sep = "_"))
field_indiv$pop_age <- as.factor(paste(field_indiv$pop, field_indiv$age, sep = "_"))

#total number of individuals measured
dim(field_indiv)
#number of tropical individuals measured
dim(subset(field_indiv, field_indiv$region == "tropical"))

rm(field, i, lats)

# pop = population
# indiv = individual plant
# season = which season (3 tropical, 2 MI)
# age = leaf age
# date_meas = date of final leaf measurement
# standing = final minus initial herbivory (range from 0-100%, but some negatives due to error)
# rate = daily herbivory rate (standing/interval from marking to measuring)
# region = tropical vs temperate (to match videos)
# pop_season = factor pasting population and season
# pop_age = factor pasting population and age
