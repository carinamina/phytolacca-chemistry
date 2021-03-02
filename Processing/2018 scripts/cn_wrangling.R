#chemistry data wrangling
#start apr 22, 2018
#combined the KBS and main campus data in Excel
library(dplyr)
library(plyr)
library(vegan)
library(lme4) 
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

lme_results <- function(lme_model)
{
  a<- hist((resid(lme_model) - mean(resid(lme_model), na.rm=T)) / sd(resid(lme_model), na.rm=T), freq=F); curve(dnorm, add = TRUE)
  b<-qqnorm(lme_model)
  c<-plot(lme_model)
  d<-summary(lme_model)
  #  e<-intervals(lme_model)[1]
  f<-anova(lme_model, type = "marginal", test = "F")
  g<-nobs(lme_model)
  
  list(a,b,c,d,f,g)
}

setwd("/Users/carina/Documents/R_working_directory")
cn = read.csv("all_CN_20180422.csv", header=TRUE)
str(cn)

# **pop** is the population  
# **line** is the maternal line ID nested within pop
# **indiv** is the plant ID letter nested within line  
# **age** is leaf age, young or mature  
# **percent_N** is the percentage of N in sample  
# **percent_C** is the percentage of C in sample  
# **C_N** is percent C/percent N for each sample
# 
# ####New variables  
# **line_ID** is a unique ID for each line with the population name, an underscore, and the line number  
# **plant_ID** is a unique ID for each plant with the population name, an underscore, the line number, an underscore, and the indiv letter  
# **lat** is latitude of each site (which is actually not necessary for this particular analysis, which uses region)  

cn$line_ID <- as.factor(paste(cn$pop,cn$line,sep="_"))
cn$plant_ID <- as.factor(paste(cn$line_ID,cn$indiv,sep="_"))
cn$pop_age <- as.factor(paste(cn$pop, cn$age, sep = "_"))

lat = read.csv("lats_long_names.csv", header = TRUE)
cn = merge(cn, lat)

#assigns regional names based on latitude. 3 tropical, 3 southernmost US, 3 central (based on 3 closest to mean of lowest and highest US) 3 northernmost US
regions <- c("tropical", "tropical", "tropical", "southern", "southern", "southern", NA, NA, NA, NA, NA, NA, NA, "northern", "northern", "northern")
regions <- cbind(regions, sort(unique(cn$lat)))
colnames(regions) <- c("region", "lat")
cn <- merge(cn, regions)

#get line-level C:N
cn_line <- ddply(cn, c("region","line_ID","pop","lat", "age"), summarise,
                C_N = mean(C_N),
                percent_N = mean(percent_N),
                percent_C = mean(percent_C) 
)

#add palatability data
cn_line$line_age <- as.factor(paste(cn_line$line_ID, cn_line$age, sep = "_"))
cn_line$line_ID <- NULL
palat = read.csv("line_level_palatability.csv", header = TRUE)
palat$line_age = paste(palat$line, palat$age, sep="_")
palat$X <- NULL
palat$pop <- NULL
palat$lat <- NULL
palat$age <- NULL
palat_cn <- merge(palat, cn_line, by = "line_age", all = TRUE)
str(palat_cn)
#remove anything that doesn't have palatability data
palat_cn <- palat_cn[!is.na(palat_cn$biomass),]
#biomass is ln(total biomass per cup) (counting dead cats as zero biomass), standardized by intitial cat count and duration
#area is ln(cumulative area consumed per cup) over course of experiment, standardized by intitial cat count and duration
#survival is number of survivors
rm(cn_line,lat,palat,regions)
