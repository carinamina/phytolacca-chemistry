#chemistry data wrangling
#start may 7, 2018 on re-aligned data
#did two steps in excel, such a cheater: deleted the pooled samples, all the normalised abundance columns and first two rows
#this time I did not subtract blanks. The Dec blanks didn't align properly with telmisartan and Tony didn't think it was a good idea.
library(dplyr)
library(plyr)
library(vegan)

#results function
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

#y = mx + b
#intercept = v1 = intercept for mature leaves
#lat = v2 = slope for mature leaves
#ageyoung = v3 = difference in intercepts between ages; intercept for young = v1+v3
#lat:ageyoung = v4 = difference in slopes between ages; slope for young = v2+v4

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

setwd("/Users/carina/Documents/R_working_directory")
raw = read.csv("20180504_export_all.csv", header=TRUE)

rownames(raw) <- raw[,1]
# #remove retention times > 7:30 (lipids that are just part of the membranes)
raw$retention <- as.numeric(read.table(text = rownames(raw), sep = "_", colClasses = "character")[,1])
raw <- droplevels(subset(raw, raw$retention < 7.5))
raw$retention <- NULL

#cut out extraneous columns, including blanks at the beginning
chem.df <- raw[,20:length(raw)]
rm(raw)

#transpose to columns as compounds
chem.df <- as.data.frame(t(chem.df))

sample <- rownames(chem.df)
month <- ifelse(grepl("CB", sample, fixed=TRUE),"Dec","Feb")
sample_month <- as.data.frame(cbind(sample, month))

feb_key <- read.csv("2018_feb_chem_key.csv", header=TRUE)
feb_key <- cbind(droplevels(subset(sample_month, sample_month$month == "Feb")), feb_key)
#check BY HAND that the LCMS sequence matches the sample ID

feb_key$LCMS.sequence <- NULL
feb_key$chem.label <- NULL
feb_key$dilution <- ifelse(feb_key$dilution == "1:01",1,20)
feb_key$dilution <- as.factor(feb_key$dilution)
#some of the bulk lines have the same number as non-bulk lines from those populations, so append a B to make them unique
feb_key$line <- ifelse(grepl("TT Bulk",feb_key$pop, fixed=TRUE) | grepl("Gav Bulk",feb_key$pop, fixed=TRUE), paste(feb_key$line, "B"),feb_key$line)
#Remove the bulk label so we know which population is which. Would be good to check someday these are no different
feb_key$pop <- revalue(feb_key$pop, c("Val" = "Tiri","Tiri Bulk" = "Tiri", "Hwy 27 bulk" = "Hwy 27", "TT Bulk" = "TT", "Dalton Bulk" = "Dalton", "Gav Bulk" = "Gav")) 
str(feb_key)

#do the same for Dec
dec_key <- read.csv("2017_dec_chem_key.csv", header=TRUE)
dec_key <- cbind(droplevels(subset(sample_month, sample_month$month == "Dec")), dec_key)
#check BY HAND that the LCMS sequence matches the sample ID

dec_key$LCMS.sequence <- NULL
dec_key$chem.label <- NULL
dec_key$region <- NULL
dec_key$dilution <- ifelse(dec_key$dilution == "1:1",1,20)
dec_key$dilution <- as.factor(dec_key$dilution)
dec_key$pop <- revalue(dec_key$pop, c("Val" = "Tiri"))
str(dec_key)

#merge Dec and Feb into master key
key <- rbind(dec_key, feb_key)
rm(dec_key,feb_key, sample_month, month, sample)
lat = read.csv("lats_long_names.csv", header = TRUE)
key <- merge(key, lat)

#add regional names
regions <- c("tropical", "tropical","tropical", "subtropical", "subtropical", "subtropical",NA,NA,NA,NA,NA,NA,NA, "temperate", "temperate", "temperate")
regions <- cbind(regions, sort(unique(key$lat)))
colnames(regions) <- c("region", "lat")
key <- merge(key, regions, by = "lat")
rm(regions,lat)

#replace rownames of chem data with a useful id: pos_age_dilution which matches something in master key
key$id <- as.factor(paste(key$pos, key$age, key$dilution, sep = "_"))
chem.df$sample <- rownames(chem.df)
chem.df <- merge(key, chem.df, by = "sample")
rownames(chem.df) <- chem.df$id

###############
#now we are done with identifying samples and their categorical characteristics. We need to normalize by telmisartin

#check for low telmisartin peaks
hist(chem.df$"4.78_514.2360n")
chem.df[which(chem.df$"4.78_514.2360n"<20000),1:11]
#33_mature_20 is lower than normal by an order of magnitude. Tony calculated the abundance according to MassLynx of this sample and the 33_young_20, and we used that ratio to determine that the telmisartan peak should be 77,900 for this sample. (We checked another peak and it seemed to have aligned normally, so it's not clear why Progenesis messed this one up but it doesn't seem to be a problem with the entire sample)
chem.df[which(chem.df$id == "33_mature_20"),]$"4.78_514.2360n" = 77900

#normalize by telmisartin: divide all peaks in each sample by telmisartan of that sample. It's 4.78_514.2360n
normal <- cbind(chem.df[,1:11],  chem.df[,12:ncol(chem.df)]/chem.df$"4.78_514.2360n")
normal$"4.78_514.2360n"

#look at distributions of total abundance
abund <- rowSums(normal[,12:nrow(normal)])
totalabund <- data.frame(normal[,1:11], abund)
hist(totalabund$abund)
hist(subset(totalabund, totalabund$dilution == "20")$abund)
hist(subset(totalabund, totalabund$dilution == "1")$abund)
summary(subset(totalabund, totalabund$dilution == "20")$abund)
summary(subset(totalabund, totalabund$dilution == "1")$abund)
#these seem ok this time, yay!

#reshape totalabund to wide, with a column for each dilution, and check that ratio of undiluted is ~~20x diluted. Should catch any really weird ratios
totalabund$pos_age <- as.factor(paste(totalabund$pos, totalabund$age, sep = "_"))
ratio <- reshape(totalabund,  timevar ="dilution", idvar = "pos_age", direction = "wide")
ratio$ratio <- ratio$abund.1/ratio$abund.20
hist(ratio$ratio)
#there is one with a pretty low ratio, 76_young. It seems like the problem is with the total abundance of the nondiluted sample; it is the lowest of all samples, even though the chromatograms look normal. So I think it's ok because we don't care about nondiluted.
rm(totalabund,ratio,abund,chem.df)


################
#to review, at this point we have removed all the membrane lipids (>7:30 retention), divided all samples by their own internal standard peak (normalized by telmisartin)
#we need to normalize by leaf mass and remove compounds that are zero for all samples.

#normalize by leaf mass
leaf <- read.csv("2017_chem_leaf_mass.csv", header=TRUE)
leaf$pos_age <- as.factor(paste(leaf$pos, leaf$age, sep = "_"))
normal$pos_age <- as.factor(paste(normal$pos, normal$age, sep = "_"))
normal <- merge(leaf[c("pos_age", "mass_mg")], normal, by = "pos_age")
#check there are no NAs in leaf mass
sum(is.na(normal$mass_mg))
normal$pos_age <- NULL
rownames(normal) <- normal$id
normal <- cbind(normal[,1:12],normal[,13:ncol(normal)]/normal$mass_mg)
rm(leaf)

#remove compounds that are zero across all samples, and then some extraneous columns
dim(normal)
str(normal[,1:13])
#counts the zeros in normal. first get a df that removes all the sample info
normal.trimmed<-normal[,c(13:ncol(normal))]
dim(normal.trimmed)
#i don't know what this line of code does but Marge added it when helping me
normal.trimmed<-data.frame(apply(normal.trimmed,2,as.numeric))
length(colSums(normal.trimmed)==0)
sum(colSums(normal.trimmed)==0)

#says 9
#now removes the zeros in normal
normal.nozero <- normal.trimmed[, which(colSums(normal.trimmed) != 0)]
ncol(normal.trimmed)-ncol(normal.nozero)
#it removed 9
#and let's count zeros in normal.nozero now
sum(colSums(normal.nozero)==0)

normal.new<-data.frame(normal[,1:12],normal.nozero)
rm(normal, normal.trimmed, normal.nozero)

#remove some columns we don't need anymore
normal.new$indiv <- NULL
normal.new$month <- NULL
normal.new$line <- as.factor(paste(normal.new$pop, normal.new$line, sep = "_"))
normal.new$"X4.78_514.2360n" <- NULL
normal.new$mass_mg <-NULL
normal.new$pos <- as.factor(normal.new$pos)
str(normal.new[,1:10])
#now there are 9 columns of information

###############
#we have now normalized by leaf mass, so peaks are "peak abundance/internal standard/mg dry weight of leaf" and we've removed compounds that are zero across all samples

#####
#**************
#add palatability data
normal.new$line_age <- as.factor(paste(normal.new$line, normal.new$age, sep = "_"))
palat = read.csv("line_level_palatability.csv", header = TRUE)
palat$line_age = paste(palat$line, palat$age, sep="_")
palat$X <- NULL
palat$pop <- NULL
palat$line <- NULL
palat$lat <- NULL
palat$age <- NULL
palat_chem <- merge(palat, normal.new, by = "line_age", all = TRUE)
str(palat_chem[,1:14])
#remove anything that has data for palatability but not chemistry. But keep vice versa bc there may be analyses of chemistry that don't need palatability.
palat_chem <- palat_chem[!is.na(palat_chem$sample),]
#dataset has chemical abundance and diversity, and palatability metrics on a per-maternal line level. #biomass is ln(total biomass per cup) (counting dead cats as zero biomass), standardized by intitial cat count and duration
#area is ln(cumulative area consumed per cup) over course of experiment, standardized by intitial cat count and duration
#survival is number of survivors


###########################
#subsetting the big dataframe that has all chemistry and palatability!
###########################
#1:20 is the desirable dilution for the first round of analyses bc it focuses on differentiating abundant compounds.
str(palat_chem[1:14])
#14th column is where compounds start
dilute <- droplevels(subset(palat_chem, palat_chem$dilution == "20"))
dilute.trimmed<-dilute[,14:ncol(dilute)]
dilute.trimmed<-data.frame(apply(dilute.trimmed,2,as.numeric))
sum(colSums(dilute.trimmed)==0)
#now 399 are zero
x <- dilute.trimmed[, which(colSums(dilute.trimmed) != 0)]
sum(colSums(x)==0)

#now relative abundances of x (each peak divided by sum of all peaks for that sample) (test code first)
# y <- x[1:5,1:7]
# rel.y <- y/rowSums(y)
# rel.y2 <- decostand(y, method = "total", margin = 1)

rel.x <- x/rowSums(x)

hist(colSums(rel.x))
#lots of peaks have very low relative abundance across all samples

length(rel.x[, which(colwise(max)(rel.x) >= 0.05)])
#26 compounds are >= 5% relative abundance in at least one sample
length(rel.x[, which(colwise(max)(rel.x) >= 0.01)])
#110 compounds are >= 1% relative abundance in at least one sample
length(rel.x[, which(colwise(max)(rel.x) >= 0.005)])
#207 compounds are >= 0.5% relative abundance in at least one sample
length(rel.x[, which(colwise(max)(rel.x) >= 0.001)])
#793 compounds are >= 0.1% relative abundance in at least one sample

#lists of the abundant compounds
shortlist5 <- colnames(rel.x[, which(colwise(max)(rel.x) >= 0.05)])
shortlist1 <- colnames(rel.x[, which(colwise(max)(rel.x) >= 0.01)])
shortlist.5 <- colnames(rel.x[, which(colwise(max)(rel.x) >= 0.005)])
shortlist.1 <- colnames(rel.x[, which(colwise(max)(rel.x) >= 0.001)])

#RELATIVE
# #dataframe of 1:20, only compounds >= 1% abundance for at least one sample: RELATIVE ABUND
# abundant1 <- data.frame(dilute[,1:13],rel.x[, which(colwise(max)(rel.x) >= 0.01)])
# #dataframe of 1:20, only compounds >= 5% abundance for at least one sample: RELATIVE ABUND
# abundant5 <- data.frame(dilute[,1:13],rel.x[, which(colwise(max)(rel.x) >= 0.05)])
# #dataframe of 1:20, all compounds: RELATIVE ABUND
# all_rel <- data.frame(dilute[,1:13],rel.x)

#NOT RELATIVE
# #dataframe of 1:20, only compounds >= 1% abundance for at least one sample: NOT RELATIVE
abundant1 <- data.frame(dilute[,1:13],dilute.trimmed[, colnames(dilute.trimmed)%in%shortlist1])
abundant.5 <- data.frame(dilute[,1:13],dilute.trimmed[, colnames(dilute.trimmed)%in%shortlist.5])
abundant.1 <- data.frame(dilute[,1:13],dilute.trimmed[, colnames(dilute.trimmed)%in%shortlist.1])
# #dataframe of 1:20, only compounds >= 5% abundance for at least one sample: NOT RELATIVE
#abundant5 <- data.frame(dilute[,1:13],dilute.trimmed[, colnames(dilute.trimmed)%in%shortlist5])
#rankedabund <- data.frame(sort(colSums(abundant5[,14:ncol(abundant5)])))
#write.csv(abundant1, file = "20180507_chem_abundant1percent.csv")
write.csv(abundant1, file = "20180531_chem_abundant1percent.csv")
write.csv(abundant.5, file = "20180531_chem_abundant.5percent.csv")
write.csv(abundant.1, file = "20180531_chem_abundant.1percent.csv")
# write.csv(abundant5, file = "20180424_chem_abundant5percent.csv")
# write.csv(all_rel, file = "20180424_chem_all_dilute.csv")

# par(mfrow=c(3, 3))
# for(i in 14:nrow(abundant5)){
#   plot(abundant5[,i] ~ abundant5$lat, ylab = colnames(abundant5[i]), col = abundant5$age)
# }
# par(mfrow=c(1,1))
#young are red

rm(x, rel.x, dilute.trimmed, dilute, normal.new, palat, palat_chem, shortlist1, shortlist5)

#eventually maybe we want to calculate relative mass defect (see lab notebook p 47) of the compounds that emerge as interesting--?