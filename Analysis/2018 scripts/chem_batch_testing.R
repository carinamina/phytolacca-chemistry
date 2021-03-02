#start jul 24, 2018
#testing whether there is a batch effect on chemistry (from Dec and Feb batches), and if so whether it goes away when using an absolute rather than relative abundance cutoff

#problem: some figures in chapter 3 that show the whole PHAM gradient for chemical traits (in particular NMDS1, richness, and abundance) show a curious spatial pattern: the middle seven populations look different than the three northern and three southernmost. This is alarming because two separate LC-MS batches were run (I didn't realize batch effects were possible): 3/5 of maternal lines were run in Dec for the three southern and northernmost PHAM (and for PHRI). The remaining 2/5 of these 9 populations were run in February, along with all five lines of the seven "intermediate" PHAM populations. I'm used to the 3 FL populations showing unusual patterns compared to the rest of PHAM, but it's not normal that the 3 northern populations would too, so that means maybe there was a batch effect. 
#It's possible said effect is contingent on the method of reducing the dataset to a manageable number of compounds. I used a relative abundance cutoff: data is included for any compound that is found at >1% abundance in at least one sample. I've tried .5% and .1% as well, and it doesn't drastically change results. But the relative abundance could be a problem if there are more RARE compounds in tropical plants (which is probably true), thus "penalizing" tropical plants more by cutting out more compounds, which would affect measures like abundance, richness, diversity.

#GOALS:
#keep some batch or date tag in dataset
#reduce the vast chem dataset to the 9 populations that had runs in both batches
#cut out compounds that are all zero
#test whether there is a batch*age*region effect (hopefully there's enough data) on different subsets of data: all of it, some relative cutoffs, some absolute cutoffs. Focus on richness and abundance before getting into diversity and NMDS (which are contingent on not messing up abundance and richness)
#be sure to compare visualizations of these things with and without cutoffs to understand how cutoffs affect them

#basis: chem_data_wrangling
#blanks are not subtracted
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
#rm(raw)

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

#before normalizing...how many peaks have abundance of >10,000? >1,000?
chem.df.trim <- droplevels(subset(chem.df, is.na(chem.df$region)==FALSE & chem.df$dilution == 20))[,12:ncol(chem.df)]
length(chem.df.trim[, which(colwise(max)(chem.df.trim) >= 1000)])
#1250
length(chem.df.trim[, which(colwise(max)(chem.df.trim) >= 10000)])
#170
shortlist.abs <- colnames(chem.df.trim[, which(colwise(max)(chem.df.trim) >= 10000)])
#oddly, later the compound names gain an X prefix...need to add it to these to be usable for later subsetting
shortlist.abs <- paste("X",shortlist.abs,sep="")

#normalize by telmisartin: divide all peaks in each sample by telmisartan of that sample. It's 4.78_514.2360n
normal <- cbind(chem.df[,1:11],  chem.df[,12:ncol(chem.df)]/chem.df$"4.78_514.2360n")
normal$"4.78_514.2360n"

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
#counts the zeros in normal. first reduce to only those populations that were in the Dec sample, and only 1:20. then get a df that removes all the sample info
normal<- droplevels(subset(normal, is.na(normal$region)==FALSE & normal$dilution == 20))
normal.trimmed<-normal[,c(13:ncol(normal))]
dim(normal.trimmed)
#i don't know what this line of code does but Marge added it when helping me
normal.trimmed<-data.frame(apply(normal.trimmed,2,as.numeric))
length(colSums(normal.trimmed)==0)
sum(colSums(normal.trimmed)==0)

#says 640 are all zero
#now removes the zeros in normal
normal.nozero <- normal.trimmed[, which(colSums(normal.trimmed) != 0)]
ncol(normal.trimmed)-ncol(normal.nozero)
#it removed 640
#and let's count zeros in normal.nozero now
sum(colSums(normal.nozero)==0)

normal.new<-data.frame(normal[,1:12],normal.nozero)
rm(normal, normal.trimmed, normal.nozero)

#remove some columns we don't need anymore
normal.new$indiv <- NULL
normal.new$dilution <- NULL
normal.new$line <- as.factor(paste(normal.new$pop, normal.new$line, sep = "_"))
normal.new$"X4.78_514.2360n" <- NULL
normal.new$mass_mg <-NULL
normal.new$pos <- as.factor(normal.new$pos)
str(normal.new[,1:10])
#now there are 9 columns of information

all <- normal.new
rm(normal.new)

###############
#we have now normalized by leaf mass, so peaks are "peak abundance/internal standard/mg dry weight of leaf," we've removed compounds that are zero across all samples, and we've reduced to only the three regions that should have samples in both batches

str(all[1:10])
#10th column is where compounds start
all.trimmed <- all[10:ncol(all)]

#RELATIVE ABUNDANCE CUTOFF
rel <- all.trimmed/rowSums(all.trimmed)
hist(colSums(rel))
#lots of peaks have very low relative abundance across all samples
length(rel[, which(colwise(max)(rel) >= 0.01)])
#108 compounds are >= 1% relative abundance in at least one sample
#list of the relatively abundant compounds
shortlist.rel <- colnames(rel[, which(colwise(max)(rel) >= 0.01)])
#dataframe of only compounds >= 1% abundance for at least one sample
abund.rel <- data.frame(all[,1:9],all.trimmed[, colnames(all.trimmed)%in%shortlist.rel])
abund.rel.mat <- data.matrix(abund.rel[,10:ncol(abund.rel)])
rm(rel, shortlist.rel)

#ABSOLUTE ABUNDANCE CUTOFF isn't working
#dataframe of only compounds >= 10000 "counts" for at least one sample
#see shortlist created before normalization
#abund.abs <- data.frame(all[,1:9],all.trimmed[, colnames(all.trimmed)%in%shortlist.abs])

#df with sample info then richness and abundance from each of the approaches
chem_sums <- all[,1:9]

chem_sums$rich.all <- ncol(all.trimmed) - rowSums(all.trimmed == 0)
hist(chem_sums$rich.all)
range(chem_sums$rich.all)
chem_sums$abund.all <- rowSums(all.trimmed)
hist(chem_sums$abund.all)
range(chem_sums$abund.all)

chem_sums$rich.rel <- ncol(abund.rel.mat) - rowSums(abund.rel.mat == 0)
hist(chem_sums$rich.rel)
range(chem_sums$rich.rel)
chem_sums$abund.rel <- rowSums(abund.rel.mat)
hist(chem_sums$abund.rel)
range(chem_sums$abund.rel)

#richness: ranges over a much greater magnitude using all vs. relatively abundant.
#abudance: min is not that different but max abundance is almost double using all vs. relatively abundant

boxplot(chem_sums$abund.all ~ chem_sums$region)
boxplot(chem_sums$abund.rel ~ chem_sums$region)
#but regional relationships of abundance don't change

boxplot(chem_sums$rich.all ~ chem_sums$region)
boxplot(chem_sums$rich.rel ~ chem_sums$region)
#richness might be greater for tropical when you use all compounds

boxplot(chem_sums$abund.all ~ chem_sums$month)
boxplot(chem_sums$abund.rel ~ chem_sums$month)
#well Feb has more compounds. great.

boxplot(chem_sums$rich.all ~ chem_sums$month)
boxplot(chem_sums$rich.rel ~ chem_sums$month)
#and more types of compounds. great.

library(nlme)
library(corrplot)
library(Hmisc)
library(lsmeans)

chem_sums$pop_age <- as.factor(paste(chem_sums$pop, chem_sums$age, sep="_"))

lme.abund.all <- lme(abund.all ~ age*region*month, random=~1|pop, method = "REML", data=chem_sums, weights=varIdent(form=~1|pop_age))
lme_results(lme.abund.all)
#oh yeah. month makes a big difference for abundance of all compounds.

lme.abund.rel <- lme(abund.rel ~ age*region*month, random=~1|pop, method = "REML", data=chem_sums, weights=varIdent(form=~1|pop_age))
lme_results(lme.abund.rel)
#and for abundance of relatively high abundance compounds

lme.rich.all <- lme(rich.all ~ age*region*month, random=~1|pop, method = "REML", data=chem_sums, weights=varIdent(form=~1|pop_age))
lme_results(lme.rich.all)
#and for richness of all compounds.

lme.rich.rel <- lme(rich.rel ~ age*region*month, random=~1|pop, method = "REML", data=chem_sums, weights=varIdent(form=~1|pop_age))
lme_results(lme.rich.rel)
#and for richness of relatively high abundance compounds


#########
#to conclude, there is definitely a batch effect. I have no idea whether this was actually the machine, or some glitch along the way with processing, there were so many steps where something could have gone wrong.
#i still don't know whether there could be a way to eliminate it with an "absolute cutoff"...but considering that the "relative cutoff" approach finds similar patterns between regions and months to no cutoff, I'm not sure why an absolute cutoff would help...?