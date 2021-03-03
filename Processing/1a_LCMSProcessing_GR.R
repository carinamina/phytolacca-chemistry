#chemistry data wrangling

#notes: see metadata for Telmisartan mass, which is different than last time (???); should be same abundance across all samples (check blanks)
#how to use blanks? Remove compounds that are highest in blanks.
#toss out all the "induced" data
#then remove compounds that are zero in everything (??) Stuff that was only in induced samples.
#how else can we have an abundance cutoff??
library(tidyverse)

raw = read.csv("Raw/20180727_LCMS_ReadyForR.csv", header=TRUE)
#for some reason every other column is blank. There shouldn't be NAs otherwise; indeed this step just removes half the rows
raw <- na.omit(raw)

#change the first column (compounds) to row names and then remove it
rownames(raw) <- raw[,1]
raw <- raw[,2:length(raw)]

#transpose to columns as compounds and samples as rows
raw <- as.data.frame(t(raw))

sample.df <- as.data.frame(rownames(raw))
sample.df$LCMS_order <- as.numeric(substr(sample.df[,1], 12, 14))

key <- read.csv("Raw/20180727_LCMS_SampleKey.csv", header=TRUE)
key <- left_join(sample.df, key, by="LCMS_order")
rm(sample.df)

#some of the bulk lines have the same number as non-bulk lines from those populations, so append a B to make them unique
key$line <- ifelse(grepl("Bulk",key$pop, fixed=TRUE) | grepl("bulk",key$pop, fixed=TRUE), paste(key$line, "B"),key$line)
#Remove the bulk label so we know which population is which. Would be good to check someday these are no different
key$pop <- plyr::revalue(key$pop, c("Val" = "Tiri","Tiri Bulk" = "Tiri", "Hwy 27 bulk" = "Hwy 27", "TT Bulk" = "TT", "Dalton Bulk" = "Dalton", "Gav Bulk" = "Gav", "Whitehall Bulk" = "Whitehall","KBS Bulk" = "KBS")) 
str(key)
#get latitude for each population
lat = read.csv("Raw/LatsPopsKey.csv", header = TRUE)
key <- left_join(key, lat, by = "pop")
#get regional name for each population
regions <- c("tropical", "tropical","tropical", "subtropical", "subtropical", "subtropical",NA,NA,NA,NA,NA,NA,NA, "temperate", "temperate", "temperate")
regions <- cbind(regions, sort(unique(key$lat)))
colnames(regions) <- c("region", "lat")










#stopped here...this only works if merging, but then the blanks are lost.
abc <- left_join(key, regions, by = "lat")
rm(regions,lat)

con_key <- droplevels(subset(key, key$defense == "con"))

#replace rownames of chem data with a useful id: pos_age_dilution which matches something in master key
#key$id <- as.factor(paste(key$pos, key$age, key$dilution, sep = "_"))
chem.df$sample <- rownames(chem.df)
chem.df <- left_join(con_key, chem.df, by = "sample")
#rownames(chem.df) <- chem.df$id



###############
#now we are done with identifying samples and their categorical characteristics. We need to normalize by telmisartin

#check for low telmisartin peaks
# hist(chem.df$"4.78_514.2360n") #does not exist !!! ???
# chem.df[which(chem.df$"4.78_514.2360n"<20000),1:11] #does not exist !!! ???
#33_mature_20 is lower than normal by an order of magnitude. Tony calculated the abundance according to MassLynx of this sample and the 33_young_20, and we used that ratio to determine that the telmisartan peak should be 77,900 for this sample. (We checked another peak and it seemed to have aligned normally, so it's not clear why Progenesis messed this one up but it doesn't seem to be a problem with the entire sample)
# chem.df[which(chem.df$id == "33_mature_20"),]$"4.78_514.2360n" = 77900 #does not exist !!! ???

#normalize by telmisartin: divide all peaks in each sample by telmisartan of that sample. It's 4.78_514.2360n
# normal <- cbind(chem.df[,1:11],  chem.df[,12:ncol(chem.df)]/chem.df$"4.78_514.2360n") #does not exist !!! ???
# normal$"4.78_514.2360n" #does not exist !!! ???

normal <- chem.df

#look at distributions of total abundance
abund <- rowSums(normal[,12:nrow(normal)])
totalabund <- data.frame(normal[,1:11], abund)
hist(totalabund$abund)
#not4jul hist(subset(totalabund, totalabund$dilution == "20")$abund)
#not4jul hist(subset(totalabund, totalabund$dilution == "1")$abund)
#not4jul summary(subset(totalabund, totalabund$dilution == "20")$abund)
#not4jul summary(subset(totalabund, totalabund$dilution == "1")$abund)
#these seem ok this time, yay!

# reshape totalabund to wide, with a column for each dilution, and check that ratio of undiluted is ~~20x diluted.
# Should catch any really weird ratios
totalabund$pos_age <- as.factor(paste(totalabund$pos, totalabund$age, sep = "_"))
#not4jul ratio <- reshape(totalabund,  timevar ="dilution", idvar = "pos_age", direction = "wide")
#not4jul ratio$ratio <- ratio$abund.1/ratio$abund.20
#not4jul hist(ratio$ratio)
#there is one with a pretty low ratio, 76_young. It seems like the problem is with the total abundance of
# the nondiluted sample; it is the lowest of all samples, even though the chromatograms look normal.
# So I think it's ok because we don't care about nondiluted.
#not4jul rm(totalabund,ratio,abund,chem.df)


################
#to review, at this point we have removed all the membrane lipids (>7:30 retention),
#divided all samples by their own internal standard peak (normalized by telmisartin)
#we need to normalize by leaf mass and remove compounds that are zero for all samples.

#normalize by leaf mass
leaf <- read.csv("data_in/2017_chem_leaf_mass.csv", header=TRUE)
leaf$pos_age <- as.factor(paste(leaf$pos, toupper(substr(leaf$age, 1, 1)), sep = "_"))
normal$pos_age <- as.factor(paste(normal$id, normal$age, sep = "_"))
normal <- merge(leaf[c("pos_age", "mass_mg")], normal, by = "pos_age")
#check there are no NAs in leaf mass
sum(is.na(normal$mass_mg))
#normal$pos_age <- NULL
#rownames(normal) <- normal$id
normal <- cbind(normal[,1:13],normal[,14:ncol(normal)]/normal$mass_mg)
rm(leaf)

#remove compounds that are zero across all samples, and then some extraneous columns
dim(normal)
str(normal[,1:14])
#counts the zeros in normal. first get a df that removes all the sample info
normal.trimmed<-normal[,c(14:ncol(normal))]
dim(normal.trimmed)
#i don't know what this line of code does but Marge added it when helping me
#gri: looks like it makes sure everything is numeric... 
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
str(normal.new[,1:16])
#now there are 9 columns of information

###############
#we have now normalized by leaf mass, so peaks are "peak abundance/internal standard/mg dry weight of leaf"
#and we've removed compounds that are zero across all samples

#####
#**************
#add palatability data
normal.new$line_age <- as.factor(paste(normal.new$line, normal.new$age, sep = "_"))
palat = read.csv("data_in/line_level_palatability.csv", header = TRUE)
palat$line_age = paste(palat$line, palat$age, sep="_")
palat$X <- NULL
palat$pop <- NULL
palat$line <- NULL
palat$lat <- NULL
palat$age <- NULL
palat_chem <- merge(palat, normal.new, by = "line_age", all = TRUE)
str(palat_chem[,1:17])
#remove anything that has data for palatability but not chemistry.
#But keep vice versa bc there may be analyses of chemistry that don't need palatability.
palat_chem <- palat_chem[!is.na(palat_chem$sample),]
#dataset has chemical abundance and diversity, and palatability metrics on a per-maternal line level.
#biomass is ln(total biomass per cup) (counting dead cats as zero biomass), standardized by
# intitial cat count and duration
#area is ln(cumulative area consumed per cup) over course of experiment, standardized by
# intitial cat count and duration
#survival is number of survivors


###########################
#subsetting the big dataframe that has all chemistry and palatability!
###########################
#1:20 is the desirable dilution for the first round of analyses bc it focuses on differentiating abundant compounds.
str(palat_chem[1:14])
#14th column is where compounds start
#not4jul dilute <- droplevels(subset(palat_chem, palat_chem$dilution == "20"))
#not4jul dilute.trimmed<-dilute[,14:ncol(dilute)]
#not4jul dilute.trimmed<-data.frame(apply(dilute.trimmed,2,as.numeric))
#not4jul sum(colSums(dilute.trimmed)==0)
#now 399 are zero
#not4jul x <- dilute.trimmed[, which(colSums(dilute.trimmed) != 0)]
#not4jul sum(colSums(x)==0)

#now relative abundances of x (each peak divided by sum of all peaks for that sample) (test code first)
# y <- x[1:5,1:7]
# rel.y <- y/rowSums(y)
# rel.y2 <- decostand(y, method = "total", margin = 1)

# dirty !!! ???
palat_chem.trimmed <- palat_chem[,16:ncol(palat_chem)]
palat_chem.trimmed<-data.frame(apply(palat_chem.trimmed,2,as.numeric))
x <- palat_chem.trimmed[, which(colSums(palat_chem.trimmed) != 0)]
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
abundant1 <- data.frame(palat_chem[,1:16],palat_chem.trimmed[, colnames(palat_chem.trimmed)%in%shortlist1])
abundant.5 <- data.frame(palat_chem[,1:16],palat_chem.trimmed[, colnames(palat_chem.trimmed)%in%shortlist.5])
abundant.1 <- data.frame(palat_chem[,1:16],palat_chem.trimmed[, colnames(palat_chem.trimmed)%in%shortlist.1])
# #dataframe of 1:20, only compounds >= 5% abundance for at least one sample: NOT RELATIVE
#abundant5 <- data.frame(dilute[,1:13],dilute.trimmed[, colnames(dilute.trimmed)%in%shortlist5])
#rankedabund <- data.frame(sort(colSums(abundant5[,14:ncol(abundant5)])))
#write.csv(abundant1, file = "20180507_chem_abundant1percent.csv")
write.csv(abundant1, file = "data_out/20180531_chem_abundant1percent.csv")
write.csv(abundant.5, file = "data_out/20180531_chem_abundant.5percent.csv")
write.csv(abundant.1, file = "data_out/20180531_chem_abundant.1percent.csv")
# write.csv(abundant5, file = "20180424_chem_abundant5percent.csv")
# write.csv(all_rel, file = "20180424_chem_all_dilute.csv")

# par(mfrow=c(3, 3))
# for(i in 14:nrow(abundant5)){
#   plot(abundant5[,i] ~ abundant5$lat, ylab = colnames(abundant5[i]), col = abundant5$age)
# }
# par(mfrow=c(1,1))
#young are red

rm(x, rel.x, normal.new, palat, palat_chem, shortlist1, shortlist5)

#eventually maybe we want to calculate relative mass defect (see lab notebook p 47) of the compounds that emerge as interesting--?
