#LC/MS data wrangling
library(tidyverse)

raw = read.csv("Raw/20180727_LCMS_ReadyForR.csv", header=TRUE)
#for some reason every other row is blank. There shouldn't be NAs otherwise; indeed this step just removes half the rows
raw <- na.omit(raw)

#change the first column (compounds) to row names and then remove it
rownames(raw) <- raw[,1]
raw <- raw[,2:length(raw)]
#transpose to columns as compounds and samples as rows
raw <- as.data.frame(t(raw))

##################
# WRANGLING THE SAMPLE KEY
##################
# Create the LCMS_ID, which matches row names in chemistry dataset
# Add latitude & region 
# Add leaf mass data from another dataset

#gets rownames (what I call LCMS_ID) from chemistry dataset and extracts the LCMS_order to match up with the key
sample.df <- as.data.frame(rownames(raw)) %>% dplyr::rename(LCMS_ID = "rownames(raw)" )
sample.df$LCMS_order <- as.numeric(substr(sample.df$LCMS_ID, 12, 14))
#read in the sample key, which has information based on the LCMS_order. I remove sample_order because it's confusing to have two different order numbers (sample_order is only the leaf samples of interest, while LCMS_order was based on the actual vials, including pooled and blanks and fruits etc)
key <- read.csv("Raw/20180727_LCMS_SampleKey.csv", header=TRUE) %>% select(-sample_order) %>% dplyr::rename(pos = id)
key <- left_join(sample.df, key, by="LCMS_order") %>% select(-c(LCMS_order)) 
rm(sample.df)
#break off the blanks because they don't survive upcoming mergers and subsetting
blanks <- key %>% filter(label == "blank")
key <- key %>% filter(label != "blank")

#add data on leaf mass
leaf <- read.csv("Raw/2017_LCMS_LeafMass_RAW.csv", header=TRUE)
leaf$pos_age <- as.factor(paste(leaf$pos, toupper(substr(leaf$age, 1, 1)), sep = "_"))
key$pos_age <- as.factor(paste(key$pos, key$age, sep = "_"))
key <- left_join(key, leaf[c("pos_age", "mass_mg")], by = "pos_age")
#check there are no NAs in leaf mass
sum(is.na(key$mass_mg))
rm(leaf)

#use only constitutive defense samples. Count how many samples per population
con_key <- key %>% filter(defense == "con") %>% select(-defense)
con_key %>% count(pop)

# cleaning up population names, adding latitude and region
#some of the bulk lines have the same number as non-bulk lines from those populations, so I will remove them as pseudoreplicates. since not every line is represented in the constitutive datasheet, this might not be a problem. These are the ones from the toughness dataset "Tiri 12", "KBS 8","TT 5","Whitehall 2","Gav 1". Only Gav and TT have both the bulk and non-bulk
sort(unique(paste(con_key$pop,con_key$line,sep="_")))
con_key <- filter(con_key, !pop %in% c("TT Bulk","Gav Bulk"))
#combine Valverde with Tiri (nearby) and remove bulk label from three others. I am a little nervous about these bulked plants because I don't remember exactly how it was done--did we just use fruits that we found on the plants? How did we know who the father was? But at worst, if we used fruits without knowing whether they were crossed by stray wasps, we would be eroding population differentiation and adding noise to the data.
con_key$pop <- plyr::revalue(con_key$pop, c("Val" = "Tiri","Tiri Bulk" = "Tiri", "Hwy 27 bulk" = "Hwy 27", "Dalton Bulk" = "Dalton"))
str(con_key)
#get latitude for each population
lat = read.csv("Raw/LatsPopsKey.csv", header = TRUE)
con_key <- left_join(con_key, lat, by = "pop")
rm(lat)

#add blanks back
con_key <- add_row(con_key, LCMS_ID = blanks[,1], label = blanks[,2])
rm(key,blanks)
###############
# subset chemistry to constitutive samples; remove internal standard; remove compounds with zero abundance across all samples; perform basic checks

con <- raw[which(rownames(raw) %in% con_key$LCMS_ID == TRUE),]

#check that telmisartan (see metadata: 4.75_514.2383n) is equivalent across all samples (should already be normalized by it), and then remove that peak
length(unique(con$"4.75_514.2383n")) == 1

con <- con %>% select(-"4.75_514.2383n")

#look at distributions of total abundance
#for each sample: many are quite low, I assume the blanks are the ones even lower
hist(rowSums(con))
#for each compound: the vast majority are so tiny that you can hardly see the tail. a double square-root is easier to see
hist(sqrt(sqrt(colSums(con))))
hist(log(colSums(con)+1))
#the log histogram is very normal...so most compounds are intermediate abundance

#remove compounds with 0 abundance across all samples
con <- con[,colSums(con) > 0]
#63 compounds. Not much of a reduction, but I'll take it

#remove compounds only found in 1-2 samples, which can't possibly be useful in a model
countpresence <- apply(con,2,function(x) sum(x > 0))
length(countpresence[countpresence < 3])
#that's another 48...it goes up by about 25 per number as threshold is raised to 6
names(countpresence[countpresence < 3])
con <- con %>% select(names(countpresence[countpresence > 2]))
rm(countpresence)

############################
#examine some properties of compound abundances and blanks to try and figure out criteria to reduce the list
#e.g. remove compounds with their max abundance in the blanks

#how many compounds are >0 in the blanks?
sum(con[nrow(con)-1,] > 0)
sum(con[nrow(con),] > 0)
#similar, >500

#histograms for the two blanks vs max abund for each compound
hist(as.numeric(con[nrow(con),]))
hist(as.numeric(con[nrow(con)-1,]))
maxabund <- apply(con,2,max)
hist(maxabund)
#max in the blanks is an order of magnitude lower than max in the samples

#make a dataframe with compounds as rows to see presence counts, max abundance, summed abundance
sumabund <- apply(con,2,sum)
countnonzero <- apply(con,2,function(x) sum(x > 0)) #this is the same as countpresence but now on a slightly smaller dataset so i named it differently
summary(countnonzero)
hist(countnonzero)
#a significant fraction are present in <1/3 of samples, and a lot are found in ~all samples
abund <- as.data.frame(cbind(maxabund,sumabund,countnonzero))
abund$compound <- rownames(abund)
rm(maxabund,sumabund,countnonzero)

#find the sample that has the maximum abundance of each compound
abund$LCMS_ID_max <- NA
for(i in 1:nrow(abund)){
  abund$LCMS_ID_max[i] <- rownames(con[con[,i] == abund$maxabund[i],])
}
rm(i)
#SHAME! SHAME! A LOOP! I gave up on finding a non-loop solution, it was hard.

abund <- left_join(abund, con_key[c("LCMS_ID","age","pop","region","lat")], by = c("LCMS_ID_max" = "LCMS_ID"))

#how many compounds are highest in the blanks?
sum(abund$LCMS_ID_max == "XS2_072718_001" | abund$LCMS_ID_max == "XS2_072718_321")
#144! That's something. List them out:
(highestblanks <- abund[which(abund$LCMS_ID_max == "XS2_072718_001" | abund$LCMS_ID_max == "XS2_072718_321"),c("compound")])
#interesting, many are quite high-retention time. Do I remember correctly that these are lipids and they can gradually gunk up the column?

#if there's anything else we could do with the blanks, it would be here, but I'm not sure. Simply subtract their abundances from everything? That seems kind of like a normalization step which we already did with an internal standard. 
#What about using the average abundance of peaks in the blanks (after contaminant removal) to set a threshold of "above this is good stuff, below is noise"?
#this step removes compounds with highest peaks in the blanks, and then removes the blanks from the chemistry dataset since we don't need that information anymore
LCMS_ID_blanks <- con_key[con_key$label == "blank",]$LCMS_ID
con <- con %>% select(-all_of(highestblanks)) %>% filter(rownames(con) %in% LCMS_ID_blanks == F)

#################################
#a few plots of max and sum abund against latitude and presence count
#also interesting to compare to the abundance dataset without the compounds highest in the blanks
abund2 <- abund %>% filter(compound %in% highestblanks == F)

#max abundance vs latitude
plot(log(abund$maxabund) ~ jitter(abund$lat))
#hard to see anything useful in this plot
plot(log(abund$maxabund) ~ abund$countnonzero)
plot(log(abund2$maxabund) ~ abund2$countnonzero)
#that's interesting, maxabund can be very low and the compound is still found in a large proportion of samples. And vice versa, maxabund can be quite high and the compound is only found in a handful of samples
plot(log(abund$sumabund) ~ abund$countnonzero)
plot(log(abund2$sumabund) ~ abund2$countnonzero)
#perhaps more as expected, compounds with high summed abundance tend to be found in more samples
#comparing the latter two plots with and without the compounds highest in the blanks, there are a surpising number of dots that disappear that were in a lot of samples. I still think it's a good criteria!

rm(highestblanks,abund,abund2,LCMS_ID_blanks)

################
#normalize peaks by leaf mass
#create new df called normal that divides peak areas by leaf mass: peak area/mg dry weight of leaf
normal <- cbind(rownames(con),con)
normal <- right_join(select(con_key, c(LCMS_ID,mass_mg)), normal, by = c("LCMS_ID" = "rownames(con)"))
normal <- cbind(normal[,1:2],normal[,3:ncol(normal)]/normal$mass_mg) %>% select(-mass_mg)

# some checks on the normalization step
# abc <- cbind(normal[,1:2],normal[,3:ncol(normal)]/normal$mass_mg)
# normal[,6]/normal$mass_mg == abc[,6]
# normal[,670]/normal$mass_mg == abc[,670]

##################################
#join the chemistry data with the key and export
con_key$age <- plyr::revalue(con_key$age, c("Y" = "young", "M" =  "mature"))
full <- inner_join(con_key[c("LCMS_ID","pos","pop","line","age","lat","region")], normal, by = "LCMS_ID") %>% select(-LCMS_ID) %>% mutate(pos_age = paste(pos,age,sep="_"),.before=pos) %>% mutate(line_age=paste(paste(pop,line,sep="_"),age,sep="_"),.before=pos)

write.csv(full, "Processing/1a_out_LCMS_Indiv.csv",row.names=F)

rm(list=ls())
