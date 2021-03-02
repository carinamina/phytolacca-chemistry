#started 4-23-18 during jury duty
#combining all defensive traits into one df for random forest analysis
#line level means
#ID: line_age (can have other characteristics and then remove into a key at the end?)
#leaf age
#palatability
#C:N
#C
#N
#toughness
#chemical richness
#chemical diversity
#10,000 chemical compounds

source("/Users/carina/Documents/dissertation/Phytolacca writing/ch3 defensive traits/ch3 defensive traits analysis/cn_wrangling.R")

##########start with the relatively clean C:N data
#make a new variable that is conversion rate of biomass per leaf area consumed
palat_cn$conv <- palat_cn$biomass/palat_cn$area
all <- palat_cn[c("line_age","pop","line","age","lat","region","conv","C_N","percent_C","percent_N")]
all$region <- revalue(all$region, c("southern" = "subtropical","northern" = "temperate"))
rm(cn, palat_cn,lme_results,lme_slopes)

################add leaf toughness
source("/Users/carina/Documents/dissertation/Phytolacca writing/ch3 defensive traits/ch3 defensive traits analysis/toughness_wrangling.R")

all <- merge(all, tough_line, by = "line_age", all = TRUE)
rm(tough_line)

################
#add chemistry data. Now this is real (not relative) abundance of all compounds that were > 1% relative abundance in at least one sample.
library(vegan)
chem <- read.csv("20180507_chem_abundant1percent.csv", header = TRUE)
chem$X <- NULL
str(chem[,1:14])

chem.mat <- data.matrix(chem[,14:length(chem)])
chem_sums <- chem[c("line_age")]
chem_sums$shann <- diversity(chem.mat)  #Shannon index
hist(chem_sums$shann)
#counts the number of zeroes in each row and subtracts from total to calculate richness
chem_sums$richness <- ncol(chem.mat) - rowSums(chem.mat == 0)
hist(chem_sums$richness)
chem_sums$abund <- rowSums(chem.mat)
hist(chem_sums$abund)
# #test code
# test <- chem.mat[1:17,1:4]
# rowSums(test== 0)

all <- merge(all, chem_sums, by = "line_age", all = TRUE)
all <- merge(all, chem[c(1,14:ncol(chem))], by = "line_age", all = TRUE)
#can't use na.omit because region has NA for many populations
#removes rows that have chemistry but not other traits
all <- all[!is.na(all$pop),]
#removes rows that have other traits but not chemistry
all <- all[!is.na(all$shann),]
all <- all[!is.na(all$tough),]
all <- all[!is.na(all$conv),]

key <- all[,1:6]
str(all[,1:15])
rm(chem, chem_sums, chem.mat)

#for future analyses, we'd want to keep population and make lists of the four different analyses to run: population only, chemistry peaks only, other traits, all traits. The latter three do not have population. But for now let's just do all traits, without population.

#will write tab-delimited datasets to my HPCC directory
setwd("/Volumes/baskettc/Documents/")

#lists of populations only, chemsummary, chemicals only, non-chemical traits only
#nochem: C_N, C, N, toughness
write.table(colnames(all[,8:11]), "nochem_list.txt", sep="\t", row.names = FALSE, quote=FALSE, col.names=FALSE)
#chemsum: shann, richness, abundance
write.table(colnames(all[,12:14]), "chemsum_list.txt", sep="\t", row.names = FALSE, quote=FALSE, col.names=FALSE)
#chemical compounds
write.table(colnames(all[,15:ncol(all)]), "chem_list.txt", sep="\t", row.names = FALSE, quote=FALSE, col.names=FALSE)
#export line_age and pop to turn into a binary matrix
write.table(all[,1:2], "all_pop", sep="\t", row.names = FALSE, quote=FALSE)
pop_binary <- read.csv("all_pop_binary.matrix.txt", sep = "\t", header=TRUE)
#write the population binary list
write.table(colnames(pop_binary[,2:ncol(pop_binary)]), "pop_list.txt", sep="\t", row.names = FALSE, quote=FALSE, col.names=FALSE)
#re-merge populations with all, now as a binary factor
all <- merge(all, pop_binary, by = "line_age")

mature <- subset(all, all$age == "mature")
mature <- mature[c(1,7:ncol(mature))]
write.table(mature, "mature_20180507", sep="\t", row.names = FALSE)

young <- subset(all, all$age == "young")
young <- young[c(1,7:ncol(young))]
write.table(young, "young_20180507", sep="\t", row.names = FALSE)

setwd("/Users/carina/Documents/R_working_directory")
