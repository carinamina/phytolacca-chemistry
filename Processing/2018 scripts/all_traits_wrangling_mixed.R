#started 5-10-18 
#combining all defensive traits into one df for mixed model analysis
#difference from RF wrangling: here it's ok to have NAs for some traits, as they are analyzed separately
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
setwd("/Users/carina/Documents/R_working_directory")
source("/Users/carina/Documents/dissertation/Phytolacca writing/ch3 defensive traits/ch3 defensive traits analysis/cn_wrangling.R")

##########start with the relatively clean C:N data
#make a new variable that is conversion rate of biomass per leaf area consumed
palat_cn$conv <- ifelse(palat_cn$area > 0, palat_cn$biomass/palat_cn$area,0)
all <- palat_cn[c("line_age","conv","C_N","percent_C","percent_N")]
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
str(all)
rm(chem, chem_sums, chem.mat)

#add NMDS axes for each leaf age after RF analysis
#would like to just source the code but it was long enough that I used R markdown and that's weird to source, so I'll export those and import to here.
nmds_scores <- read.csv("20180531_nmds_scores.csv", header = TRUE)
nmds_scores$lat <- NULL
nmds_scores$region <- NULL
nmds_scores$age <- NULL
all <- merge(all, nmds_scores, by = "line_age", all = TRUE)
rm(nmds_scores)
str(all)

#get the population, line, leaf age, and region in a key
all_key <- read.table(text = all$line_age, sep = "_", colClasses = "factor")
colnames(all_key) <- c("pop","line","age")
all_key$line <- as.factor(paste(all_key$pop, all_key$line, sep="_"))
all_key$line_age <- as.factor(paste(all_key$line, all_key$age, sep = "_"))
lats = read.csv("lats_long_names.csv", header=TRUE)
all_key <- merge(all_key, lats, by = "pop")
regions <- c("tropical", "tropical","tropical", "subtropical", "subtropical", "subtropical",NA,NA,NA,NA,NA,NA,NA, "temperate", "temperate", "temperate")
regions <- cbind(regions, sort(unique(all_key$lat)))
colnames(regions) <- c("region", "lat")
all_key <- merge(all_key, regions, by = "lat")

all <- merge(all_key, all, by = "line_age")
rm(all_key, lats, regions)
