#output will be the following traits merged by line_age, with information on latitude and region:
# palatability (script 1e)
# chemical abundance, richness, diversity (1b)
# toughness (1c)
# carbon, nitrogen, C:N (1d)
# ~2000 chemical abundance peaks (1a)
library(tidyverse)

#chemsums and compounds have latitude etc, the others all just line_age. I think the only way to guarantee no NA in pop,line,age is to extract it from line_age once all the data is combined

palat <- read.csv("Processing/1e_out_Palat_Line.csv", header=T)
chemsums <- read.csv("Processing/1b_out_ChemSummaries_Indiv.csv",header=T) %>% select(-c(pos_age,pos,pop,line,age,lat,region))
tough <- read.csv("Processing/1c_out_Toughness_Line.csv",header=T)
cn <- read.csv("Processing/1d_out_CarbonNitrogen_Line.csv",header=T)
compounds <- read.csv("Processing/1a_out_LCMS_Indiv.csv",header=T) %>% select(-c(pos_age,pos,pop,line,age,lat,region))

all.traits <- full_join(tough,cn,by="line_age") %>% full_join(palat,by="line_age") %>% full_join(chemsums,by="line_age") %>% full_join(compounds,by="line_age")

key <- as.data.frame(all.traits$line_age)
colnames(key) <- "line_age"
key <- cbind(key, read.table(text = as.character(key$line_age), sep = "_", colClasses = "character"))
colnames(key)  <- c("line_age","pop","line","age")
key$line <- paste(key$pop, key$line, sep="_")
key <- merge(key, read.csv("Raw/LatsPopsKey.csv",header=T),by="pop")

all.traits <- left_join(key,all.traits,by="line_age")

area_chem_only <- all.traits[,c(13,23:1939)]

write.csv(all.traits, "Processing/2_out_AllTraits.csv",row.names=F)
#write.csv(all.traits, "Processing/2_out_AllTraitsTab.csv",row.names=F,sep="\t")
write.table(drop_na(all.traits), "Processing/2_out_AllTraitsTab.csv",row.names=F,sep="\t")

write.table(drop_na(area_chem_only), "Processing/2_out_AreaChemTab.csv",row.names=F,sep="\t")

rm(list=ls())
