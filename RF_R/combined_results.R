#exploring RF results

library(tidyverse)

###############
#comparing lists of compounds found in young and mature best models
###############
#best mature has 9 features and young has 24

#can't figure out how to do these in one step because the ranking refers to the row names of the table and it has to exist already
young <- read.csv("RF_R/youngfine_24_imp",sep = "\t", header = T) %>% rename(feature=X,imp=mean_imp) 
young <- young %>% add_column(final_rank = as.numeric(row.names(young)), age= "young")
mature <- read.csv("RF_R/maturefine_9_imp",sep = "\t", header = T) %>% rename(feature=X,imp=mean_imp) 
mature <- mature %>% add_column(final_rank = as.numeric(row.names(mature)), age= "mature")

#number of unique features
length(unique(bind_rows(mature,young)$feature))

#features on both lists
(overlap = intersect(young$feature,mature$feature))

#2 compounds and latitude overlap
mature[mature$feature %in% overlap==T,]
young[young$feature %in% overlap==T,]
#two of the top 3 in mature (not #1) are also found to be important chemicals in young leaves
rm(overlap)
###############
#how did the non-chemical features rank in the maximal model? And how did the chemicals in the final model rank in the original model? I was going to ask when did the non-chemical features drop out during model selection, but realized that maybe answering this easier question is sufficient. Especially when you see how little the rankings actually change--the original model is amazingly good at picking up signal through the noise!
###############
nonchem <- c("tough","percent_N","C_N","log.abund","richness","diversity","reg.north_temperate","reg.subtropical","reg.temperate","reg.tropical","spp.PHRI")

orig <- bind_rows(mature,young) %>% select(-imp) %>% add_row(feature = nonchem[2:length(nonchem)], final_rank = NA, age= "young") %>% add_row(feature = nonchem, final_rank = NA, age= "mature") %>% arrange(age) %>% add_column(orig_rank = NA)

maturemax <- read.csv("RF_R/mature_max_imp",sep = "\t", header = T)
youngmax <- read.csv("RF_R/young_max_imp",sep = "\t", header = T)
for(i in 1:nrow(orig)){
  if (orig$age[i] == "mature") {
    orig$orig_rank[i] = as.numeric(row.names(maturemax[maturemax$X == orig$feature[i],]))
  } else {
    orig$orig_rank[i] = as.numeric(row.names(youngmax[youngmax$X == orig$feature[i],]))
  }
}
rm(maturemax,youngmax,i,nonchem)

orig <- orig %>% select(age,feature,orig_rank,final_rank) %>% mutate(change = orig_rank - final_rank) %>% arrange(orig_rank) %>% arrange(age) %>% rename("Leaf age" = age, Predictor = feature, "Original Rank" = orig_rank, "Final Rank" = final_rank, "Change in Rank" = change) %>% add_column(Class = "", Notes = "")
#the final 9 features for mature were in the top 33 in the original model (top 2%). The final 24 features for young were in the top 208 (11%); 23 of those were in the top 127 (7%)

write.csv(orig,"FiguresTables/Table_RF_FirstLastModels.csv",row.names=F)
rm(list=ls())
