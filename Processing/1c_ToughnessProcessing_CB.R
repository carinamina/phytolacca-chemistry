#toughness wrangling
library(tidyverse)
raw = read.csv("Raw/20171005_ToughnessAge_RAW.csv", header=TRUE)
str(raw)
raw$pop <- plyr::revalue(raw$pop, c("Valverde" = "Tiri","Tiri Bulk" = "Tiri", "Hwy 27 bulk" = "Hwy 27", "TT Bulk" = "TT", "Dalton Bulk" = "Dalton", "Gav Bulk" = "Gav", "KBS Bulk" = "KBS", "Whitehall Bulk" = "Whitehall")) 

raw$young = (raw$young.1+raw$young.2)/2
raw$mature = (raw$mature.1+raw$mature.2)/2

lats <- read.csv("Raw/LatsPopsKey.csv", header = TRUE)

#reshaping data from wide to long format
df <- reshape(raw,
              varying = c("young", "mature"),
              v.names = "tough",
              timevar = "age",
              times = c("young", "mature"),
              direction = "long") %>% select(pos,pop,line,age,tough) %>% mutate(line_age = paste(paste(pop,line,sep="_"),age,sep="_"),.before=pos) %>% mutate(pos_age = paste(pos,age,sep="_"),.before=line_age) %>% left_join(lats,by="pop")

#export indiv means
write.csv(df[,c(1:6,8,9,7)], "Processing/1c_out_Toughness_Indiv.csv",row.names=F)

#line-level means
tough_line <- df %>% group_by(line_age) %>% summarise(tough = mean(tough))

write.csv(tough_line,"Processing/1c_out_Toughness_Line.csv",row.names=F)

rm(list=ls())
