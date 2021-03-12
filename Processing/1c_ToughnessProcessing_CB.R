#toughness wrangling
library(tidyverse)
raw = read.csv("Raw/20171005_ToughnessAge_RAW.csv", header=TRUE) %>% mutate(young = (young.1+young.2)/2, mature = (mature.1+mature.2)/2) %>% select(-c(young.1,young.2,mature.1,mature.2))
str(raw)

#remove some lines where we have both bulk and non-bulk from the same maternal line, I don't really think these are a good idea
raw <- filter(raw, !pop %in% c("Tiri Bulk", "KBS Bulk","TT Bulk","Whitehall Bulk","Gav Bulk"))
#combine Valverde with Tiri (nearby) and remove bulk label from two others. I am a little nervous about these bulked plants because I don't remember exactly how it was done--did we just use fruits that we found on the plants? How did we know who the father was? But at worst, if we used fruits without knowing whether they were crossed by stray wasps, we would be eroding population differentiation and adding noise to the data.
raw$pop <- plyr::revalue(raw$pop, c("Valverde" = "Tiri", "Hwy 27 bulk" = "Hwy 27", "Dalton Bulk" = "Dalton")) 

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
