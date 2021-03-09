#carbon:nitrogen data wrangling
#combined the KBS and main campus data in Excel
library(tidyverse)

cn = read.csv("Raw/20180422_CarbonNitrogen_RAW.csv", header=TRUE) %>% mutate(indiv = paste(paste(pop,line,sep="_"),indiv,sep="_"))
str(cn)
#check how many observations per individual
length(unique(cn$indiv))*2 == nrow(cn)
count(cn,age)
#two individuals didn't have young leaves. But we don't have to take individual means, there's just the usual one young and one mature leaf measured per individual

#I will make up a pos just to keep consistency with other datasets
cn <- cn %>% mutate(pos = seq(1:nrow(cn)), pos_age = paste(pos,age,sep="_"), line_age = paste(paste(pop,line,sep="_"),age,sep="_"))

lat = read.csv("Raw/LatsPopsKey.csv", header = TRUE)
cn = merge(cn, lat, by= "pop")

write.csv(cn[,c(9,10,8,1,2,4,11,12,5:7)],"Processing/1d_out_CarbonNitrogen_Indiv.csv",row.names=F)

cn_line <- cn %>% group_by(line_age) %>% summarise(percent_N = mean(percent_N), percent_C = mean(percent_C), C_N = mean(C_N))

write.csv(cn_line, "Processing/1d_out_CarbonNitrogen_Line.csv",row.names=F)

rm(list=ls())

