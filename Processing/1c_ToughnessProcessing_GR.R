#toughness wrangling
raw = read.csv("Raw/20171005_ToughnessAge_RAW.csv", header=TRUE)
library(tidyverse)
raw$pop <- plyr::revalue(raw$pop, c("Valverde" = "Tiri","Tiri Bulk" = "Tiri", "Hwy 27 bulk" = "Hwy 27", "TT Bulk" = "TT", "Dalton Bulk" = "Dalton", "Gav Bulk" = "Gav", "KBS Bulk" = "KBS", "Whitehall Bulk" = "Whitehall")) 

raw$young = (raw$young.1+raw$young.2)/2
raw$mature = (raw$mature.1+raw$mature.2)/2

#reshaping data from wide to long format
df <- reshape(raw,
              varying = c("young", "mature"),
              v.names = "tough",
              timevar = "age",
              times = c("young", "mature"),
              direction = "long")
df$mature.1 <- NULL
df$mature.2 <- NULL
df$young.1 <- NULL
df$young.2 <- NULL
df$id <- NULL
df$line <- as.factor(paste(df$pop, df$line, sep="_"))
df$line_age <- as.factor(paste(df$line, df$age, sep="_"))

#line-level means
tough_line <- plyr::ddply(df, c("line_age"), summarise,
                 tough =mean(tough)
)

write.csv(tough_line,"Processing/1c_out_ToughnessProcessed_CB.csv",row.names=F)

rm(list=ls())
