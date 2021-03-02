#started 4-23-28
#toughness wrangling
#setwd("/Users/carina/Documents/R_working_directory")
raw = read.csv("data_in/20171005_toughness_age.csv", header=TRUE)
raw = na.omit(raw)
library(plyr)
raw$pop <- revalue(raw$pop, c("Valverde" = "Tiri","Tiri Bulk" = "Tiri", "Hwy 27 bulk" = "Hwy 27", "TT Bulk" = "TT", "Dalton Bulk" = "Dalton", "Gav Bulk" = "Gav", "KBS Bulk" = "KBS", "Whitehall Bulk" = "Whitehall")) 

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
tough_line <- ddply(df, c("line_age"), summarise,
                 tough =mean(tough)
)
rm(df, raw)
