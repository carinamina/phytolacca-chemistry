#going waaaaay back to the raw palatability data...I don't like that there are already calculated columns in the dataset used for the tropical paper, I want to start over with the raw data

library(tidyverse)

#read in data, remove controls (fed artificial diet) and some experimental unit info we don't need
raw = read.csv("Raw/2016_PalatabilitySurvArea_RAW.csv", header=TRUE) %>% filter(pop != "control") %>% select(-c(line_pos,rep,tray))
#paste pop and line
#raw$line <- paste(raw$pop, raw$line, sep="_")
str(raw) 

#create new variables:
#surv is number of survivors on the last day (day 8 for mature, day 9 for young)
#initial is initial caterpillar count (5 for mature, 3 for young)
#duration is duration of experiment (8 days for mature, 9 days for young)
#cons is total leaf area consumed, there are some NAs (I think just mistakes where leaves were accidentally tossed before measuring). I want those cups to have NA for the sum, but when taking line means we can use na.rm=T so the whole line isn't NA. Young and mature leaves were measured on different days. Young on odd days, mature every day except day9
#area converts cons to square milimeters (it was measured by counting small quilting squares)
#eventually we need to log-transform stuff, but that happens later, after taking mean per line.

working <- raw %>% mutate(line = paste(pop,line,sep="_"), 
                          line_age = paste(line,age,sep="_"), 
                          surv = ifelse(age=="mature",surv_d8,surv_d9), 
                          initial = ifelse(age=="mature",5,3), 
                          duration = ifelse(age=="mature",8,9),
                          cons = ifelse(age=="mature",
                                        cons_d1+cons_d2+cons_d3+cons_d4+cons_d5+cons_d6+cons_d7+cons_d8,
                                        cons_d1+cons_d3+cons_d5+cons_d7+cons_d9),
                          area = cons*10.080625
                          ) %>% select(-c(cons_d1,cons_d2,cons_d3,cons_d4,cons_d5,cons_d6,cons_d7,cons_d8,cons_d9,surv_d1,surv_d2,surv_d3,surv_d4,surv_d5,surv_d6,surv_d7,surv_d8,surv_d9,cons))

sum(is.na(working$surv))
sum(is.na(working$area))
#every cup has survival and only 10 are missing area due to missing leaf measurements on one day

#add latitude and region, re-order
working <- working %>% left_join(read.csv("Raw/LatsPopsKey.csv",header=T),by="pop") %>% select(line_age, cupID, pop, line, age, lat, region, initial, duration, surv, area)

#export cup-level data for geographic analysis
#area clearly has an exponential distribution, which makes sense because cats can eat exponentially more as they grow
hist(working$area)
#so we need to log-transform, but it's a bit complicated. First we want area consumed per caterpillar, which is still exponentially distributed because each caterpillar eats in an exponential manner
hist(working$area/working$initial)
#we log-transform that, and THEN standardize by the experiment duration. Unfortunately, in the Ecology Letters paper, I did log(area)/initial/duration. This might be a problem because initial differs by leaf age, which was an explanatory variable of interest. So I will include a column with the "old calculation" to compare results later during analysis.
#furthermore, we need to do this log-transformation and standardization after taking a mean, not before, so that's why I do it twice: here for the cup-level data and again below for line-level data.
working_logs <- working %>% mutate(log.area = log(area/initial+1)/duration, old.log.area = log(area+1)/initial/duration) %>% select(-c(area, duration, initial))
write.csv(working_logs,"Processing/1e_out_Palat_Cup.csv",row.names=F)


#take line-level means for trait analysis. I won't include the "old" calculations because I think I won't compare them for the line-level analysis. I remove NA bc missing data for one cup shouldn't invalidate a whole line
line_palat <- working %>% group_by(line_age) %>% summarise(initial = mean(initial), duration = mean(duration), surv = mean(surv), area = mean(area,na.rm=T)) %>% mutate(log.area = log(area/initial+1)/duration) %>% select(-c(area, duration, initial))

write.csv(line_palat, "Processing/1e_out_Palat_Line.csv", row.names=F)

#some visualization for fun
qplot(lat, old.log.area, data=working_logs, colour=age)
qplot(lat, log.area, data=working_logs, colour=age)
#visual comparisons of old and new plots of area vs latitude:  the old plots slightly under-estimated the leaf age effect
#mean leaf area consumed by age for old and new for PHAM (reported in Ecology Letters paper)
working_logs %>% filter(region != "tropical") %>% group_by(age) %>% summarise(old.log.area = mean(old.log.area,na.rm=T), log.area = mean(log.area,na.rm=T))
#the numbers for old are within .007 of the numbers in the paper, so it seems that I did a decent job of replicating the calculations (why are they are not identical though, I don't know). The difference between young and mature is actually reversed with the new calculations!

rm(list=ls())
