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
#log.area takes natural log of area

working <- raw %>% mutate(line = paste(pop,line,sep="_"), 
                          line_age = paste(line,age,sep="_"), 
                          surv = ifelse(age=="mature",surv_d8,surv_d9), 
                          initial = ifelse(age=="mature",5,3), 
                          duration = ifelse(age=="mature",8,9),
                          cons = ifelse(age=="mature",
                                        cons_d1+cons_d2+cons_d3+cons_d4+cons_d5+cons_d6+cons_d7+cons_d8,
                                        cons_d1+cons_d3+cons_d5+cons_d7+cons_d9),
                          area = cons*10.080625,
                          log.area = log(area+1)
                          ) %>% select(-c(cons_d1,cons_d2,cons_d3,cons_d4,cons_d5,cons_d6,cons_d7,cons_d8,cons_d9,surv_d1,surv_d2,surv_d3,surv_d4,surv_d5,surv_d6,surv_d7,surv_d8,surv_d9,cons,area))

sum(is.na(working$surv))
sum(is.na(working$area))
#every cup has survival and only 10 are missing area

#################################################
#now import caterpillar biomass data
mass <- read.csv("Raw/2016_PalatabilityBiomass_RAW.csv",header=T)
#it's a little complicated because occasionally a caterpillar was lost during measurement (those things were TINY!), so we can get the mean biomass per cup (mass_surv) and then multiply that by the number of real survivors from the other datasheet to get mass_cat (so we interpolate the few missing caterpillars as having the same biomass as their cup-mates)
#after looking at notes from old code, there are 4 caterpillars that were lost during measurement, and the cupID matches these (56,84,98,222)
sum(is.na(mass_cat))
mass[which(is.na(mass$mass_mg)==T),]
#list of cupID for cups with no survivors; these should not be found in mass
nosurv <- working[working$surv==0,]$cupID
sum(mass$cupID %in% nosurv)
mass[which(mass$cupID %in% nosurv == T),]
#ok that's odd, there is data for one caterpillar from a cup that should have zero survivors. The survivor count or tube label was a mistake! I will change the mass to NA
mass[which(mass$cupID %in% nosurv == T),4] <- NA

#calculate mean mass per survivor in each cup
mass_cat <- mass %>% group_by(cupID) %>% summarise(mass_surv = mean(mass_mg,na.rm=T))
sum(is.na(mass_cat$mass_surv))
#there are still two NA from these mistake cups (222 and 590) that had only one caterpillar and it was lost. After merging there will also be NA for cups with zero survivors.





abc <- working %>% left_join(mass_cat,by="cupID") %>% mutate(mass_cup = ifelse(surv == 0, 0, mass_surv*surv))

#now standardize the area and mass by the number of starting caterpillars and duration
#not sure where to do the log-transform: before or after standardizing by cats and duration?