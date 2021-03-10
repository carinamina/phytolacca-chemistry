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
#LogAreaCat_day: first get area consumed per caterpillar (using initial counts and ignoring caterpillar death). That should have an exponential distribution, so we log-transform it. Then standardize by the experiment duration. 
#unfortunately, in the Ecology Letters paper, I did log(area)/initial/duration. This might be a problem because initial differs by leaf age, which was an explanatory variable of interest. So I will include a column with the "old calculation" to compare results later during analysis.

working <- raw %>% mutate(line = paste(pop,line,sep="_"), 
                          line_age = paste(line,age,sep="_"), 
                          surv = ifelse(age=="mature",surv_d8,surv_d9), 
                          initial = ifelse(age=="mature",5,3), 
                          duration = ifelse(age=="mature",8,9),
                          cons = ifelse(age=="mature",
                                        cons_d1+cons_d2+cons_d3+cons_d4+cons_d5+cons_d6+cons_d7+cons_d8,
                                        cons_d1+cons_d3+cons_d5+cons_d7+cons_d9),
                          area = cons*10.080625,
                          LogAreaCat_day = log(area/initial+1)/duration,
                          old.LogArea_cat_day = log(area+1)/initial/duration
                          ) %>% select(-c(cons_d1,cons_d2,cons_d3,cons_d4,cons_d5,cons_d6,cons_d7,cons_d8,cons_d9,surv_d1,surv_d2,surv_d3,surv_d4,surv_d5,surv_d6,surv_d7,surv_d8,surv_d9,cons,area))

sum(is.na(working$surv))
sum(is.na(working$LogAreaCat_day))
#every cup has survival and only 10 are missing area

#################################################
#import caterpillar biomass data
#Note the log-transforms. For average mass of survivors, we don't need to standardize by the number of caterpillars; it's already included in the mean. Log-transform the biomass and THEN standardize by duration. For total biomass in the cup, which counts dead caterpillars as zero biomass, we calculate biomass per initial caterpillar count, log-transform that, and THEN standardize by duration. We want to log-transform the exponential growth, which happens for each caterpillar, and then standardize by the duration.
#unfortunately, in the Ecology Letters and Ecological Monographs papers, I did log(biomass)/caterpillars/duration. This might be a problem because initial differs by leaf age, which was an explanatory variable of interest. So I will include columns with the "old calculation" (earlylog) to compare results later during analysis.
mass <- read.csv("Raw/2016_PalatabilityBiomass_RAW.csv",header=T) %>% mutate(earlylogmass = log(mass_mg+1))
#occasionally a caterpillar was lost during measurement (those things were TINY!). CupIDs (from old code) were 56,84,98,222
mass[which(is.na(mass$mass_mg)==T),]
#list of cupID for cups with no survivors; these should not be found in mass
nosurv <- working[working$surv==0,]$cupID
sum(mass$cupID %in% nosurv)
mass[which(mass$cupID %in% nosurv == T),]
#ok that's odd, there is data for one caterpillar from a cup that should have zero survivors. The survivor count or tube label was a mistake! I will change the mass to NA
mass[which(mass$cupID %in% nosurv == T),c(4,5)] <- NA

#calculate mean mass per survivor in each cup and total mass of each cup (those which are missing caterpillars should be NA for the sum)
cup_mass <- mass %>% group_by(cupID) %>% summarise(mass_surv = mean(mass_mg,na.rm=T), earlylog_mass_surv = mean(earlylogmass,na.rm=T), mass_cup = sum(mass_mg), earlylog_mass_cup = sum(earlylogmass))
cup_mass[which(is.na(cup_mass$mass_cup)==T),]
#shows all the NA, which we already discussed. I'm not sure whether NaN is treated differently from NA so I will just manually change those two entries
cup_mass[which(cup_mass$cupID == 222 | cup_mass$cupID == 590),c(2,3)] <- NA

#merge the caterpillar biomass with the rest of the data and add latitude and region
working <- working %>% left_join(cup_mass,by="cupID") %>% mutate(
                            LogMassSurv_day = log(mass_surv+1)/duration,
                            old.LogMassSurv_day = earlylog_mass_surv/duration,
                            LogCupCat_day = log(ifelse(surv == 0, 0, mass_cup)/initial+1)/duration,
                            old.LogCup_cat_day = ifelse(surv == 0, 0, earlylog_mass_cup)/initial/duration
                            ) %>% select(-c(duration,mass_surv,mass_cup,earlylog_mass_surv,earlylog_mass_cup)
                            ) %>% left_join(read.csv("Raw/LatsPopsKey.csv",header=T),by="pop")

sum(is.na(working$LogAreaCat_day))
#10 missing area when leaves were not measured on one day
sum(is.na(working$LogCupCat_day))
#4 missing a cup mass (cup 590 got turned into 0 because of 0 survivors)
sum(is.na(working$LogMassSurv_day))
#50 missing survivor mass (222 and 49 cups with no survivors)

#rename responses for sake of brevity
working <- working %>% dplyr::rename(area = LogAreaCat_day, mass_surv = LogMassSurv_day, mass_cup = LogCupCat_day, old.mass_cup = old.LogCup_cat_day, old.area = old.LogArea_cat_day, old.mass_surv = old.LogMassSurv_day) %>% select(line_age, cupID, pop, line, age, lat, region, initial, surv, area, old.area, mass_surv, old.mass_surv, mass_cup, old.mass_cup)

#export cup-level data
write.csv(working,"Processing/1e_out_Palat_Cup.csv",row.names=F)

#export line-level data. I won't include the "old" calculations because I think I won't compare them for the line-level analysis. I remove NA bc missing data for one cup shouldn't invalidate a whole line
line_palat <- working %>% group_by(line_age) %>% summarise(initial = mean(initial), surv = mean(surv), area = mean(area,na.rm=T), mass_surv = mean(mass_surv,na.rm=T), mass_cup = mean(mass_cup,na.rm=T))
sum(is.na(line_palat$mass_surv))
#there are still some NaN in mass_surv because these lines had zero survivors
write.csv(line_palat, "Processing/1e_out_Palat_Line.csv", row.names=F)

#some visualization for fun
qplot(lat, old.area, data=working, colour=age)
qplot(area, mass_cup, data = working, color = age)
#visual comparisons of old and new plots of biomass vs latitude: I can hardly tell the difference, but maybe the old plots slightly under-estimated the leaf age effect
#comparing old and new plots of area vs latitude, which hasn't been published, the old calculations clearly underestimate leaf age differences.
#comparing old and new plots of mass_cup vs. area, there might be a problem, but I *think* it would just cause a shift in the intercept of the fitted relationships, but the slopes would still differ. It's really hard to tell. Definitely should re-analyze this and possibly issue a correction to the Ecology Letters paper. Luckily that plot (Fig. 4B) was hardly discussed at all.
#mean leaf area consumed by age for old and new for PHAM (reported in Ecology Letters paper)
working %>% filter(region != "tropical") %>% group_by(age) %>% summarise(old.area = mean(old.area,na.rm=T), new.area = mean(area,na.rm=T))
#the numbers for old are within .007 of the numbers in the paper, so it seems that I did a decent job of replicating the calculations (why are they are not identical though, I don't know). The difference between young and mature is actually reversed with the new calculations!

rm(list=ls())
