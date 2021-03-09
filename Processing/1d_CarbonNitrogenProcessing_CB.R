#carbon:nitrogen data wrangling
#combined the KBS and main campus data in Excel
library(tidyverse)

cn = read.csv("Raw/20180422_CarbonNitrogen_RAW.csv", header=TRUE)
str(cn)

# **pop** is the population  
# **line** is the maternal line ID nested within pop
# **indiv** is the plant ID letter nested within line  
# **age** is leaf age, young or mature  
# **percent_N** is the percentage of N in sample  
# **percent_C** is the percentage of C in sample  
# **C_N** is percent C/percent N for each sample
# 
# ####New variables  
# **line_ID** is a unique ID for each line with the population name, an underscore, and the line number  
# **plant_ID** is a unique ID for each plant with the population name, an underscore, the line number, an underscore, and the indiv letter  
# **lat** is latitude of each site (which is actually not necessary for this particular analysis, which uses region)  

cn$line_ID <- as.factor(paste(cn$pop,cn$line,sep="_"))
cn$plant_ID <- as.factor(paste(cn$line_ID,cn$indiv,sep="_"))
cn$pop_age <- as.factor(paste(cn$pop, cn$age, sep = "_"))

lat = read.csv("Raw/LatsPopsKey.csv", header = TRUE)
cn = merge(cn, lat, by= "pop")

#assigns regional names based on latitude. 3 tropical, 3 southernmost US, 3 central (based on 3 closest to mean of lowest and highest US) 3 northernmost US
regions <- c("tropical", "tropical", "tropical", "southern", "southern", "southern", NA, NA, NA, NA, NA, NA, NA, "northern", "northern", "northern")
regions <- cbind(regions, sort(unique(cn$lat)))
colnames(regions) <- c("region", "lat")
cn <- merge(cn, regions, by = "lat")

#get line-level C:N
cn_line <- plyr::ddply(cn, c("region","line_ID","pop","lat", "age"), summarise,
                C_N = mean(C_N),
                percent_N = mean(percent_N),
                percent_C = mean(percent_C) 
)

#add palatability data
cn_line$line_age <- as.factor(paste(cn_line$line_ID, cn_line$age, sep = "_"))
cn_line$line_ID <- NULL
palat = read.csv("Raw/LinePalatability.csv", header = TRUE)
palat$line_age = paste(palat$line, palat$age, sep="_")
palat$X <- NULL
palat$pop <- NULL
palat$lat <- NULL
palat$age <- NULL
palat_cn <- merge(palat, cn_line, by = "line_age", all = TRUE)
str(palat_cn)
#remove anything that doesn't have palatability data
palat_cn <- palat_cn[!is.na(palat_cn$biomass),]
#biomass is ln(total biomass per cup) (counting dead cats as zero biomass), standardized by intitial cat count and duration
#area is ln(cumulative area consumed per cup) over course of experiment, standardized by intitial cat count and duration
#survival is number of survivors
rm(cn_line,lat,palat,regions)

