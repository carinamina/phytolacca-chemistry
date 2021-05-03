#create a heatmap of compounds
library(ggplot2)
library(tidyverse)
library(reshape)
library(viridis)
library(ggh4x)
library(ggmap)
library(maps)
library(mapdata)

#import chemical data
chem_raw <- read.csv("Processing/1a_out_LCMS_indiv.csv",header=T)
#import the list of populations and their coordinates and assign letters because this will be used in multiple plots
coord <- read.csv("Raw/LatsPopsKey.csv",header=T) %>% filter(pop %in% unique(chem_raw$pop)) %>% arrange(-lat) %>% add_column(letter = LETTERS[seq(1, 16)])

#################
# Panel 1: map
#################
world <- map_data("world")

basemap = ggplot(data = world) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = NA, color = "grey") + 
  guides(fill=FALSE) + theme_classic() + coord_fixed(xlim = c(-87, -76),  ylim = c(8.5, 43), ratio = 1.3) 
basemap

sitemap = basemap + geom_text(data = coord, aes(x = long, y = lat, label = recode(letter, O="O-P",P=""), color = region), size = 4) + scale_color_manual(values = c("steelblue1","navyblue","gray44","maroon2")) + xlab("Longitude") + ylab("Latitude") + theme(legend.position = "none")
sitemap

# setEPS()
# postscript("20190122map.eps", height = 7)
# sitemap
# dev.off()

#################
# Panel 2: palatability
#################
#import data. There might be a problem with this that I just realized: I took the log-area eaten for each cup and then a line mean, which was used in RF analysis. Should I have done the log-transformation after taking the means? Shit. Here also I am taking the mean of a log-transformed variable.

palat <- read.csv("Processing/2_out_AllTraits.csv",header=T) %>% select(pop,age,area) %>% group_by(pop,age) %>% summarise(area=mean(area,na.rm=T)) %>% left_join(coord,by="pop")

dots <- ggplot(palat, aes(x=1,y=letter)) +
  geom_point(aes(size=area)) +
  facet_grid(age~.,scales="free_y") +
  theme_bw() +theme(panel.spacing = unit(.1, "lines"), strip.background = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank()) +
  xlab(label ="Palatability") +
  ylab(label = "Site (see map)")

dots

#export
# setEPS()
# postscript("FiguresTables/HeatMap.eps")
# dots + theme()
# dev.off()


#################
# Panel 3: heatmap
#################
#take log of chemical abundances
chem <- bind_cols(left_join(chem_raw[1:8],select(coord,c(pop,letter)),by="pop"),log(chem_raw[9:ncol(chem_raw)]+1)) %>% select(-c(pos,line,lat)) 

#%>% arrange(line_age) %>% arrange(desc(lat)) %>% mutate(lat = round(lat, digits = 2))
#we need two digits for lat to distinguish Bella and Gav
#chem$lat <- as.character(chem$lat)

#this prepares for the geom_tile. For troubleshooting, use a subset:
#chem_melt <- melt(chem[1:100])
chem_melt <- melt(chem)
head(chem_melt)
#(lat_order <- rev(levels(as.factor(as.numeric(chem_melt$lat)))))
chem_melt$age <- as.factor(chem_melt$age) %>% recode(young = "young leaves", mature = "mature leaves")

heat <- ggplot(chem_melt, aes(x=variable,y=pos_age)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_viridis(direction=-1) + 
  facet_nested(age + letter ~.,scales="free_y",switch="y",nest_line=T) + 
  theme_bw() +theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.spacing = unit(.1, "lines"), strip.background = element_blank())  +
  xlab(label="LC/MS Peak (from low to high retention time)") + 
  ylab(label="Site (see map)")

heat

#putting this in theme() rotates both the age and site, and the rotated letters have a lot of space around them
#, strip.text.y.left = element_text(angle = 0)

#export
# setEPS()
# postscript("FiguresTables/HeatMap.eps")
# heat
# dev.off()

