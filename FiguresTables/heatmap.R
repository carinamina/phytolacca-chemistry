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
rm(basemap,world)
#################
# Panel 2: palatability (dots)
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

rm(palat)
#################
# Panel 3: heatmap: all compounds
#################
#take log of chemical abundances
chem <- bind_cols(left_join(chem_raw[1:8],select(coord,c(pop,letter)),by="pop"),log(chem_raw[9:ncol(chem_raw)]+1)) %>% select(-c(pos,line,lat)) 

#this prepares for the geom_tile. For troubleshooting, use a subset [1:100}:
#chem_melt <- melt(chem[1:100])
chem_melt <- melt(chem)
head(chem_melt)
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

rm(heat, dots, sitemap, chem_melt)

#################
# Second figure: heatmap: important compounds
#################
#subset chem list to compounds found to be important from RF
young <- (read.csv("RF_R/youngfine_45_features.txt",sep = "\t", header = F) %>% filter(V1 != "lat") %>% arrange(V1))$V1
mature <- (read.csv("RF_R/mature_15_features.txt",sep = "\t", header = F) %>% arrange(V1) )$V1

#subset log-abund for each leaf age
imp_y <- select(chem, c(pos_age, line_age, pop, age, region, letter, young)) %>% filter(age == "young")
names(imp_y) <- gsub("X", "", names(imp_y))
imp_m <- select(chem, c(pos_age, line_age, pop, age, region, letter, mature)) %>% filter(age == "mature")
names(imp_m) <- gsub("X", "", names(imp_m))

#figure for mature leaves
melt_m <- melt(imp_m)
head(melt_m)
heat_m <- ggplot(melt_m, aes(x=variable,y=pos_age)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_viridis(direction=-1) + 
  facet_grid(letter ~.,scales="free_y",switch="y") + 
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(.1, "lines"), strip.background = element_blank(), strip.text.y.left = element_text(angle = 0))  +
  xlab(label="LC/MS Peak (from low to high retention time)") + 
  ylab(label="Site (see map)")

heat_m

#figure for young leaves
melt_y <- melt(imp_y)
head(melt_y)
heat_y <- ggplot(melt_y, aes(x=variable,y=pos_age)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_viridis(direction=-1) + 
  facet_grid(letter ~.,scales="free_y",switch="y") + 
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(.1, "lines"), strip.background = element_blank(), strip.text.y.left = element_text(angle = 0))  +
  xlab(label="LC/MS Peak (from low to high retention time)") + 
  ylab(label="Site (see map)")

heat_y

########## trying to make colored boxes around the site letters
#I found this example here https://github.com/tidyverse/ggplot2/issues/2096 but I can't get it to work when applying it to heat_y, because I have no idea what it's doing
p <- ggplot(mpg, aes(displ, cty)) + geom_point() + facet_grid(drv ~ cyl)
p
g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("red","green","blue","yellow","red","green","blue","yellow")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g)

g <- ggplot_gtable(ggplot_build(heat_y))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("red","green","blue","yellow","orange","purple","pink","white","red","green","blue","yellow","orange","purple","pink","white")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g)
