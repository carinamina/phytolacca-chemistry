#create a heatmap of compounds
library(ggplot2)
library(tidyverse)
library(reshape)
library(viridis)
library(ggmap)
library(maps)
library(mapdata)

#import chemical data
chem_raw <- read.csv("Processing/1a_out_LCMS_indiv.csv",header=T)

#import the list of populations and their coordinates and assign letters because this will be used in multiple plots
coord <- read.csv("Raw/LatsPopsKey.csv",header=T) %>% filter(pop %in% unique(chem_raw$pop)) %>% arrange(-lat) %>% add_column(letter = LETTERS[seq(1, 16)])

#take log of chemical abundances
chem <- bind_cols(left_join(chem_raw[1:8],select(coord,c(pop,letter)),by="pop"),log(chem_raw[9:ncol(chem_raw)]+1)) %>% select(-c(pos,line,lat)) 


#################
# Panel 1: map
#################
world <- map_data("world")

basemap = ggplot(data = world) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = NA, color = "grey") + 
  guides(fill=FALSE) + theme_classic() + coord_fixed(xlim = c(-87, -80),  ylim = c(8.5, 43), ratio = 1.3) 
basemap

sitemap = basemap + geom_text(data = coord, aes(x = long, y = lat, label = recode(letter, O="O-P",P=""), color = region), size = 4) + scale_color_manual(values = c("steelblue1","navyblue","gray44","maroon2")) + xlab("Longitude") + ylab("Latitude") + theme(legend.position = "none")
sitemap

setEPS()
postscript("FiguresTables/SamplingMap.eps", height = 10)
sitemap
dev.off()

sitemap_key = basemap + geom_point(data = coord, aes(x = long, y = lat, color = region), size = 4) + scale_color_manual("Region" ,values = c("steelblue1","gray44","navyblue","maroon2"),breaks=c("north temperate","temperate","subtropical","tropical")) + xlab("Longitude") + ylab("Latitude")
sitemap_key

setEPS()
postscript("FiguresTables/SamplingMapKey.eps", height = 10)
sitemap_key
dev.off()

rm(basemap,world)

#################
# Panel 2: heatmap: all compounds
#################
#function to color-code the site boxes
popcolor <- function(plotname)
{
  g <- ggplot_gtable(ggplot_build(plotname))
  strip_side <- which(grepl('strip-l', g$layout$name))
  fills <- c("steelblue1","steelblue1","steelblue1","gray60","gray60","gray60","gray60","gray60","gray60","gray60","navyblue","navyblue","navyblue","maroon2","maroon2","maroon2")
  k <- 1
  for (i in strip_side) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  grid::grid.draw(g)
  rm(g,i,j,k,fills,strip_side)
}

all_m <- melt(subset(chem, chem$age =="mature"))
head(all_m)
heat_all_m <- ggplot(all_m, aes(x=variable,y=pos_age)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_viridis(direction=-1,limits=c(0,10.7)) + 
  facet_grid(letter ~.,scales="free_y",switch="y") + 
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.spacing = unit(0, "lines"), strip.text.y.left = element_text(angle = 0,color="white"), axis.title=element_text(size=18))  +
  xlab(label="LC/MS Peak") + 
  ylab(label="Mature leaves") 

popcolor(heat_all_m)

setEPS()
postscript("FiguresTables/HeatMapMatureAll.eps", height = 5, width = 9)
popcolor(heat_all_m + theme(legend.position = "none"))
dev.off()

all_y <- melt(subset(chem, chem$age =="young"))
head(all_y)
heat_all_y <- ggplot(all_y, aes(x=variable,y=pos_age)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_viridis(direction=-1,limits=c(0,10.7)) + 
  facet_grid(letter ~.,scales="free_y",switch="y") + 
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.spacing = unit(0, "lines"), strip.text.y.left = element_text(angle = 0,color="white"), axis.title=element_text(size=18))  +
  xlab(label="LC/MS Peak") + 
  ylab(label="Young leaves") 

popcolor(heat_all_y)

setEPS()
postscript("FiguresTables/HeatMapYoungAll.eps", height = 5, width = 9)
popcolor(heat_all_y + theme(legend.position = "none"))
dev.off()

#################
# Second figure: heatmap: important compounds
#################
#subset chem list to compounds found to be important from RF
young <- (read.csv("RF_R/youngfine_24_features.txt",sep = "\t", header = F) %>% filter(V1 != "lat") %>% arrange(V1))$V1
mature <- (read.csv("RF_R/maturefine_9_features.txt",sep = "\t", header = F) %>% filter(V1 != "lat", V1 != "log.abund") %>% arrange(V1) )$V1

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
  scale_fill_viridis(direction=-1,limits=c(0,10.7)) + 
  facet_grid(letter ~.,scales="free_y",switch="y") + 
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 50, hjust = 1), panel.spacing = unit(0, "lines"), strip.text.y.left = element_text(angle = 0,color="white"), axis.title=element_text(size=18))  +
  xlab(label="LC/MS Peak") + 
  ylab(label="Mature leaves") 

popcolor(heat_m)

setEPS()
postscript("FiguresTables/HeatMapMatureImportant.eps", height = 5, width = 2.3)
popcolor(heat_m + theme(legend.position = "none"))
dev.off()

setEPS()
postscript("FiguresTables/HeatMapKey.eps")
heat_m + theme(legend.position = "top",legend.title = element_blank()) 
dev.off()

#figure for young leaves
melt_y <- melt(imp_y)
head(melt_y)
heat_y <- ggplot(melt_y, aes(x=variable,y=pos_age)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_viridis(direction=-1,limits=c(0,10.7)) + 
  facet_grid(letter ~.,scales="free_y",switch="y") + 
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0, "lines"), strip.text.y.left = element_text(angle = 0,color="white"), axis.title=element_text(size=18))  +
  xlab(label="LC/MS Peak") + 
  ylab(label="Young leaves")

popcolor(heat_y)

setEPS()
postscript("FiguresTables/HeatMapYoungImportant.eps", height = 5, width = 6)
popcolor(heat_y + theme(legend.position = "none"))
dev.off()









########## example code to make colored boxes around the site letters
#I found this example here https://github.com/tidyverse/ggplot2/issues/2096 and Dan Turner helpfully commented it for me:
p <- ggplot(mpg, aes(displ, cty)) + geom_point() + facet_grid(drv ~ cyl) # create an object of a ggplot
p # print the object
g <- ggplot_gtable(ggplot_build(p)) # plot the constituent parts of the 'grid graphic object' (aka grob)
strip_both <- which(grepl('strip-', g$layout$name)) # find the index in g$layout$name for each of the facet "strips" from both the top and right facet names (aka background rectangles) and put them into a vector to be iterated later on in the for loop. Check g$layout$name to see the whole list of stuff that can be modified
fills <- c("red","green","blue","yellow","red","green","blue","yellow") # make a vector of color names (here the example used more than necessary, which is okay, as long as there are seven for the plot to cycle through)
k <- 1 # start index at the first observation/cell
for (i in strip_both) { # initiate for loop for how ever many places there are to change the colors
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder)) # locate in this huge gtable object where the rectangles are to replace
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k] # replace those names with the fills given in that vector above
  k <- k+1 # add one to the index to move onto the next strip value
}
grid::grid.draw(g)
rm(g,i,j,k,p,fills,strip_both)

#Dan provided this modification to the strip vector to change the color on only one side.
p <- ggplot(mpg, aes(displ, cty)) + geom_point() + facet_grid(drv ~ cyl) # create an object of a ggplot
p # print the object
g <- ggplot_gtable(ggplot_build(p)) # plot the constituent parts of the 'grid graphic object' (aka grob)
strip_r <- which(grepl('strip-r', g$layout$name)) # this only replaces the colors on facet backgrounds at the top of the plot
fills <- c("red","green","blue","yellow","red","green","blue","yellow") # make a vector of color names (here the example used more than necessary, which is okay, as long as there are seven for the plot to cycle through)
k <- 1 # start index at the first observation/cell
for (i in strip_r) { # initiate for loop for how ever many places there are to change the colors
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder)) # locate in this huge gtable object where the rectangles are to replace
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k] # replace those names with the fills given in that vector above
  k <- k+1 # add one to the index to move onto the next strip value
}
grid::grid.draw(g)
rm(g,i,j,k,p,fills,strip_t)

#troubleshoot for my plot. Note that strip.background = element_blank() needs to be removed from theme() in the plot for this to work!
g <- ggplot_gtable(ggplot_build(heat_y))
strip_side <- which(grepl('strip-l', g$layout$name))
fills <- c("red","green","blue","yellow","orange","purple","pink","white","red","green","blue","yellow","orange","purple","pink","white")
k <- 1
for (i in strip_side) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g)
rm(g,i,j,k,p,fills,strip_side)

###################
#old version of all chemicals that uses nested faceting

#this prepares for the geom_tile. For troubleshooting, use a subset [1:100}:
#chem_melt <- melt(chem[1:100])
chem_melt <- melt(chem)
head(chem_melt)
chem_melt$age <- as.factor(chem_melt$age) %>% recode(young = "young leaves", mature = "mature leaves")

heat <- ggplot(chem_melt, aes(x=variable,y=pos_age)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_viridis(direction=-1,limits=c(0,10.7)) + 
  facet_nested(age + letter ~.,scales="free_y",switch="y",nest_line=T) + 
  theme_bw() +theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.spacing = unit(.1, "lines"), strip.background = element_blank())  +
  xlab(label="LC/MS Peak") + 
  ylab(label="Site")

heat

#putting this in theme() rotates both the age and site, and the rotated letters have a lot of space around them
#, strip.text.y.left = element_text(angle = 0)