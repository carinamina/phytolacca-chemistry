library(ecodist)
library(vegan)
library(ade4)
library(tidyverse)
library(ggplot2)

#makes scree plots (takes a few minutes)
NMDS.scree<-function(x,name) #where x is the name of the data frame variable
{ 
  plot(rep(1,10),replicate(10,metaMDS(x,autotransform=F,k=1,trymax=500)$stress),xlim=c(1,11),ylim=c(0,0.3),xlab="# of Dimensions",ylab="Stress",main=name)
  
  for (i in 1:10) {
    points(rep(i+1,10),replicate(10,metaMDS(x,autotransform=F,k=i+1)$stress))
  }
}

##Data import for NMDS: get chemicals, subset to young and mature

chem <- read.csv("Processing/1a_out_LCMSprocessed_CB.csv", header = TRUE) %>% mutate(line_age = paste(paste(pop,line,sep="_"),age,sep="_"),.before=pos)
#use pop_line_age as unique identifier; they are not repeated
length(unique(chem$line_age))
levels(chem$region)<-c(levels(chem$region),"n/a")
chem$region[is.na(chem$region)] <- "n/a"    
str(chem[,1:11])

##Both ages, all compounds
chem.mat <- data.matrix(chem[,11:length(chem)])
#use line_age for identifier
rownames(chem.mat)<- chem$line_age
key <- chem[,1:10]

# make bray curtis distance matrix from data
dist <-distance(chem.mat,"bray-curtis")

#scree plots
# setEPS()
# postscript("FiguresTables/20210308_nmds_scree.eps")
# NMDS.scree(dist, name = "young and mature leaves")
# dev.off()
#this says you want stress to be 0.05-1 https://mb3is.megx.net/gustame/dissimilarity-based-methods/nmds
#same here https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
#4 dimensions is less than 0.1 so I'll set k=4 in metaMDS call

nmds<-metaMDS(chem.mat,k=4,trymax=250,distance="bray")
stressplot(nmds)
#Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions. https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
#looks pretty good
nmds$stress
#values < 0.05 excellent, < 0.1 good; this is ok at 0.09 (if k = 3, it's .12 and the plot looks slightly messier)

#this extracts scores to a matrix
scores <- as.data.frame(scores(nmds))
scores$line_age <- rownames(scores)  # create a column of rownames (line_age)
#use key to get latitude based on pos
scores <- merge(key, scores, by = "line_age")
head(scores)

write.csv(scores, file = "Processing/1b_out_NMDSscores.csv", row.names = FALSE )

#################################
#Plots

#colored by latitude, shaped by age
plot = ggplot(scores,aes(x=NMDS1, y=NMDS2, colour = lat, shape = age))+  geom_point(size=4, stroke = 2) + scale_colour_gradient2(low="gold",mid = "red", high="blue", midpoint = 26) + scale_shape_manual(values = c(19,1))
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 24))
vv
# setEPS()
# postscript("nmds1_all.eps")
# vv + theme(legend.position = "none")
# dev.off()
# setEPS()
# postscript("20180517_nmds1_all_legend.eps")
# vv + labs(shape="Leaf age", color = "Latitude")
# dev.off()
plot = ggplot(scores,aes(x=NMDS3, y=NMDS4, colour = lat, shape = age))+  geom_point(size=4, stroke = 2) + scale_colour_gradient2(low="gold",mid = "red", high="blue", midpoint = 26) + scale_shape_manual(values = c(19,1))
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 24))
vv
# setEPS()
# postscript("nmds3_all.eps")
# vv + theme(legend.position = "none")
# dev.off()

#re-order levels for key, also makes listing colors easier
scores$region <- factor(scores$region, levels =c("temperate", "n/a","subtropical", "tropical") )
levels(scores$region)

#colored by region
plot = ggplot(scores,aes(x=NMDS1, y=NMDS2, colour = region, shape = age))+  geom_point(size=4, stroke = 2) + scale_colour_manual(values=c("steelblue1","grey","navyblue","maroon2"))+ scale_shape_manual(values = c(19,1))
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 20))
vv
# setEPS()
# postscript("nmds1_all.eps")
# vv + theme(legend.position = "none")
# dev.off()
# setEPS()
# postscript("nmds1_all_legend.eps")
# vv
# dev.off()
#axis 1 is mostly separating by latitude (but also age for tropical); for mature leaves, it's tropical/subtropical vs everything north of FL; for young leaves, it's similar but tropical and subtropical are separating out more. axis 2 separates by leaf age for PHAM but not PHRI. THIS IS SO COOL because it shows that more divergence is happening in young leaves! And it's crazy that FL is grouping with CR!
plot = ggplot(scores,aes(x=NMDS3, y=NMDS4, colour = region, shape = age))+  geom_point(size=4, stroke = 2) + scale_colour_manual(values=c("steelblue1","grey","navyblue","maroon2"))+ scale_shape_manual(values = c(19,1))
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 20))
vv
# setEPS()
# postscript("20180531_nmds2_all.eps")
# vv + theme(legend.position = "none")
# dev.off()
#for the most part leaf ages still cluster geographically, but there is more geographic overlap in general. Axis 3 is separating young from mature, to differing degrees depending on region. Axis 4 separates PHAM from PHRI (more or less)

rm(list=ls())
