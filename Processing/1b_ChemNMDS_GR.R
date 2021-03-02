##---
##title: "NMDS"
##author: "Carina Baskett"
##date: "5/13/2018"
##output: html_document
##---

##```{r setup, include=FALSE}
##knitr::opts_chunk$set(echo = TRUE)
library(ecodist)
library(vegan)
library(ade4)
library(plyr)
library(ggplot2)

#makes scree plots (takes a few minutes)
NMDS.scree<-function(x,name) #where x is the name of the data frame variable
{ 
  plot(rep(1,10),replicate(10,metaMDS(x,autotransform=F,k=1,trymax=500)$stress),xlim=c(1,11),ylim=c(0,0.3),xlab="# of Dimensions",ylab="Stress",main=name)
  
  for (i in 1:10) {
    points(rep(i+1,10),replicate(10,metaMDS(x,autotransform=F,k=i+1)$stress))
  }
}

##```

##Data import for NMDS: get chemicals, subset to young and mature, further subset to lists of important compounds we got from RF
##```{r}

#setwd("/Users/carina/Documents/R_working_directory")
chem <- read.csv("data_out/20180507_chem_abundant1percent.csv", header = TRUE)
chem$X <- NULL
levels(chem$region)<-c(levels(chem$region),"n/a")
chem$region[is.na(chem$region)] <- "n/a"    
str(chem[,1:14])

chem.y <- droplevels(subset(chem, chem$age == "young"))
chem.m <- droplevels(subset(chem, chem$age == "mature"))

#import lists of important compounds
imp.y <- read.csv("data_out/20180517_imp_young.csv", header = TRUE)
imp.list.y <- as.vector(imp.y$V1)
imp.m <- read.csv("data_out/20180517_imp_mature.csv", header = TRUE)
imp.list.m <- as.vector(imp.m$V1)

#subset by those compounds
chem.y.RF <- data.frame(chem.y[,1:13],chem.y[,colnames(chem.y)%in%imp.list.y])
chem.m.RF <- data.frame(chem.m[,1:13],chem.m[,colnames(chem.m)%in%imp.list.m])
rm(imp.m, imp.y, imp.list.m, imp.list.y,chem.m,chem.y)

##```

##Both ages, all compounds
##```{r}
chem.mat <- data.matrix(chem[,14:length(chem)])
rownames(chem.mat)<- chem$line_age
key <- chem[c("line_age","lat","region","age")]

# make bray curtis distance matrix from data
dist <-distance(chem.mat,"bray-curtis") 
#Following was in Marge's code but I'm not sure I need any of it
# range(distance.y)
# is.euclid(distance.y)
# is.euclid(sqrt(distance.y))
# distance.y<-sqrt(distance.y)

#scree plots
# setEPS()
# postscript("20180517_nmds_scree.eps")
# NMDS.scree(dist, name = "young and mature leaves")
# dev.off()
#4 dimensions seems good, so I'll set k=4 in metaMDS call

nmds<-metaMDS(chem.mat,k=4,trymax=250,distance="bray")
stressplot(nmds)
#looks good
nmds$stress
#values < 0.05 excellent, < 0.1 good; this is good at 0.58 (if k = 2, it's .12 and the plot looks messier)

#this extracts scores to a matrix
scores <- as.data.frame(scores(nmds))
scores$line_age <- rownames(scores)  # create a column of rownames (line_age)
#use key to get latitude based on line_age
scores <- merge(key, scores, by = "line_age")
head(scores)

#colored by latitude, shaped by age
plot = ggplot(scores,aes(x=NMDS1, y=NMDS2, colour = lat, shape = age))+  geom_point(size=4, stroke = 2) + scale_colour_gradient2(low="gold",mid = "red", high="blue", midpoint = 26) + scale_shape_manual(values = c(19,1))
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 24))
vv
# setEPS()
# postscript("20180517_nmds1_all.eps")
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
# postscript("20180517_nmds3_all.eps")
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
# postscript("20180531_nmds1_all.eps")
# vv + theme(legend.position = "none")
# dev.off()
# setEPS()
# postscript("20180531_nmds1_all_legend.eps")
# vv
# dev.off()
#axis 1 is mostly separating by latitude; for mature leaves, it's tropical/subtropical vs everything north of FL; for young leaves, it's similar but tropical and subtropical are separating out more. axis 2 separates by leaf age, with more separation for tropical and subtropical. THIS IS SO COOL because it shows that more divergence is happening in young leaves! And it's crazy that FL is grouping with CR!
plot = ggplot(scores,aes(x=NMDS3, y=NMDS4, colour = region, shape = age))+  geom_point(size=4, stroke = 2) + scale_colour_manual(values=c("steelblue1","grey","navyblue","maroon2"))+ scale_shape_manual(values = c(19,1))
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 20))
vv
# setEPS()
# postscript("20180531_nmds2_all.eps")
# vv + theme(legend.position = "none")
# dev.off()
#for the most part leaf ages still cluster geographically, but there is more geographic overlap in general. Axis 3 seems to separate PHAM young from mature, but not so much for PHRI. Axis 4 separates geography, but interestingly northern PHAM are now between FL and tropical. Wild speculation here, but imagine back to a common ancestor of all 3...FL and PHRI have diverged a lot from each other, while north PHAM hasn't evolved a whole lot along this axis. Very cool.

#again, using Jaccard (presence/absence)
#nmds2<-metaMDS(chem.mat,k=4,trymax=250,distance="jaccard", binary = TRUE)

write.csv(scores, file = "20180531_nmds_scores.csv", row.names = FALSE )
##```

##Young leaves, important compounds
##```{r}
rm(plot, vv, chem.mat,dist,key,nmds,scores, chem)
chem.y.RF.mat <- data.matrix(chem.y.RF[,14:length(chem.y.RF)])
rownames(chem.y.RF.mat)<- chem.y.RF$line_age
key.y.RF <- chem.y.RF[c("line_age","lat","region")]

# make bray curtis distance matrix from data
distance.y.RF<-distance(chem.y.RF.mat,"bray-curtis") 
range(distance.y.RF)
#Following was in Marge's code but I'm not sure I need any of it
# is.euclid(distance.y.RF)
# is.euclid(sqrt(distance.y.RF))
# distance.y.RF<-sqrt(distance.y.RF)

#scree plots
# setEPS()
# postscript("20180517_nmds_scree.y.RF.eps")
# NMDS.scree(distance.y.RF,name="young")
# dev.off()
#2 dimensions gets to about 0.5, so I'll set k=2 in metaMDS call

nmds.y.RF<-metaMDS(chem.y.RF.mat,k=2,trymax=250,distance="bray")
stressplot(nmds.y.RF)
#looks good
nmds.y.RF$stress
#values < 0.05 excellent, < 0.1 good; this is good

#this extracts scores to a matrix
scores.y.RF <- as.data.frame(scores(nmds.y.RF))
scores.y.RF$line_age <- rownames(scores.y.RF)  # create a column of rownames (line_age)
#use key to get latitude based on line_age
scores.y.RF <- merge(key.y.RF, scores.y.RF, by = "line_age")
head(scores.y.RF)

plot = ggplot(scores.y.RF,aes(x=NMDS1, y=NMDS2, colour = lat))+  geom_point(size=4, stroke = 2, shape = 1) + scale_colour_gradient2(low="gold",mid = "red", high="blue", midpoint = 26)
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 24))
vv
# setEPS()
# postscript("20180517_nmds1.y.RF.eps")
# vv + theme(legend.position = "none")
# dev.off()

#re-order levels for key, also makes listing colors easier
scores.y.RF$region <- factor(scores.y.RF$region, levels =c("temperate", "n/a","subtropical", "tropical") )
levels(scores.y.RF$region)

#colored by region
plot = ggplot(scores.y.RF,aes(x=NMDS1, y=NMDS2, colour = region))+  geom_point(size=4, stroke=2,shape=1) + scale_colour_manual(values=c("steelblue1","grey","navyblue","maroon2"))
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 20))
vv
# setEPS()
# postscript("20180531_nmds1.y.RF.eps")
# vv + theme(legend.position = "none")
# dev.off()

##```

##Mature leaves, important compounds
##```{r}
rm(plot, vv, chem.y.RF.mat,distance.y.RF,key.y.RF,chem.y.RF,scores.y.RF)
chem.m.RF.mat <- data.matrix(chem.m.RF[,14:length(chem.m.RF)])
rownames(chem.m.RF.mat)<- chem.m.RF$line_age
key.m.RF <- chem.m.RF[c("line_age","lat","region")]

#Following was in Marge's code but I'm not sure I need any of it
# make bray curtis distance matrix from data
distance.m.RF<-distance(chem.m.RF.mat,"bray-curtis") 
# range(distance.m.RF)
# is.euclid(distance.m.RF)
# is.euclid(sqrt(distance.m.RF))
# distance.m.RF<-sqrt(distance.m.RF)

#scree plots
# NMDS.scree<-function(x,name) {#where x is the name of the data frame variable
#   
#   plot(rep(1,10),replicate(10,metaMDS(x,autotransform=F,k=1,trymax=500)$stress),xlim=c(1,4),ylim=c(0,0.3),xlab="# of Dimensions",ylab="Stress",main=name)
#   
#   for (i in 1:4) {
#     points(rep(i+1,10),replicate(10,metaMDS(x,autotransform=F,k=i+1)$stress))
#   }
# }
# setEPS()
# postscript("20180517_nmds_scree.m.RF.eps")
# NMDS.scree(distance.m.RF,name="mature")
# dev.off()
#2 dimensions seems good, so I'll set k=2 in metaMDS call

nmds.m.RF<-metaMDS(chem.m.RF.mat,k=2,trymax=250,distance="bray")
stressplot(nmds.m.RF)
#looks good
nmds.m.RF$stress
#values < 0.05 excellent, < 0.1 good; this is good
(nmds.m.RF$stress + nmds.y.RF$stress)/2

#this extracts scores to a matrix
scores.m.RF <- as.data.frame(scores(nmds.m.RF))
scores.m.RF$line_age <- rownames(scores.m.RF)  # create a column of rownames (line_age)
#use key to get latitude based on line_age
scores.m.RF <- merge(key.m.RF, scores.m.RF, by = "line_age")
head(scores.m.RF)

plot = ggplot(scores.m.RF,aes(x=NMDS1, y=NMDS2, colour = lat))+  geom_point(size=4, stroke=2,shape=19) + scale_colour_gradient2(low="gold",mid = "red", high="blue", midpoint = 26)
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 24))
vv
# setEPS()
# postscript("20180517_nmds1_all.m.RF.eps")
# vv + theme(legend.position = "none")
# dev.off()

#re-order levels for key, also makes listing colors easier
scores.m.RF$region <- factor(scores.m.RF$region, levels =c("temperate", "n/a","subtropical", "tropical") )
levels(scores.m.RF$region)

#colored by region
plot = ggplot(scores.m.RF,aes(x=NMDS1, y=NMDS2, colour = region))+  geom_point(size=4, stroke=2) + scale_colour_manual(values=c("steelblue1","grey","navyblue","maroon2"))
vv=plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 20))
vv
# setEPS()
# postscript("20180531_nmds1.m.RF.eps")
# vv + theme(legend.position = "none")
# dev.off()

rm(plot, vv, chem.m.RF.mat,distance.m.RF,key.m.RF,nmds.m.RF,scores.m.RF,nmds.y.RF,chem.m.RF)
rm(NMDS.scree)
##```
