#summarizing LC/MS data into NMDS axes, chemical abundance, chemical richness, and chemical diversity

library(ecodist)
library(vegan)
library(ade4)
library(tidyverse)
library(ggplot2)

##Data import for NMDS: get chemicals, subset to young and mature

chem <- read.csv("Processing/1a_out_LCMS_Indiv.csv", header = TRUE)
#use pop_line_age as unique identifier; they are not repeated
length(unique(chem$line_age))
str(chem[,1:9])

##Both ages, all compounds
chem.mat <- data.matrix(chem[,9:length(chem)])
#use line_age for identifier
rownames(chem.mat)<- chem$line_age
key <- chem[,1:8]

#from chemistry matrix, calculate Shannon diversity, richness, and abundance for each sample
chem.sums <- data.frame(row.names = row.names(chem.mat))
#abundance is sum of all compound abundances. We really want log.abund, which is beautifully normally distributed
chem.sums$abund <- rowSums(chem.mat)
hist(chem.sums$abund)
chem.sums$log.abund <- log(chem.sums$abund)
hist(chem.sums$log.abund)
#richness is number of compounds present. counts the number of zeroes in each row and subtracts from total to calculate richness
chem.sums$richness <- ncol(chem.mat) - rowSums(chem.mat == 0)
hist(chem.sums$richness)
#test code for richness calculation
# abc <- chem.mat[1:17,1:4]
# rowSums(abc== 0)
#Shannon diversity index
chem.sums$diversity <- diversity(chem.mat)
hist(chem.sums$diversity)


#NMDS
# make bray curtis distance matrix from data
dist <-distance(chem.mat,"bray-curtis")

#function to make scree plot
NMDS.scree<-function(x,name) #where x is the name of the data frame variable
{ 
  plot(rep(1,10),replicate(10,metaMDS(x,autotransform=F,k=1,trymax=250)$stress),xlim=c(1,11),ylim=c(0,0.3),xlab="# of Dimensions",ylab="Stress",main=name)
  
  for (i in 1:10) {
    points(rep(i+1,10),replicate(10,metaMDS(x,autotransform=F,k=i+1)$stress))
  }
}

#scree plot (takes a few minutes)
#I tried running this with the chem.mat instead of the distance matrix, because it seems like the metaMDS makes the distance matrix, and the plot looked exactly the same. I don't understand why stress is around 0.15 in this plot for k=2, yet it's 0.18 when I calculate it below. Is this plot really necessary??
# setEPS()
# postscript("FiguresTables/FigS1_nmds_scree.eps")
# NMDS.scree(dist, name = "young and mature leaves")
# dev.off()
#this says you want stress to be 0.05-1 https://mb3is.megx.net/gustame/dissimilarity-based-methods/nmds
#same here https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
#2 dimensions is less than 0.2 so I'll set k=2 in metaMDS call (see exact value in a few lines)

nmds<-metaMDS(chem.mat,k=2,trymax=250,distance="bray")
stressplot(nmds)
#Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions. https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
#looks pretty good
nmds$stress
#values < 0.05 excellent, < 0.1 great, < 0.2 is good/ok; this is ok at 0.18. I think since I'm pretty much just using this for visualization it's better to stick with two axes--way easier to visualize!

#this extracts scores to a matrix
scores <- as.data.frame(scores(nmds))

#use key to get other info based on line_age and merge chem.sums with nmds scores
scores$line_age <- rownames(scores)  # create a column of rownames (line_age)
chem.sums$line_age <- rownames(chem.sums)
scores <- merge(chem.sums,scores, by = "line_age")
scores <- merge(key, scores, by = "line_age")
head(scores)
scores$abund <- NULL

write.csv(scores[c(2,1,3:ncol(scores))], file = "Processing/1b_out_ChemSummaries_Indiv.csv", row.names = FALSE )

#For each NMDS axis, add a variable I'll call "ontogenetic similarity" or onto.1, onto.2 etc which is the difference in NMDS scores for young and mature leaves on the same plant (absolute value, because NMDS axes are arbitrary)
length(unique(scores$pos))
plant <- scores %>% select(pos,age,NMDS1,NMDS2) %>% reshape(timevar = "age", idvar = "pos", direction = "wide") %>% mutate(onto.1 = abs(NMDS1.mature - NMDS1.young), onto.2 = abs(NMDS2.mature - NMDS2.young)) %>% right_join(distinct(select(scores,c(pos,pop,line,lat,region))),by="pos") %>% select(pos, pop, line,lat,region, onto.1, onto.2)

write.csv(plant, "Processing/1b_out_OntoSimilarity_Indiv.csv", row.names = FALSE)


#some visualization for fun
qplot(lat, log.abund, data=scores, colour=age)
qplot(lat, onto.2, data=plant)

rm(list=ls())
