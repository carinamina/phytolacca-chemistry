# results tables and figures that combine info from young and mature leaves
library(tidyverse)
library(ggplot2)
#################
# compare fit of models with and without population
#################

rsq <- data.frame(matrix(0, nrow = 6, ncol = 4))
names(rsq) <- c("model","age","R2","SE")
rsq$model <- rep(c("population","traits","traits + population"),2)
rsq$age <- c(rep("mature",3),rep("young",3))

#functions to extract R2 and SE from RF results
res <- function(results_file)
{
  results <- read.csv(results_file, sep = "\t", header=T)
  (c(as.numeric(results[24,1]),as.numeric(results[26,1])))
}

rsq[1,3:4] = res("RF_R/mature_pop_results.txt")
rsq[2,3:4] = res("RF_R/maturefine_9_results.txt")
rsq[3,3:4] = res("RF_R/mature_popchem_results.txt")
rsq[4,3:4] = res("RF_R/young_pop_results.txt")
rsq[5,3:4] = res("RF_R/youngfine_24_results.txt")
rsq[6,3:4] = res("RF_R/young_popchem_results.txt")

#plot R^2
rsq$model <- gsub(" ", "\n", rsq$model)
plot = ggplot(rsq,aes(x=model, y=R2, group = sort(age, decreasing = TRUE))) + geom_bar(stat = "identity", position = position_dodge(width =.9), aes(fill = age)) + scale_fill_manual(values=c("gray28","gray87")) + geom_errorbar(position = position_dodge(width =.9), aes(ymin=R2-SE, ymax=R2+SE),width=.2,size=.4,color="black") 
vv=plot + xlab("Factors in model")+ylab(expression("Model "~R^2)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 24)) + scale_x_discrete(limits=c("population","traits","traits\n+\npopulation"))
vv
setEPS()
postscript("FiguresTables/Fig_RF_pop.eps")
vv + theme(legend.position = "none")
dev.off()
#this version has legend
setEPS()
postscript("FiguresTables/Fig_RF_pop_legend.eps")
vv + labs(fill="Leaf age")
dev.off()

rm(rsq,res,vv,plot)
#################
# tables with important compounds: importance scores
#################
imp.y <- read.csv("RF_R/young_48_imp", sep = "\t", header=T)
sum(imp.y$mean_imp)
names(imp.y) <- c("compound","imp.y")

imp.m <- read.csv("RF_R/mature_21_imp", sep = "\t", header=T)
sum(imp.m$mean_imp)
names(imp.m) <- c("compound","imp.m")

imp <- full_join(imp.y,imp.m,by="compound")

#reshape table from wide to long to plot young and mature points separately
imp.wide <- na.omit(reshape(imp,
                    varying = c("imp.m", "imp.y"),
                    v.names = "importance",
                    timevar = "age",
                    times = c("mature","young"),
                    direction = "long"))
imp.wide$id <- NULL
#break off latitude because I want to treat retention time as a numeric
imp_lat <- imp.wide[(imp.wide$compound == "lat"),]
imp.wide <- imp.wide[(imp.wide$compound != "lat"),]
# get retention time from the compound name
imp.wide$time <- as.numeric(str_remove(unlist(strsplit(imp.wide$compound,"_"))[seq(1:nrow(imp.wide))*2-1],"[X]"))

#plot importance scores by retention time
plot = ggplot(imp.wide,aes(x=time, y=importance, shape = age)) + geom_point(size=3, stroke = 1) + scale_shape_manual(values = c(19,1)) + geom_abline(slope=0, intercept= imp_lat[2,3], linetype = "dashed") + geom_abline(slope=0, intercept= imp_lat[1,3])
vv=plot + xlab("Retention time")+ylab("Importance score") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 24)) #+ coord_cartesian(xlim=c(0,7.5))
vv
# setEPS()
# postscript("FiguresTables/Fig_RFresults.eps")
# vv + theme(legend.position = "none")
# dev.off()
# setEPS()
# postscript("FiguresTables/Fig_RFresults_legend.eps")
# vv + labs(shape="Leaf age")
# dev.off()

rm(plot,vv,imp,imp_lat,imp.wide)

#################
# correlations of important compounds with palatability and latitude
#################
#tricky part is keeping the 71 lines used for young and 74 for mature because na.omit() on the original datasheet takes out too many and on the selected columns takes out too few.
indiv.y <- read.csv("RF_R/RF_young_tab.csv",sep="\t",header=T) %>% select(line) %>% mutate(line_age = paste(line,"young",sep="_"))
indiv.m <- read.csv("RF_R/RF_mature_tab.csv",sep="\t",header=T) %>% select(line) %>% mutate(line_age= paste(line,"mature",sep="_"))
indiv <- c(indiv.y$line_age,indiv.m$line_age)

traits <- read.csv("Processing/2_out_AllTraits.csv",header=T) %>% mutate(line_age = paste(line,age,sep="_")) %>% filter(line_age %in% indiv) %>% select(c(pop,age,lat,area,unique(c(imp.y$compound,imp.m$compound)))) 
rm(indiv.y,indiv.m,indiv)

#use population means for correlations so we don't pseudoreplicate. separate into young and mature.
pop.y <- subset(aggregate(. ~ pop + age, data=traits, mean), age == "young")
pop.m <- subset(aggregate(. ~ pop + age, data=traits, mean), age == "mature")

#create table with correlation coefficients for each compound: young
cor.y <- select(imp.y,compound) %>% filter(compound != "lat")
cor.y$age <- "young"

#testing: this snippet gets abundance from pop.y of the compound found in row 13 of cor.y (which happens to be column 17 in pop.y). compare:
#pop.y[,colnames(pop.y)==cor.y[13,1]]
#cor.test(pop.y$lat, pop.y[,17], method="pearson")
#cor.test(pop.y$lat, pop.y[,colnames(pop.y)==cor.y[13,1]], method="pearson")
#cor.test(pop.y$lat, pop.y[,13],method="spearman",exact=F)$p.value
#cor.test(pop.y$lat, pop.y[,13],method="spearman",exact=F)$estimate

#i know this is very inefficient, oh well. loops through and gets correlation coefficients, p-values, for both pearson and spearman. also 95% c.i. for pearson.
for(i in 1:nrow(cor.y)){
    #run spearman correlation with latitude for each compound
  spearman <- cor.test(pop.y$lat, pop.y[,colnames(pop.y)==cor.y[i,1]], method="spearman",exact=F)
  #returns correlation slope
  cor.y$lat.rho[i] <- as.numeric(spearman$estimate)
  #returns correlation p-value
  cor.y$lat.p[i] <- as.numeric(spearman$p.value)
  cor.y$lat.sig[i] <- ifelse(spearman$p.value < 0.05/(nrow(imp.y)-1), "yes", "no")
  
  #run pearson correlation with palatability for each compound
  pearson <- cor.test(pop.y$area, pop.y[,colnames(pop.y)==cor.y[i,1]], method="pearson")
  #returns correlation slope
  cor.y$palat.r[i] <- as.numeric(pearson$estimate)
  #returns correlation p-value
  cor.y$palat.p[i] <- as.numeric(pearson$p.value)
  cor.y$palat.sig[i] <- ifelse(pearson$p.value < 0.05/(nrow(imp.y)-1), "yes", "no")
  
}
rm(i,pearson,spearman)

#create table with correlation coefficients for each compound: young
cor.m <- select(imp.m,compound) %>% filter(compound != "lat")
cor.m$age <- "mature"

#i know this is very inefficient, oh well. loops through and gets correlation coefficients, p-values, for both pearson and spearman. also 95% c.i. for pearson.
for(i in 1:nrow(cor.m)){
  #run spearman correlation with latitude for each compound
  spearman <- cor.test(pop.m$lat, pop.m[,colnames(pop.m)==cor.m[i,1]], method="spearman",exact=F)
  #returns correlation slope
  cor.m$lat.rho[i] <- as.numeric(spearman$estimate)
  #returns correlation p-value
  cor.m$lat.p[i] <- as.numeric(spearman$p.value)
  cor.m$lat.sig[i] <- ifelse(spearman$p.value < 0.05/(nrow(imp.m)-1), "yes", "no")
  
  #run pearson correlation with palatability for each compound
  pearson <- cor.test(pop.m$area, pop.m[,colnames(pop.m)==cor.m[i,1]], method="pearson")
  #returns correlation slope
  cor.m$palat.r[i] <- as.numeric(pearson$estimate)
  #returns correlation p-value
  cor.m$palat.p[i] <- as.numeric(pearson$p.value)
  cor.m$palat.sig[i] <- ifelse(pearson$p.value < 0.05/(nrow(imp.m)-1), "yes", "no")
  
}
rm(i,pearson,spearman)

#plot correlation coefficients by retention time against lat and palat
cor <- bind_rows(cor.m,cor.y)
rm(cor.y,cor.m)
cor$time <- as.numeric(str_remove(unlist(strsplit(cor$compound,"_"))[seq(1:nrow(cor))*2-1],"[X]"))

cor$lat.sign <- ifelse(cor$lat.rho>0, "+", "-")
cor$palat.sign <- ifelse(cor$palat.r>0, "+", "-")
count(cor, age, lat.sig)
count(cor,age,lat.sig,lat.sign)
count(cor, age, palat.sig)
count(cor, age, palat.sig,palat.sign)
cor$bothyes <- ifelse(cor$lat.sig == "yes" & cor$palat.sig == "yes","yes","no" )
count(cor,age,bothyes)
count(cor,age,bothyes,lat.sign)


#plot correlations with latitude
plot = ggplot(cor,aes(x=time, y=lat.rho, shape = age, color = lat.sig)) + geom_point(size=3, stroke = 1) + scale_shape_manual(values = c(19,1)) + scale_color_manual(values = c("grey","black")) + geom_abline(slope=0, intercept= 0)
vv=plot + xlab("Retention time")+ylab("Spearman corr coeff with lat") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 24))
vv
# setEPS()
# postscript("FiguresTables/Fig_RFchemlat.eps")
# vv + theme(legend.position = "none")
# dev.off()

#plot correlations with palatability
plot = ggplot(cor,aes(x=time, y=palat.r, shape = age, color = palat.sig)) + geom_point(size=3, stroke = 1) + scale_shape_manual(values = c(19,1)) + scale_color_manual(values = c("grey","black")) + geom_abline(slope=0, intercept= 0)
vv=plot + xlab("Retention time")+ylab("Pearson corr coeff with palatability") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 24))
vv
# setEPS()
# postscript("FiguresTables/Fig_RFchempalat.eps")
# vv + theme(legend.position = "none")
# dev.off()

rm(list=ls())
