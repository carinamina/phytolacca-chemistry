#testing how cutoffs in chemical abundance affect regional comparisons of abundance, richness, and diversity

####################
#importing data with different cutoffs (1%, 0.5%, 0.1% relative abundance) and get shannon diversity, total abundance, and richness
####################
library(nlme)
library(corrplot)
library(Hmisc)
library(lsmeans)
lme_results <- function(lme_model)
{
  a<- hist((resid(lme_model) - mean(resid(lme_model), na.rm=T)) / sd(resid(lme_model), na.rm=T), freq=F); curve(dnorm, add = TRUE)
  b<-qqnorm(lme_model)
  c<-plot(lme_model)
  d<-summary(lme_model)
  #  e<-intervals(lme_model)[1]
  f<-anova(lme_model, type = "marginal", test = "F")
  g<-nobs(lme_model)
  
  list(a,b,c,d,f,g)
}


library(vegan)
#first, the 1% cutoff
chem <- read.csv("20180531_chem_abundant1percent.csv", header = TRUE)
chem$X <- NULL
str(chem[,1:14])

chem.mat <- data.matrix(chem[,14:length(chem)])
chem_sums <- chem[c("line_age")]
chem_sums$shann1 <- diversity(chem.mat)  #Shannon index
hist(chem_sums$shann1)
chem_sums$richness1 <- ncol(chem.mat) - rowSums(chem.mat == 0)
hist(chem_sums$richness1)
chem_sums$abund1 <- rowSums(chem.mat)
hist(chem_sums$abund1)
rm(chem, chem.mat)

#now 0.5% cutoff
chem <- read.csv("20180531_chem_abundant.5percent.csv", header = TRUE)
chem$X <- NULL
str(chem[,1:14])

chem.mat <- data.matrix(chem[,14:length(chem)])
chem_sums$shann.5 <- diversity(chem.mat)  #Shannon index
hist(chem_sums$shann.5)
chem_sums$richness.5 <- ncol(chem.mat) - rowSums(chem.mat == 0)
hist(chem_sums$richness.5)
chem_sums$abund.5 <- rowSums(chem.mat)
hist(chem_sums$abund.5)
rm(chem, chem.mat)

#now 0.1% cutoff
chem <- read.csv("20180531_chem_abundant.1percent.csv", header = TRUE)
chem$X <- NULL
str(chem[,1:14])

chem.mat <- data.matrix(chem[,14:length(chem)])
chem_sums$shann.1 <- diversity(chem.mat)  #Shannon index
hist(chem_sums$shann.1)
chem_sums$richness.1 <- ncol(chem.mat) - rowSums(chem.mat == 0)
hist(chem_sums$richness.1)
chem_sums$abund.1 <- rowSums(chem.mat)
hist(chem_sums$abund.1)
rm(chem, chem.mat)

#get the population, line, leaf age, and region in a key
chem_key <- read.table(text = as.character(chem_sums$line_age), sep = "_", colClasses = "factor")
colnames(chem_key) <- c("pop","line","age")
chem_key$line <- as.factor(paste(chem_key$pop, chem_key$line, sep="_"))
chem_key$line_age <- as.factor(paste(chem_key$line, chem_key$age, sep = "_"))
lats = read.csv("lats_long_names.csv", header=TRUE)
chem_key <- merge(chem_key, lats, by = "pop")
regions <- c("tropical", "tropical","tropical", "subtropical", "subtropical", "subtropical",NA,NA,NA,NA,NA,NA,NA, "temperate", "temperate", "temperate")
regions <- cbind(regions, sort(unique(chem_key$lat)))
colnames(regions) <- c("region", "lat")
chem_key <- merge(chem_key, regions, by = "lat")

chem <- merge(chem_key, chem_sums, by = "line_age")
rm(chem_key, lats, regions, chem_sums)


(trait_list <- colnames(chem)[7:length(colnames(chem))])

#the df was called all in previous code and let's keep that
all <- droplevels(subset(chem, is.na(chem$region) == FALSE))

####################
#visualization
#####################

#check for heterogeneity of variance in each trait
#young
for (i in 1:length(trait_list)){
  plotdata.y <- subset(all, all$age == "young")
  plotdata.y <- cbind(plotdata.y[, colnames(plotdata.y)%in%trait_list[i]], plotdata.y[3])
  boxplot(plotdata.y[,1]~plotdata.y[,2],
          data=na.omit(plotdata.y),
          main=paste(trait_list[i],"by pop_young"), 
          xlab="population" ,
          ylab = trait_list[i])
}
for (i in 1:length(trait_list)){  
  plotdata.m <- subset(all, all$age == "mature")
  plotdata.m <- cbind(plotdata.m[, colnames(plotdata.m)%in%trait_list[i]], plotdata.m[3])
  boxplot(plotdata.m[,1]~plotdata.m[,2],
          data=na.omit(plotdata.m),
          main=paste(trait_list[i],"by pop_mature"), 
          xlab="population" ,
          ylab = trait_list[i])
}
#richness and abundance are not that variable, but diversity is

#####################
#latitude and age as predictors of trait variation
#####################
#all the traits we want to test
trait_list
par(mfrow=c(1,1))
#### first, check assumptions for all models
#shows the three usual assumption-checking figures for a given trait

all$pop_age <- as.factor(paste(all$pop, all$age, sep = "_"))

ass_check <- function(trait){
  fmla <- as.formula(paste(paste(trait),"~ region*age"))
  lme_model <- lme(fmla, random=~1|pop, method = "REML", data=all, na.action = na.exclude)
  a<- hist((resid(lme_model) - mean(resid(lme_model), na.rm=T)) / sd(resid(lme_model), na.rm=T), freq=F, main = paste(trait,sep="_"),xlab="residuals"); curve(dnorm, add = TRUE)
  b<-qqnorm(lme_model, main = paste(trait,sep="_"))
  c<-plot(lme_model, main = paste(trait,sep="_"))
  list(a,b,c)
  #, weights=varIdent(form=~1|pop_age)
}

ass_check_weight <- function(trait){
  fmla <- as.formula(paste(paste(trait),"~ region*age"))
  lme_model <- lme(fmla, random=~1|pop, method = "REML", data=all, na.action = na.exclude, weights=varIdent(form=~1|pop_age))
  a<- hist((resid(lme_model) - mean(resid(lme_model), na.rm=T)) / sd(resid(lme_model), na.rm=T), freq=F, main = paste(trait,sep="_"),xlab="residuals"); curve(dnorm, add = TRUE)
  b<-qqnorm(lme_model, main = paste(trait,sep="_"))
  c<-plot(lme_model, main = paste(trait,sep="_"))
  list(a,b,c)
}

#again, can't figure out how to condense these figures, so I'll have to do it one at a time

#are weights needed?
ass_check(trait_list[1])
ass_check_weight(trait_list[1])
#no
ass_check(trait_list[2])
ass_check_weight(trait_list[2])
#no
ass_check(trait_list[3])
ass_check_weight(trait_list[3])
#yes; needs pop_age
ass_check(trait_list[4])
ass_check_weight(trait_list[4])
#no
ass_check(trait_list[5])
ass_check_weight(trait_list[5])
#no
ass_check(trait_list[6])
ass_check_weight(trait_list[6])
##yes; needs pop_age
ass_check(trait_list[7])
ass_check_weight(trait_list[7])
#no
ass_check(trait_list[8])
ass_check_weight(trait_list[8])
#no
ass_check(trait_list[9])
ass_check_weight(trait_list[9])
#yes; needs pop_age

#in sum, only the abundance-related measures (3,6,9) need weights (pop_age)

#big tables with results of testing effects of region and age on each trait
###################
lme_table_int <- function(lme_model)
{
  f = anova(lme_model)$F[2]
  numDF = anova(lme_model)$numDF[2]
  denDF = anova(lme_model)$denDF[2]
  p = anova(lme_model)$p[2]
  coef <- as.vector(c(f,numDF,denDF,p))
  names(coef) <- c("F","numDF","denDF","p")
  table <- as.data.frame(t(coef))
  #rownames(table)<-rownames(anova(lme_model))[2]
  return(table)
}

pred_table <- data.frame(  F=numeric(),
                           numDF=numeric(),
                           denDF=numeric(),
                           p=numeric(),
                           trait=character(),
                           stringsAsFactors=FALSE)

#tst <- lme(NMDS1 ~ region*age, random=~1|pop, method = "REML", data=all, na.action = na.exclude)

all$region_age <- interaction(all$region, all$age)
letters <- data.frame(matrix(0, nrow = 6, ncol = 0))

#no weights
trait_list_sub <- trait_list[c(1:2,4:5,7:8)]
for (i in 1:length(trait_list_sub)){
  fmla <- as.formula(paste(paste(trait_list_sub[i]),"~ region*age"))
  mod <- lme(fmla, random=~1|pop, method = "REML", data=all, na.action = na.exclude)
  results <- cbind(as.data.frame(anova(mod))[2:4,],trait_list_sub[i])
  pred_table <- rbind(pred_table,results)
  
  fmla2 <- as.formula(paste(paste(trait_list_sub[i]),"~ region_age"))
  mod2 <- lme(fmla2 , random=~1|pop, method = "REML", data=all, na.action = na.exclude)
  tuk <- glht(mod2, linfct=mcp(region_age="Tukey"))
  summary(tuk)  
  tuk.cld <- as.data.frame(cld(tuk)$mcletters$Letters)
  colnames(tuk.cld) <- trait_list_sub[i]
  letters <- cbind(letters,tuk.cld)
}
#pop_age weights
trait_list_sub <- trait_list[c(3,6,9)]
for (i in 1:length(trait_list_sub)){
  fmla <- as.formula(paste(paste(trait_list_sub[i]),"~ region*age"))
  mod <- lme(fmla, random=~1|pop, method = "REML", data=all, na.action = na.exclude, weights=varIdent(form=~1|pop_age))
  results <- cbind(as.data.frame(anova(mod))[2:4,],trait_list_sub[i])
  pred_table <- rbind(pred_table,results)
  
  fmla2 <- as.formula(paste(paste(trait_list_sub[i]),"~ region_age"))
  mod2 <- lme(fmla2 , random=~1|pop, method = "REML", data=all, na.action = na.exclude, weights=varIdent(form=~1|pop_age))
  tuk <- glht(mod2, linfct=mcp(region_age="Tukey"))
  summary(tuk)  
  tuk.cld <- as.data.frame(cld(tuk)$mcletters$Letters)
  colnames(tuk.cld) <- trait_list_sub[i]
  letters <- cbind(letters,tuk.cld)
}
pred_table
letters

#there are significant region effects on diversity and richness for the 1% and .5%, but not .1%. There are significant age and region*age effects on abundance for all three cutoffs. The abundance letters do not change depending on cutoff.

letters2 <- data.frame(matrix(0, nrow = 3, ncol = 0))
trait_list_sub <- trait_list[c(1,2,4,5)]
for (i in 1:length(trait_list_sub)){
  fmla2 <- as.formula(paste(paste(trait_list_sub[i]),"~ region*age"))
  mod2 <- lme(fmla2 , random=~1|pop, method = "REML", data=all, na.action = na.exclude)
  tuk <- glht(mod2, linfct=mcp(region="Tukey"))
  summary(tuk)  
  tuk.cld <- as.data.frame(cld(tuk)$mcletters$Letters)
  colnames(tuk.cld) <- trait_list_sub[i]
  letters2 <- cbind(letters2,tuk.cld)
}
letters2

###################
#plots of regional means.
###################

plotdata <- data.frame(all[c("age","region","region_age")], all[, colnames(all)%in%trait_list])

trait_list <- sort(trait_list)
betternames <- trait_list

#i can't get the usual lsmeans to work. so i'll just do actual regional means and se.
reg_mean <- ddply(plotdata,.(region,age,region_age), colwise(mean,na.rm=T))
reg_sd <- ddply(plotdata,.(region,age,region_age), colwise(sd,na.rm=T))
reg_N <- ddply(plotdata,.(region,age,region_age), colwise(length))
reg_se <- data.frame(reg_sd[1:3] ,reg_sd[4:length(reg_sd)]/sqrt(reg_N[4:length(reg_N)]))

for (i in 1:length(trait_list)){
  df <- data.frame(reg_mean[1:3],reg_mean[, colnames(reg_mean)%in%trait_list[i]],reg_se[, colnames(reg_mean)%in%trait_list[i]])
  colnames(df) <- c("region","age","age_reg","y_mean","y_se")
  plot = ggplot(df,aes(x=region, y=y_mean, group = age)) + geom_bar(stat = "identity", position = position_dodge(width =.9), aes(fill = age_reg, color = age_reg)) + scale_fill_manual(values=c("navyblue","steelblue1","maroon2","white","white","white")) + scale_color_manual(values=c("navyblue","steelblue1","maroon2","navyblue","steelblue1","maroon2")) + geom_errorbar(position = position_dodge(width =.9), aes(ymin=y_mean-y_se, ymax=y_mean+y_se),width=.2,size=.4,color="black") 
  vv = plot  +ylab(betternames[i]) + theme_classic() + theme(text = element_text(size = 8),legend.position = "none") + scale_x_discrete(limits=c("tropical","subtropical", "temperate")) 
  print(vv)
  # pdf(paste("20180518_trait_region",i,".pdf"),height=2.25,width=2.5)
  # print(vv)
  # dev.off()
}

##################
#conclusions
##################
#changing the cutoff to include 0.5% relative abundance (207 peaks) or 0.1% relative abundance (793 peaks) does not change the results for regional analyses of overall abundance at all, versus 1% cutoff (110 peaks). For richness and diversity, the results are similar for 1% and 0.5%, while differences start to become erased at 0.1%