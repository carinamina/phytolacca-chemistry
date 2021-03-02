#analysis of traits as a function of region for PHRI

#setup
####################
source("/Users/carina/Documents/dissertation/Phytolacca writing/ch3 defensive traits/ch3 defensive traits analysis/all_traits_wrangling_mixed.R")
str(all)
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

trait_list <- colnames(all[8:ncol(all)])
#remove percent carbon and diversity
trait_list <- trait_list[-2]
trait_list <- trait_list[-4]
trait_list <- trait_list[-1]
trait_list

all <- data.frame(all[,1:6],all[,colnames(all)%in%trait_list])
all <- droplevels(subset(all, is.na(all$region) == FALSE))

#get each trait as a standardized z-score (subtract mean and divide by standard deviation) so we can put coefficients in a table and the slopes are comparable
# for(i in 1:length(trait_list)){
#   all[,colnames(all)==trait_list[i]] <- scale(all[,colnames(all)==trait_list[i]], center=TRUE, scale = TRUE)
# }

####################
#visualization
#####################

#check for heterogeneity of variance in each trait
par(mfrow=c(3,3))
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
#not as bad as normal

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
  lme_model <- lme(fmla, random=~1|pop, method = "REML", data=all, na.action = na.exclude, weights=varIdent(form=~1|pop))
  a<- hist((resid(lme_model) - mean(resid(lme_model), na.rm=T)) / sd(resid(lme_model), na.rm=T), freq=F, main = paste(trait,sep="_"),xlab="residuals"); curve(dnorm, add = TRUE)
  b<-qqnorm(lme_model, main = paste(trait,sep="_"))
  c<-plot(lme_model, main = paste(trait,sep="_"))
  list(a,b,c)
}

#again, can't figure out how to condense these figures, so I'll have to do it one at a time

ass_check(trait_list[1])
ass_check_weight(trait_list[1])
#yes
ass_check(trait_list[2])
ass_check_weight(trait_list[2])
#yes; needs pop_age
ass_check(trait_list[3])
ass_check_weight(trait_list[3])
#yes
ass_check(trait_list[4])
ass_check_weight(trait_list[4])
#yes, needs pop_age
ass_check(trait_list[5])
ass_check_weight(trait_list[5])
#no
ass_check(trait_list[6])
ass_check_weight(trait_list[6])
#no
ass_check(trait_list[7])
ass_check_weight(trait_list[7])
#no
ass_check(trait_list[8])
ass_check_weight(trait_list[8])
#yes

#in sum, these do not need weights argument: NMDS1-3
#these need it (pop): N, maybe richness, NMDS4
#these need it (pop_age): toughness and abund


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
trait_list_sub <- trait_list[c(3,5:7)]
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
#pop weights
trait_list_sub <- trait_list[c(1,8)]
for (i in 1:length(trait_list_sub)){
  fmla <- as.formula(paste(paste(trait_list_sub[i]),"~ region*age"))
  mod <- lme(fmla, random=~1|pop, method = "REML", data=all, na.action = na.exclude, weights=varIdent(form=~1|pop))
  results <- cbind(as.data.frame(anova(mod))[2:4,],trait_list_sub[i])
  pred_table <- rbind(pred_table,results)
  
  fmla2 <- as.formula(paste(paste(trait_list_sub[i]),"~ region_age"))
  mod2 <- lme(fmla2 , random=~1|pop, method = "REML", data=all, na.action = na.exclude, weights=varIdent(form=~1|pop))
  tuk <- glht(mod2, linfct=mcp(region_age="Tukey"))
  summary(tuk)  
  tuk.cld <- as.data.frame(cld(tuk)$mcletters$Letters)
  colnames(tuk.cld) <- trait_list_sub[i]
  letters <- cbind(letters,tuk.cld)
}
#pop_age weights
trait_list_sub <- trait_list[c(2,4)]
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
#these two traits don't have significant region*age interaction
letters$richness <- NA
letters$NMDS4 <- NA

letters2 <- data.frame(matrix(0, nrow = 3, ncol = 0))
#NMDS4
    fmla2 <- as.formula(paste(paste(trait_list[8]),"~ region*age"))
    mod2 <- lme(fmla2 , random=~1|pop, method = "REML", data=all, na.action = na.exclude, weights=varIdent(form=~1|pop))
    tuk <- glht(mod2, linfct=mcp(region="Tukey"))
    summary(tuk)  
    tuk.cld <- as.data.frame(cld(tuk)$mcletters$Letters)
    colnames(tuk.cld) <- trait_list[8]
    letters2 <- cbind(letters2,tuk.cld)
    
#richness
    fmla2 <- as.formula(paste(paste(trait_list[3]),"~ region*age"))
    mod2 <- lme(fmla2 , random=~1|pop, method = "REML", data=all, na.action = na.exclude, weights=varIdent(form=~1|pop))
    tuk <- glht(mod2, linfct=mcp(region="Tukey"))
    summary(tuk)  
    tuk.cld <- as.data.frame(cld(tuk)$mcletters$Letters)
    colnames(tuk.cld) <- trait_list[3]
    letters2 <- cbind(letters2,tuk.cld)

letters2

###################
#plots of regional means.
###################

plotdata <- data.frame(all[c("age","region","region_age")], all[, colnames(all)%in%trait_list])

betternames <- c("% Nitrogen\n","leaf toughness\n(g force to pierce)","chem. richness\n","chem. abundance\n(peak area/IS/mgdw)","chem. NMDS1 scores\n","chem. NMDS2 scores\n","chem. NMDS3 scores\n","chem. NMDS4 scores\n")

#i can't get the usual lsmeans to work. so i'll just do actual regional means and se.
reg_mean <- ddply(plotdata,.(region,age,region_age), colwise(mean,na.rm=T))
reg_sd <- ddply(plotdata,.(region,age,region_age), colwise(sd,na.rm=T))
reg_N <- ddply(plotdata,.(region,age,region_age), colwise(length))
reg_se <- data.frame(reg_sd[1:3] ,reg_sd[4:length(reg_sd)]/sqrt(reg_N[4:length(reg_N)]))

for (i in 1:length(trait_list)){
  df <- data.frame(reg_mean[1:3],reg_mean[, colnames(reg_mean)%in%trait_list[i]],reg_se[, colnames(reg_mean)%in%trait_list[i]])
  colnames(df) <- c("region","age","age_reg","y_mean","y_se")
  plot = ggplot(df,aes(x=region, y=y_mean, group = age)) + geom_bar(stat = "identity", position = position_dodge(width =.9), aes(fill = age_reg, color = age_reg)) + scale_fill_manual(values=c("navyblue","steelblue1","maroon2","white","white","white")) + scale_color_manual(values=c("navyblue","steelblue1","maroon2","navyblue","steelblue1","maroon2")) + geom_errorbar(position = position_dodge(width =.9), aes(ymin=y_mean-y_se, ymax=y_mean+y_se),width=.2,size=.4,color="black") 
  vv = plot  +ylab(betternames[i]) + theme_classic() + theme(text = element_text(size = 10),legend.position = "none") + scale_x_discrete(limits=c("tropical","subtropical", "temperate")) 
  print(vv)
  pdf(paste("20180518_trait_region",i,".pdf"),height=2.25,width=2.5)
  print(vv)
  dev.off()
}

letters <- letters[c(4,7,5,8,1,2,3,6)]
letters <- as.data.frame(t(letters))
letters <- letters[c(3,6,1,4,2,5)]
letters2
