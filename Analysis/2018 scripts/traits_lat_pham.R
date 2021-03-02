#analysis of traits as a function of latitude

#setup
####################
source("/Users/carina/Documents/dissertation/Phytolacca writing/ch3 defensive traits/ch3 defensive traits analysis/all_traits_wrangling_mixed.R")
str(all)
library(nlme)
library(corrplot)
library(Hmisc)
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
all <- droplevels(subset(all, all$lat > 20))

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
  plotdata.y <- cbind(all[, colnames(all)%in%trait_list[i]], all[3])
  boxplot(plotdata.y[,1]~plotdata.y[,2],
          data=na.omit(plotdata.y),
          main=paste(trait_list[i],"by pop_young"), 
          xlab="population" ,
          ylab = trait_list[i])
}
for (i in 1:length(trait_list)){  
  plotdata.m <- subset(all, all$age == "mature")
  plotdata.m <- cbind(all[, colnames(all)%in%trait_list[i]], all[3])
  boxplot(plotdata.m[,1]~plotdata.m[,2],
          data=na.omit(plotdata.m),
          main=paste(trait_list[i],"by pop_mature"), 
          xlab="population" ,
          ylab = trait_list[i])
}
#nearly all have pretty bad heterogeneity of variance, so maybe the default should be include the weights argument, and check that results are not changed for those traits that probably don't need it (if changed, see if they do need it)

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
  fmla <- as.formula(paste(paste(trait),"~ lat*age"))
  lme_model <- lme(fmla, random=~1|pop, method = "REML", data=all, na.action = na.exclude)
  a<- hist((resid(lme_model) - mean(resid(lme_model), na.rm=T)) / sd(resid(lme_model), na.rm=T), freq=F, main = paste(trait,sep="_"),xlab="residuals"); curve(dnorm, add = TRUE)
  b<-qqnorm(lme_model, main = paste(trait,sep="_"))
  c<-plot(lme_model, main = paste(trait,sep="_"))
  list(a,b,c)
  #, weights=varIdent(form=~1|pop_age)
}

ass_check_weight <- function(trait){
  fmla <- as.formula(paste(paste(trait),"~ lat*age"))
  lme_model <- lme(fmla, random=~1|pop, method = "REML", data=all, na.action = na.exclude, weights=varIdent(form=~1|pop))
  a<- hist((resid(lme_model) - mean(resid(lme_model), na.rm=T)) / sd(resid(lme_model), na.rm=T), freq=F, main = paste(trait,sep="_"),xlab="residuals"); curve(dnorm, add = TRUE)
  b<-qqnorm(lme_model, main = paste(trait,sep="_"))
  c<-plot(lme_model, main = paste(trait,sep="_"))
  list(a,b,c)
}

#again, can't figure out how to condense these figures, so I'll have to do it one at a time

ass_check(trait_list[1])
ass_check_weight(trait_list[1])
#yes but histogram is a little tall. I tried log and sqrt transform (before standardizing) and checked models, and normality improves a little bit but results are the same, so I think don't worry about it
ass_check(trait_list[2])
ass_check_weight(trait_list[2])
#yes; needs pop_age
ass_check(trait_list[3])
ass_check_weight(trait_list[3])
#maybe
ass_check(trait_list[4])
ass_check_weight(trait_list[4])
#yes, needs pop_age
ass_check(trait_list[5])
ass_check_weight(trait_list[5])
#no
ass_check(trait_list[6])
ass_check_weight(trait_list[6])
#yes
ass_check(trait_list[7])
ass_check_weight(trait_list[7])
#yes
ass_check(trait_list[8])
ass_check_weight(trait_list[8])
#yes

#in sum, these do not need weights argument: NMDS1
#these need it (pop): N, maybe richness, NMDS2-4
#these need it (pop_age): toughness and abund


#big tables with results of testing effects of lat and age on each trait
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
lme_slopes <- function(lme_model,trait)
{
  b_mature = summary(lme_model)$coef$fixed[1]
  m_mature = summary(lme_model)$coef$fixed[2]
  b_young = summary(lme_model)$coef$fixed[1]+summary(lme_model)$coef$fixed[3]
  m_young = summary(lme_model)$coef$fixed[2]+summary(lme_model)$coef$fixed[4]
  coef <- as.vector(c(b_mature,m_mature,b_young,m_young))
  names(coef) <- c("b_mature","m_mature","b_young","m_young")
  slopes <- as.data.frame(coef)
  colnames(slopes) <- trait
  return(slopes)
}

pred_table <- data.frame(  F=numeric(),
                           numDF=numeric(),
                           denDF=numeric(),
                           p=numeric(),
                         trait=character(),
                           stringsAsFactors=FALSE)
slopes <- data.frame(matrix(0, nrow = 4, ncol = 0))

#no weights
trait_list_sub <- trait_list[5]
for (i in 1:length(trait_list_sub)){
  fmla <- as.formula(paste(paste(trait_list_sub[i]),"~ lat*age"))
  mod <- lme(fmla, random=~1|pop, method = "REML", data=all, na.action = na.exclude)
  results <- cbind(as.data.frame(anova(mod))[2:4,],trait_list_sub[i])
  pred_table <- rbind(pred_table,results)
  slopes <- cbind(slopes,lme_slopes(mod,trait_list_sub[i]))
}
#pop weights
trait_list_sub <- trait_list[c(1,3,6:8)]
for (i in 1:length(trait_list_sub)){
  fmla <- as.formula(paste(paste(trait_list_sub[i]),"~ lat*age"))
  mod <- lme(fmla, random=~1|pop, method = "REML", data=all, na.action = na.exclude, weights=varIdent(form=~1|pop))
  results <- cbind(as.data.frame(anova(mod))[2:4,],trait_list_sub[i])
  pred_table <- rbind(pred_table,results)
  slopes <- cbind(slopes,lme_slopes(mod,trait_list_sub[i]))
}
#pop_age weights
trait_list_sub <- trait_list[c(2,4)]
for (i in 1:length(trait_list_sub)){
  fmla <- as.formula(paste(paste(trait_list_sub[i]),"~ lat*age"))
  mod <- lme(fmla, random=~1|pop, method = "REML", data=all, na.action = na.exclude, weights=varIdent(form=~1|pop_age))
  results <- cbind(as.data.frame(anova(mod))[2:4,],trait_list_sub[i])
  pred_table <- rbind(pred_table,results)
  slopes <- cbind(slopes,lme_slopes(mod,trait_list_sub[i]))
}
pred_table
slopes

#remove slopes if there was no effect of lat or lat*age
slopes$richness[2] = NA
slopes$richness[4] = NA
slopes$abund[2] = NA
slopes$abund[4] = NA
slopes$NMDS1[2] = NA
slopes$NMDS1[4] = NA
slopes$NMDS3[2] = NA
slopes$NMDS3[4] = NA

###################
#plots of those with significant relationships. Here, population means.
###################

plotdata <- data.frame(all[c("pop","age","lat")], all[, colnames(all)%in%trait_list])
  
betternames <- c("% Nitrogen\n","leaf toughness\n(g force to pierce)","chem. richness\n","chem. abundance\n(peak area/IS/mgdw)","chem. NMDS1 scores\n","chem. NMDS2 scores\n","chem. NMDS3 scores\n","chem. NMDS4 scores\n")

pop_mean <- ddply(plotdata,.(pop,age), colwise(mean,na.rm=T))
pop_sd <- ddply(plotdata,.(pop,age), colwise(sd,na.rm=T))
pop_N <- ddply(plotdata,.(pop,age), colwise(length))
pop_se <- data.frame(pop_sd[1:3] ,pop_sd[4:length(pop_sd)]/sqrt(pop_N[4:length(pop_N)]))

for (i in 1:length(trait_list)){
  b_mature <- slopes[1,colnames(slopes)%in%trait_list[i]]
  m_mature <- slopes[2,colnames(slopes)%in%trait_list[i]]
  b_young <- slopes[3,colnames(slopes)%in%trait_list[i]]
  m_young <- slopes[4,colnames(slopes)%in%trait_list[i]]
  
  df <- data.frame(pop_mean[2],pop_mean[3], pop_mean[, colnames(pop_mean)%in%trait_list[i]],pop_se[, colnames(pop_mean)%in%trait_list[i]])
  colnames(df) <- c("age","lat","y_mean","y_se")
  df$lat <- ifelse(df$age == "young",df$lat-0.3,df$lat)
  plot = ggplot(df,aes(x=lat, y=y_mean, group = age)) + geom_point(size = 2, aes(shape = age))  + scale_shape_manual(values=c(19,1)) + geom_errorbar(aes(ymin=y_mean-y_se, ymax=y_mean+y_se),size=.4,width=0) + xlab("") +ylab(betternames[i])  + geom_abline(slope=m_young, intercept= b_young,linetype="dashed") + geom_abline(slope=m_mature, intercept= b_mature,linetype="solid")
  vv = plot + theme_classic() + theme(text = element_text(size = 10),legend.position = "none")
  print(vv)
  pdf(paste("20180518_trait_lat",i,".pdf"),height=2.25,width=2.5)
  print(vv)
  dev.off()
}
