#analysis of traits predicting palatability using mixed models

source("/Users/carina/Documents/dissertation/Phytolacca writing/ch3 defensive traits/ch3 defensive traits analysis/all_traits_wrangling_mixed.R")
original <- all
rm(key)
trait_list <- c("percent_N","tough","richness","abund","NMDS1","NMDS2","NMDS3","NMDS4")
all <- data.frame(all[c("line_age","lat","pop","age","conv")], all[,colnames(all)%in%trait_list])
all <- na.omit(all)
str(all)
young <- subset(all, all$age == "young")
mature <- subset(all, all$age == "mature")

library(nlme)
library(corrplot)
library(Hmisc)
library(MuMIn)
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

#check for heterogeneity of variance in caterpillar conv
#all pops, mature
boxplot(conv~pop,
        data=subset(all, all$age == "mature"), 
        main="mature conv by pop", 
        xlab="population", 
        ylab="conv")
#yes

#all pops young
boxplot(conv~pop,
        data=subset(all, all$age == "young"), 
        main="young conv by pop", 
        xlab="population", 
        ylab="conv")
#yes

#will use the weights argument for all analyses to be consistent, it will probably be necessary

#####################
#correlation matrices
#####################

traits.y <- as.matrix(subset(all, all$age == "young")[,colnames(all)%in%trait_list])
colnames(traits.y) <- c("% nitrogen","leaf toughness","chemical richness","chemical abundance","chemical NMDS1","chemical NMDS2","chemical NMDS3","chemical NMDS4")
cor.y <-rcorr(traits.y)
cor.matrix.y <- cor.y$r
corrplot(cor.y$r, type = "upper", tl.col = "black", tl.srt = 45, p.mat = cor.y$P, sig.level = 0.05, insig = "blank", order = "original")
# setEPS()
# postscript("20180531_corplot.y.eps")
# corrplot(cor.y$r, type = "upper", tl.col = "black", tl.srt = 45, p.mat = cor.y$P, sig.level = 0.05, insig = "blank", order = "original")
# dev.off()
rm(traits.y,cor.y,cor.matrix.y)

traits.m <- as.matrix(subset(all, all$age == "mature")[,colnames(all)%in%trait_list])
colnames(traits.m) <- c("% nitrogen","leaf toughness","chemical richness","chemical abundance","chemical NMDS1","chemical NMDS2","chemical NMDS3","chemical NMDS4")
cor.m <-rcorr(traits.m)
cor.matrix.m <- cor.m$r
corrplot(cor.m$r, type = "upper", tl.col = "black", tl.srt = 45, p.mat = cor.m$P, sig.level = 0.05, insig = "blank", order = "original")
# setEPS()
# postscript("20180531_corplot.m.eps")
# corrplot(cor.m$r, type = "upper", tl.col = "black", tl.srt = 45, p.mat = cor.m$P, sig.level = 0.05, insig = "blank", order = "original")
# dev.off()
rm(traits.m,cor.m,cor.matrix.m)

#####################
#predicting palatability
#####################

#YOUNG
model<-lme(conv~percent_N + tough + richness + abund + NMDS1 + NMDS2 + NMDS3 + NMDS4, random=~1|pop, method = "REML", data=young, weights=varIdent(form=~1|pop) )
summary(model)
dredge.results<-dredge(model)
head(dredge.results)
#write.csv(as.data.frame(dredge.results), “dredge_young.csv”)

#lots of non-convergence warnings, so let's try the other optimizer and see if that helps. but there was a pretty clear winner! yay!
ctrl <- lmeControl(opt = 'optim')
model2<-lme(conv~percent_N + tough + richness + abund + NMDS1 + NMDS2 + NMDS3 + NMDS4, random=~1|pop, method = "REML", control = ctrl, data=young, weights=varIdent(form=~1|pop) )
summary(model2)
dredge.results2<-dredge(model2)
head(dredge.results2)
#ok this time only had 10 non-convergence, that's way better than the >50 last time. and results are exactly the same.
write.csv(as.data.frame(dredge.results2), "20180531dredge_young.csv")

##############FINAL MODEL YOUNG
mod.y<-lme(conv~ NMDS2 + abund , random=~1|pop, method = "REML", data=young, weights=varIdent(form=~1|pop) )
lme_results(mod.y)
#check for multicolinearity
mod.X <- model.matrix(mod.y)
eigen.x <- eigen(t(mod.X) %*%mod.X)
eigen.x$val # eigenvalues from the design matrix; none should be very close to zero; these are ok
sqrt(eigen.x$val[1]/eigen.x$val) # condition numbers; shouldn't be over 30; ok









#MATURE
model3<-lme(conv~percent_N + tough + richness + abund + NMDS1 + NMDS2 + NMDS3 + NMDS4, random=~1|pop, method = "REML", data=mature, weights=varIdent(form=~1|pop) )
summary(model3)
dredge.results3<-dredge(model3)
head(dredge.results3)
write.csv(as.data.frame(dredge.results3), "20180531dredge_mature.csv")
#only had 6 warnings

##########FINAL MODEL MATURE
mod.m<-lme(conv~  NMDS2 + NMDS4 , random=~1|pop, method = "REML", data=mature, weights=varIdent(form=~1|pop) )
lme_results(mod.m)
#check for multicolinearity
mod.X <- model.matrix(mod.m)
eigen.x <- eigen(t(mod.X) %*%mod.X)
eigen.x$val # eigenvalues from the design matrix; none should be very close to zero; these are ok
sqrt(eigen.x$val[1]/eigen.x$val) # condition numbers; shouldn't be over 30; ok






###################
#plots of correlations between traits and palatability
###################

#parameters for all the plots

param_list <- list( ylab(""),
                      scale_y_continuous(limits = c(0,1.8)) ,
                      theme_classic(),
                      theme(text = element_text(size = 10))
)

#############   YOUNG PLOTS
summary(mod.y)
b.y <- summary(mod.y)$coef$fixed[1]
m_NMDS2.y <- summary(mod.y)$coef$fixed[2]
m_abund.y <- summary(mod.y)$coef$fixed[3]
plotdata.y <- subset(all, all$age == "young")[c("pop","conv","NMDS2","abund")]

#indiv datapoints, not pop means
df <- data.frame(plotdata.y$conv,plotdata.y$NMDS2)
colnames(df) <- c("conv_mean","x_mean")
vv = ggplot(df,aes(x=x_mean, y=conv_mean)) +  xlab("chemical NMDS2 scores\n") + geom_abline(slope=m_NMDS2.y, intercept= b.y) + geom_point(size = 2, shape=1) + param_list
print(vv)
pdf(paste("20180531_traits1.pdf"),height=2.1,width=2.1)
print(vv)
dev.off()

df <- data.frame(plotdata.y$conv,plotdata.y$abund)
colnames(df) <- c("conv_mean","x_mean")
vv = ggplot(df,aes(x=x_mean, y=conv_mean)) +  xlab("chemical abundance\n(peak area/IS/mgdw)") + geom_abline(slope=m_abund.y, intercept= b.y) + geom_point(size = 2, shape=1) + param_list
print(vv)
pdf(paste("20180531_traits2.pdf"),height=2.1,width=2.1)
print(vv)
dev.off()

#############   MATURE PLOTS
summary(mod.m)
b.m <- summary(mod.m)$coef$fixed[1]
m_NMDS2.m <- summary(mod.m)$coef$fixed[2]
m_NMDS4.m <- summary(mod.m)$coef$fixed[3]
plotdata.m <- subset(all, all$age == "mature")[c("pop","conv","NMDS2","NMDS4")]

#indiv datapoints, not pop means
df <- data.frame(plotdata.m$conv,plotdata.m$NMDS2)
colnames(df) <- c("conv_mean","x_mean")
vv = ggplot(df,aes(x=x_mean, y=conv_mean)) +  xlab("chemical NMDS2 scores\n") + geom_abline(slope=m_NMDS2.m, intercept= b.m) + geom_point(size = 2, shape=19)+ param_list
print(vv)
pdf(paste("20180531_traits3.pdf"),height=2.1,width=2.1)
print(vv)
dev.off()

df <- data.frame(plotdata.m$conv,plotdata.m$NMDS4)
colnames(df) <- c("conv_mean","x_mean")
vv = ggplot(df,aes(x=x_mean, y=conv_mean)) +  xlab("chemical NMDS4 scores\n") + geom_abline(slope=m_NMDS4.m, intercept= b.m) + geom_point(size = 2, shape=19)+ param_list
print(vv)
pdf(paste("20180531_traits4.pdf"),height=2.1,width=2.1)
print(vv)
dev.off()









# param_list_se <- list(geom_point(size = 2, shape=19) , 
#                       geom_errorbar(aes(ymin=conv_mean-conv_se, ymax=conv_mean+conv_se),size=.4,width=0),
#                       geom_errorbarh(aes(xmin=x_mean-x_se, xmax=x_mean+x_se),size=.4,height=0),
#                       ylab(""),
#                       scale_y_continuous(limits = c(0,1.8)) ,
#                       theme_classic(),
#                       theme(text = element_text(size = 8))
# )


# #young, population means
# mean.y <- ddply(plotdata.y,.(pop), colwise(mean,na.rm=T))
# sd.y <- ddply(plotdata.y,.(pop), colwise(sd,na.rm=T))
# N.y <- ddply(plotdata.y,.(pop), colwise(length))
# se.y <- data.frame(sd.y[1] ,sd.y[2:length(sd.y)]/sqrt(N.y[2:length(N.y)]))
# 
# df <- data.frame(mean.y$conv,se.y$conv,mean.y$NMDS2,se.y$NMDS2)
# colnames(df) <- c("conv_mean","conv_se","x_mean","x_se")
# vv = ggplot(df,aes(x=x_mean, y=conv_mean)) +  xlab("chemical NMDS2 scores\n") + geom_abline(slope=m_NMDS2.y, intercept= b.y) + param_list_se
# print(vv)
# 
# df <- data.frame(mean.y$conv,se.y$conv,mean.y$abund,se.y$abund)
# colnames(df) <- c("conv_mean","conv_se","x_mean","x_se")
# vv = ggplot(df,aes(x=x_mean, y=conv_mean)) +  xlab("chemical abundance\n(peak area/IS/mgdw)") + geom_abline(slope=m_abund.y, intercept= b.y) + param_list_se
# print(vv)



# #mature, population means
# mean.m <- ddply(plotdata.m,.(pop), colwise(mean,na.rm=T))
# sd.m <- ddply(plotdata.m,.(pop), colwise(sd,na.rm=T))
# N.m <- ddply(plotdata.m,.(pop), colwise(length))
# se.m <- data.frame(sd.m[1] ,sd.m[2:length(sd.m)]/sqrt(N.m[2:length(N.m)]))
# 
# df <- data.frame(mean.m$conv,se.m$conv,mean.m$NMDS2,se.m$NMDS2)
# colnames(df) <- c("conv_mean","conv_se","x_mean","x_se")
# vv = ggplot(df,aes(x=x_mean, y=conv_mean)) + xlab("chemical NMDS2 scores\n") + geom_abline(slope=m_NMDS2.m, intercept= b.m) + param_list_se
# print(vv)
# 
# df <- data.frame(mean.m$conv,se.m$conv,mean.m$NMDS4,se.m$NMDS4)
# colnames(df) <- c("conv_mean","conv_se","x_mean","x_se")
# vv = ggplot(df,aes(x=x_mean, y=conv_mean)) +  xlab("chemical NMDS4 scores\n") + geom_abline(slope=m_NMDS4.m, intercept= b.m) + param_list_se
# print(vv)

