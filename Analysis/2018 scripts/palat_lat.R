#Analysis of palatability by lat and region
#5/12/18

###############
#setup
###############
setwd("/Users/carina/Documents/R_working_directory")

library(afex)
library(pbkrtest)
library(MuMIn)
library(arm)
library(effects)
library(lsmeans)
library(car)
library(nlme)
library(bbmle)
library(multcomp)
library(plyr)
library(ggplot2)
library(googlesheets)

lme_results <- function(lme_model)
{
  a<- hist((resid(lme_model) - mean(resid(lme_model), na.rm=T)) / sd(resid(lme_model), na.rm=T), freq=F); curve(dnorm, add = TRUE)
  b<-qqnorm(lme_model)
  c<-plot(lme_model)
  d<-summary(lme_model)
  #e<-intervals(lme_model)[1]
  f<-anova(lme_model, type = "marginal", test = "F")
  g<-nobs(lme_model)
  
  list(a,b,c,d,f,g)
}

lme_slopes <- function(lme_model)
{
  b_mature = summary(lme_model)$coef$fixed[1]
  m_mature = summary(lme_model)$coef$fixed[2]
  b_young = summary(lme_model)$coef$fixed[1]+summary(lme_model)$coef$fixed[3]
  m_young = summary(lme_model)$coef$fixed[2]+summary(lme_model)$coef$fixed[4]
  coef <- as.vector(c(b_mature,m_mature,b_young,m_young))
  names(coef) <- c("b_mature","m_mature","b_young","m_young")
  slopes <- as.data.frame(t(coef))
  return(slopes)
}

palat <- read.csv("palatability_PHRI_20180327.csv", header = TRUE)
sum(is.na(palat))
palat$X <- NULL
str(palat)
# cupID = cup identifier (experimental unit)
# line = maternal line of diet eaten in that cup
# lat = latitude of origin of plant diet
# pop = population of origin of plant diet
# age = leaf age
# cons_init_day = ln(leaf area consumed+1) per initial caterpillar per day
# mass_init_day = ln(final larval biomass+1) per initial caterpillar per day

regions <- c("tropical", "tropical","tropical", "subtropical", "subtropical", "subtropical",NA,NA,NA,NA,NA,NA,NA, "temperate", "temperate", "temperate")
regions <- cbind(regions, sort(unique(palat$lat)))
colnames(regions) <- c("region", "lat")
palat <- merge(regions, palat, by = "lat", all = TRUE)
palat$conv <- ifelse(palat$cons_init_day >0, palat$mass_init_day/palat$cons_init_day,0)
#conv = converstion rate of mass per leaf area consumed
hist(palat$conv)
rm(regions)
palat$lat <- as.numeric(paste(palat$lat))
#boxplots
#####################
boxplot(conv~pop*age,
        data= palat, 
        main="palatability by pop", 
        xlab="population", 
        ylab="larval biomass")
#yes heterogeneity of variance

boxplot(conv~region*age,
        data= palat, 
        main="palatability by region", 
        xlab="region", 
        ylab="larval biomass")
#yes heterogeneity of variance, we can assume the weights argument is needed

#gradient analysis of PHAM
#####################
palat$pop_age <- as.factor(paste(palat$pop, palat$age, sep = "_"))

pham <- lme(conv ~ lat*age , random=~1|pop/line, method = "REML", data=subset(palat, palat$lat>20), weights=varIdent(form=~1|pop), na.action = na.exclude)
lme_results(pham)
#wow there's still fan-shaped residuals after the usual fix. maybe try estimating variance by each age-pop combo

ctrl <- lmeControl(opt = 'optim')
pham2 <- lme(conv ~ lat*age , random=~1|pop/line, method = "REML", data=subset(palat, palat$lat>20), weights=varIdent(form=~1|pop_age), control = ctrl, na.action = na.exclude)
lme_results(pham2)
(s <- lme_slopes(pham2))
#this one needs some kind of override to converge (using the "old optimizer")
#https://stats.stackexchange.com/questions/40647/lme-error-iteration-limit-reached
#assumption-fitting is really great. Effects of lat and age but not interaction

#regional analysis with tropical plants
#####################

# trop <- lme(conv ~ region*age , random=~1|pop/line, method = "REML", data=palat, weights=varIdent(form=~1|pop), na.action = na.exclude)
# lme_results(trop)
#same thing with fan-shaped residuals after estimating variance for each population separately. Need to do pop_age

trop2 <- lme(conv ~ region*age , random=~1|pop/line, method = "REML", data=palat, weights=varIdent(form=~1|pop_age), control = ctrl, na.action = na.exclude)
lme_results(trop2)

#tukey test for interaction
palat$reg_age <- interaction(palat$region, palat$age)
mixed2 <- lme(conv ~ reg_age, random = ~1|pop/line, method = "REML", data=palat, weights=varIdent(form=~1|pop_age), control = ctrl, na.action = na.exclude) 
tuk <- glht(mixed2, linfct=mcp(reg_age="Tukey"))
summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
opar <- par(mai=c(2,1,1.5,1))
plot(tuk.cld,las=2)
par(opar)

#as with biomass as palatability, tropical young is sig. lower than southern young

region_summary <- ddply(palat, c("region","age"), summarise,
                       mean = mean(conv)
)
region_summary

#subtrop young/tropical young
region_summary[2,3]/region_summary[6,3]
#northern young/tropical young
region_summary[4,3]/region_summary[6,3]
#trop&north mature/southern mature
(region_summary[5,3]+region_summary[3,3])/2/region_summary[1,3]


######################
#plots 
#population means with SE (here I can't do LS means bc not all pops were in each model)
######################

summary_stats1 <- ddply(palat, c("region","pop","lat","line","age"), summarise,
                        N = length(line),
                        total = mean(conv)
)

summary_stats <- ddply(summary_stats1, c("region","pop","lat","age"), summarise,
                       N = length(lat),
                       mean = mean(total),
                       sd = sd(total),
                       se = sd / sqrt(N)
)
summary_stats
summary_stats$region <- addNA(summary_stats$region)
#no lines for tropical part
plot = ggplot(summary_stats,aes(x=lat, y=mean,group=region))+  geom_point(size=4,stroke=2,aes(shape = age, color=region)) + scale_shape_manual(values=c(19,1)) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,size=.4)+ scale_colour_manual(values=c("navyblue","steelblue1","maroon2","grey")) 
vv=plot + xlab("") + ylab("")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 20), legend.position = "none") #+ylab(bquote('ln(biomass (mg))/ln(leaf consumed ('~mm^2*'))')) 
vv
# setEPS()
# postscript("20180531_palat_lat_phri.eps", width = 10, height = 6)
# vv
# dev.off()

#lines for PHAM
plot = ggplot(summary_stats,aes(x=lat, y=mean,group=region))+  geom_point(size=4,stroke=2,aes(shape = age, color=region)) + scale_shape_manual(values=c(19,1)) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,size=.4)+ scale_colour_manual(values=c("navyblue","steelblue1","maroon2","grey")) + geom_abline(slope=s$m_young, intercept= s$b_young, linetype = "dashed") + geom_abline(slope=s$m_mature, intercept= s$b_mature)
vv=plot + xlab("") + ylab("")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 20), legend.position = "none") #+ylab(bquote('ln(biomass (mg))/ln(leaf consumed ('~mm^2*'))')) 
vv
# setEPS()
# postscript("20180531_palat_lat_pham.eps", width = 10, height = 6)
# vv
# dev.off()



# #plot of lsmeans as bars like the tropical paper
# lsmeans <- as.data.frame(summary(lsmeans(trop2, c("region","age"))))
# lsmeans
# lsmeans$age_reg <- as.factor(paste(lsmeans$age, lsmeans$region, sep = "_"))
# 
# plot = ggplot(lsmeans,aes(x=region, y=lsmean, group = age)) + geom_bar(stat = "identity", position = position_dodge(width =.9), aes(fill = age_reg, color = age_reg)) + scale_fill_manual(values=c("navyblue","steelblue1","maroon2","white","white","white")) + scale_color_manual(values=c("navyblue","steelblue1","maroon2","navyblue","steelblue1","maroon2")) + geom_errorbar(position = position_dodge(width =.9), aes(ymin=lsmean-SE, ymax=lsmean+SE),width=.2,size=.4,color="black") 
# vv=plot + xlab("Region of origin")+ylab("ln(Total larval biomass (mg))") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 20)) + scale_x_discrete(limits=c("tropical","subtropical", "temperate")) 
# vv + theme(legend.position = "none")
# # setEPS()
# # postscript("20180417_palat.eps")
# # vv + theme(legend.position = "none")
# # dev.off()

rm(mixed2, opar, palat, trop2, tuk, tuk.cld, region_summary)
rm(plot, summary_stats, summary_stats1, vv)

