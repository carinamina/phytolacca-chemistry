#Analysis of palatability by lat and region
#5/12/21

###############
#setup
###############
library(nlme)
library(multcomp)
library(ggplot2)
library(tidyverse)
library(ggh4x)
#library(lsmeans)

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
  slopes = data.frame(species = c("PHRI","PHRI","PHAM","PHAM"),age = c("mature","young"), m = c(0,0,NA,NA), b = c(1,1,NA,NA))
  slopes$species = factor(slopes$species, levels=c("PHRI","PHAM"))
  #mature slope
  slopes[3,3] = summary(lme_model)$coef$fixed[2]
  #mature intercept
  slopes[3,4] = summary(lme_model)$coef$fixed[1]
  #young slope
  slopes[4,3] = summary(lme_model)$coef$fixed[2]+summary(lme_model)$coef$fixed[4]
  #young intercept
  slopes[4,4] = summary(lme_model)$coef$fixed[1]+summary(lme_model)$coef$fixed[3]
  return(slopes)
}

###############
#import data
###############
palat <- na.omit(read.csv("Processing/1e_out_Palat_Cup.csv", header = TRUE))
#omitting NA removes 10 lines
str(palat)

#boxplots
#####################
boxplot(log.area~pop*age,
        data= palat, 
        main="palatability by pop", 
        xlab="population", 
        ylab="leaf area consumed")
#yes heterogeneity of variance by population and age

boxplot(log.area~region*age,
        data= palat, 
        main="palatability by region", 
        xlab="region", 
        ylab="leaf area consumed")
#yes heterogeneity of variance, we can assume the weights argument is needed

#gradient analysis of PHAM
#####################
pham <- lme(log.area ~ lat*age , random=~1|pop/line, method = "REML", data=palat[palat$lat>20,], weights=varIdent(form=~1|pop))
lme_results(pham)
#assumption fitting is ok, a bit of a long tail in the normality of residuals but this assumption is fairly robust. using pop_age in weights argument doesn't change it much
(s <- lme_slopes(pham))
#very significant lat*age and slopes are quite different

#Note May 19: I belatedly realized I could probably transform area in the formula (log(area/initial+1)/duration). Results are virtually identical to above. I thought this would be helpful because I could get lsmeans and not worry about the log-mean issue, but turns out you can't get lsmeans for a blocking factor, only fixed categorical effects, so I decided to not change anything.

rm(pham)
#regional analysis with tropical plants
#####################

trop <- lme(log.area ~ region*age , random=~1|pop/line, method = "REML", data=palat[palat$region!="temperate",], weights=varIdent(form=~1|pop))
lme_results(trop)
#similar assumption fit as above; pop_age again doesn't change anything

#tukey test for interaction
palat$reg_age <- interaction(palat$region, palat$age)
trop_int <- lme(log.area ~ reg_age, random = ~1|pop/line, method = "REML", data=droplevels(palat[palat$region!="temperate",]), weights=varIdent(form=~1|pop)) 
tuk <- glht(trop_int, linfct=mcp(reg_age="Tukey"))
summary(tuk)          # standard display
(tuk.cld <- cld(tuk))   # letters
#this plot isn't showing all the letters
# opar <- par(mai=c(2,2,2,2))
# plot(tuk.cld,las=2)
# par(opar)
rm(trop, tuk, trop_int)

######################
#plots 
#it's very confusing, this whole business about log(mean) versus mean(log). 
#Katka Bodova, who is much better than math than me, gave me this advice
#"This is something I did multiple times. If your natural variable is log(area/initial)/duration then you compute std as you did. But you are saying that your natural variable is log(mean_area/initial)/duration so to get a std on this variable you need to use one of the methods that estimates variation of this complex variable. The best methods to use are leave-one-out (or leave-p-out) https://en.wikipedia.org/wiki/Cross-validation_(statistics) or jackknife https://en.wikipedia.org/wiki/Jackknife_resampling "
#So I should be calculating means and standard error from log-transformed data for these plots, therefore I should use log.area from the cup-level data and calculate mean/SE as needed
######################
#first calculate population means and regional means
#palat_line <- na.omit(read.csv("Processing/2_out_AllTraits.csv",header=T) %>% select(pop,line_age,line,age,lat,region,species,log.area))

palat$species = ifelse(palat$region == "tropical","PHRI","PHAM")

#comparing means for Results
#north temperate young vs. tropical&subtropical young
palat %>% filter(region=="north temperate", age == "young") %>% summarise(mean = mean(log.area))/
  ((palat %>% filter(region=="subtropical", age == "young") %>% summarise(mean = mean(log.area)) +
      palat %>% filter(region=="tropical", age == "young") %>% summarise(mean = mean(log.area)))/2)
#tropical and north temperate mature vs. subtropical mature
((palat %>% filter(region=="north temperate", age == "mature") %>% summarise(mean = mean(log.area)) +
      palat %>% filter(region=="tropical", age == "mature") %>% summarise(mean = mean(log.area)))/2)/
  palat %>% filter(region=="subtropical", age == "mature") %>% summarise(mean = mean(log.area))
#tropical mature vs young
palat %>% filter(region=="tropical", age == "mature") %>% summarise(mean = mean(log.area))/
  palat %>% filter(region=="tropical", age == "young") %>% summarise(mean = mean(log.area))
#subtropical mature vs young
palat %>% filter(region=="subtropical", age == "mature") %>% summarise(mean = mean(log.area))/
  palat %>% filter(region=="subtropical", age == "young") %>% summarise(mean = mean(log.area))
#temperate mature vs young
palat %>% filter(region=="north temperate", age == "mature") %>% summarise(mean = mean(log.area))/
  palat %>% filter(region=="north temperate", age == "young") %>% summarise(mean = mean(log.area))

palat_pop <- palat %>% mutate(pop_age = paste(pop,age,sep="_")) %>% group_by(pop_age,region,age,species) %>% summarise(area_mean = mean(log.area,na.rm=T), area_sd = sd(log.area,na.rm=T), area_n = length(pop_age), area_se = area_sd/sqrt(area_n), lat = mean(lat)) %>% select(-c(area_n,area_sd))
palat_pop$species = factor(palat_pop$species, levels=c("PHRI","PHAM"))

palat_region <- palat %>% mutate(region.age = paste(region,age,sep=".")) %>% group_by(region.age,region,age,species) %>% summarise(area_mean = mean(log.area,na.rm=T), lat.min = min(lat), lat.max = max(lat), lat.middle = mean(c(lat.min,lat.max)))


#convert the Tukey letter results to a dataframe for plotting letters
tukey <- as.data.frame(tuk.cld$mcletters$Letters)
colnames(tukey) <- "letter"
tukey <-tukey %>% add_column(region.age = rownames(tukey)) %>% right_join(palat_region,by="region.age") %>% mutate(plot_y = area_mean+0.11)
tukey$species = factor(tukey$species, levels=c("PHRI","PHAM"))
tukey[tukey$region=="temperate",c("lat.middle","lat.min","lat.max","plot_y")] = NA
tukey

rm(tuk.cld,palat_region)

######################
#faceting
fac = ggplot(palat_pop,aes(x=lat, y=area_mean,group=region))+  
  geom_point(size=3,stroke=1,aes(shape = age, color=region)) + 
  geom_errorbar(aes(ymin=area_mean-area_se, ymax=area_mean+area_se),width=.2,size=.4) +
  scale_shape_manual(values=c(19,1),name="Leaf age") + 
  scale_colour_manual(values=c("steelblue1","navyblue","grey","maroon2")) +
  geom_text(data = tukey, aes(x = lat.middle, y = plot_y, label = letter), size = 6) +
  geom_segment(data = tukey, aes(x = lat.min, xend = lat.max, y = plot_y-.025, yend = plot_y-.025)) +
  facet_grid(.~species,scales="free_x",switch="x") +
  force_panelsizes(cols = c(0.3, 1), rows = c(1.3), respect = TRUE) +
  facetted_pos_scales(x = list(scale_x_continuous(breaks = c(9, 10), limits = c(8.4,10.7))),NULL) +
  geom_abline(data=s, aes(slope=m,intercept=b,linetype = age)) +
  scale_linetype_manual(values=c("solid", "dashed"),name="Leaf age") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.spacing = unit(0.2, "lines"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 20), strip.background = element_blank(),strip.placement="outside") + 
  ylab(bquote('ln(Leaf area consumed ('~mm^2*'))')) + 
  xlab(expression("Latitude (\u00B0N)"))
fac

setEPS()
postscript("FiguresTables/PalatLat.eps", height = 5, width = 5)
fac + theme(legend.position = "none")
dev.off()

setEPS()
postscript("FiguresTables/PalatLatLegend.eps", height = 5, width = 5)
fac
dev.off()

#simple version for playing with settings
######################
vv = ggplot(palat_pop,aes(x=lat, y=area_mean,group=region))+  
  geom_point(size=4,stroke=2,aes(shape = age, color=region)) + 
  geom_errorbar(aes(ymin=area_mean-area_se, ymax=area_mean+area_se),width=.2,size=.4) +
  scale_shape_manual(values=c(19,1)) + 
  geom_text(data = tukey, aes(x = lat.middle, y = plot_y, label = letter,color = region), size = 6) +
  scale_colour_manual(values=c("steelblue1","navyblue","grey","maroon2")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 20), legend.position = "none") + 
  ylab(bquote('ln(leaf area consumed ('~mm^2*'))')) + 
  xlab("Latitude (N)")
vv

rm(list=ls())

