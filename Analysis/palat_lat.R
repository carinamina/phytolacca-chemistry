#Analysis of palatability by lat and region
#5/12/21

###############
#setup
###############
# library(afex)
# library(pbkrtest)
# library(MuMIn)
# library(arm)
# library(effects)
# library(lsmeans)
# library(car)
library(nlme)
# library(bbmle)
library(multcomp)
# library(plyr)
library(ggplot2)
library(tidyverse)
library(ggh4x)

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
rm(palat, trop, tuk, trop_int)

######################
#plots 
#it's very confusing, this whole business about log(mean) versus mean(log). I want to plot the population means and SE because jittering points that represent lines or cups makes the latitude distribution look much more continuous than it really is. I tried to get population mean from the original, untransformed area, but then how do you log-transform and then standardize SE by initial and duration?! I think the best will be to calculate mean and standard error of the line-level data (?!)
######################
#first calculate population means and regional means
palat_line <- na.omit(read.csv("Processing/2_out_AllTraits.csv",header=T) %>% select(pop,line_age,line,age,lat,region,species,log.area))

palat_pop <- palat_line %>% mutate(pop_age = paste(pop,age,sep="_")) %>% group_by(pop_age,region,age,species) %>% summarise(area_mean = mean(log.area,na.rm=T), area_sd = sd(log.area,na.rm=T), area_n = length(pop_age), area_se = area_sd/sqrt(area_n), lat = mean(lat)) %>% select(-c(area_n,area_sd))
palat_pop$species = factor(palat_pop$species, levels=c("PHRI","PHAM"))

palat_region <- palat_line %>% mutate(region.age = paste(region,age,sep=".")) %>% group_by(region.age,region,age,species) %>% summarise(area_mean = mean(log.area,na.rm=T), lat.mean = mean(lat), lat.min = min(lat), lat.max = max(lat))

#here would be the time to compare means between regions


#convert the Tukey letter results to a dataframe for plotting letters
tukey <- as.data.frame(tuk.cld$mcletters$Letters)
colnames(tukey) <- "letter"
tukey <-tukey %>% add_column(region.age = rownames(tukey)) %>% right_join(palat_region,by="region.age") %>% mutate(plot_y = area_mean+0.11)
tukey$species = factor(tukey$species, levels=c("PHRI","PHAM"))
tukey[tukey$region=="temperate",c("lat.mean","lat.min","lat.max","plot_y")] = NA
tukey

rm(tuk.cld,palat_line,palat_region)

######################
#faceting
fac = ggplot(palat_pop,aes(x=lat, y=area_mean,group=region))+  
  geom_point(size=3,stroke=1,aes(shape = age, color=region)) + 
  geom_errorbar(aes(ymin=area_mean-area_se, ymax=area_mean+area_se),width=.2,size=.4) +
  scale_shape_manual(values=c(19,1),name="Leaf age") + 
  scale_colour_manual(values=c("steelblue1","navyblue","grey","maroon2")) +
  geom_text(data = tukey, aes(x = lat.mean, y = plot_y, label = letter,color = region), size = 6) +
  geom_segment(data = tukey, aes(x = lat.min, xend = lat.max, y = plot_y-.025, yend = plot_y-.025, color = region)) +
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
  geom_text(data = tukey, aes(x = lat.mean, y = plot_y, label = letter,color = region), size = 6) +
  scale_colour_manual(values=c("steelblue1","navyblue","grey","maroon2")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 20), legend.position = "none") + 
  ylab(bquote('ln(leaf area consumed ('~mm^2*'))')) + 
  xlab("Latitude (N)")
vv

rm(list=ls())

