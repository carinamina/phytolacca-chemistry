#analysis and figures for analyzing latitudinal distributions of leaf traits
#started 18 May 2021, modeled after palat_lat and 2018 scripts traits_lat_pham and traits_lat_phri


#libraries
###############
library(nlme)
library(multcomp)
library(ggplot2)
library(tidyverse)
library(ggh4x)

##############
# Notes and structure of analysis
##############
#three datasets for the six traits of interest: toughness, %N, C:N, and chemical abund, richness, and diversity. Chemistry was measured on one individual per line, the others on multiple. For each trait, for a latitude and regional model, we need to make sure the nesting is correct, check assumptions to decide if weights of pop or pop_age are needed, extract model results to a big table, and then extract slopes and/or Tukey comparisons for plotting if significant. Then calculate population and regional means for plotting.
#I'll try and write functions for everything to minimize pasting big blocks.

###############
#import data, check out boxplots, calculate population and regional means
###############
#need to start over with this and actually combine the datasets in long format (with a "trait" and "value" column) before doing anything else. It's going to make everything more efficient (getting means, making summary tables, plotting, running models) to just have one big dataset. But it's a big investment!









all.traits <- read.csv("Processing/2_out_AllTraits.csv",header=T)[1:16] %>% select(-c(surv,log.area,percent_C))
trait_list <- names(all.traits)[8:ncol(all.traits)]
all.traits$pop_age <- as.factor(paste(all.traits$pop, all.traits$age, sep = "_"))

#get min, middle, and max latitude by region for plotting purposes
lat_regions <- all.traits %>% group_by(region) %>% 
  summarise(lat.min = min(lat), 
            lat.max = max(lat),
            lat.middle = mean(c(lat.min,lat.max)))

#check for heterogeneity of variance in each trait
# par(mfrow=c(2,3))
# #young
# for (i in 1:length(trait_list)){
#   plotdata.y <- subset(all.traits, all.traits$age == "young")
#   plotdata.y <- cbind(all.traits[, colnames(all.traits)%in%trait_list[i]], all.traits[1])
#   boxplot(plotdata.y[,1]~plotdata.y[,2],
#           data=na.omit(plotdata.y),
#           main=paste(trait_list[i],"by pop_young"), 
#           xlab="population" ,
#           ylab = trait_list[i])
#   rm(plotdata.y,i)
# }
# for (i in 1:length(trait_list)){  
#   plotdata.m <- subset(all.traits, all.traits$age == "mature")
#   plotdata.m <- cbind(all.traits[, colnames(all.traits)%in%trait_list[i]], all.traits[1])
#   boxplot(plotdata.m[,1]~plotdata.m[,2],
#           data=na.omit(plotdata.m),
#           main=paste(trait_list[i],"by pop_mature"), 
#           xlab="population" ,
#           ylab = trait_list[i])
#   rm(plotdata.m,i)
# }
#percent_N, C_N, diversity look like there will be problems with heterogeneity of variance. tough, log.abund, richness might be ok.
rm(all.traits)
#par(mfrow=c(1,1))

chem <- read.csv("Processing/1b_out_ChemSummaries_Indiv.csv", header=T)
length(unique(chem$line_age)) == nrow(chem)
chem$reg_age <- interaction(chem$region, chem$age)
chem$species = ifelse(chem$region == "tropical","PHRI","PHAM")

chem_region <- chem %>% mutate(region.age = paste(region,age,sep=".")) %>% 
  group_by(region.age) %>% 
  summarise(log.abund = mean(log.abund), 
            richness = mean(richness),
            diversity = mean(diversity))

chem_pop <- chem %>% mutate(pop_age = paste(pop,age,sep="_")) %>% 
  group_by(pop_age,region,age,species) %>% 
  summarise(lat = mean(lat),
            log.abund_mean = mean(log.abund), 
            log.abund_sd = sd(log.abund), 
            log.abund_n = length(pop_age), 
            log.abund_se = log.abund_sd/sqrt(log.abund_n), 
            richness_mean = mean(richness), 
            richness_sd = sd(richness), 
            richness_n = length(pop_age), 
            richness_se = richness_sd/sqrt(richness_n),
            diversity_mean = mean(diversity), 
            diversity_sd = sd(diversity), 
            diversity_n = length(pop_age), 
            diversity_se = diversity_sd/sqrt(diversity_n)) %>% 
  select(-c(log.abund_n,log.abund_sd,richness_n,richness_sd,diversity_n,diversity_sd))


tough <- read.csv("Processing/1c_out_Toughness_Indiv.csv", header=T)
length(unique(tough$line_age)) == nrow(tough) #multiple indiv per line
length(unique(tough$pos_age)) == nrow(tough)
tough$reg_age <- interaction(tough$region, tough$age)
tough$species = ifelse(tough$region == "tropical","PHRI","PHAM")

tough_region <- tough %>% mutate(region.age = paste(region,age,sep=".")) %>% 
  group_by(region.age) %>% 
  summarise(tough = mean(tough))

tough_pop <- tough %>% mutate(pop_age = paste(pop,age,sep="_")) %>% 
  group_by(pop_age) %>% 
  summarise(tough_mean = mean(tough), 
            tough_sd = sd(tough), 
            tough_n = length(pop_age), 
            tough_se = tough_sd/sqrt(tough_n)) %>% 
  select(-c(tough_n,tough_sd))

cn <- read.csv("Processing/1d_out_CarbonNitrogen_Indiv.csv", header=T)
length(unique(cn$line_age)) == nrow(cn) #multiple indiv per line
length(unique(cn$pos_age)) == nrow(cn)
cn$reg_age <- interaction(cn$region, cn$age)
cn$species = ifelse(cn$region == "tropical","PHRI","PHAM")

cn_region <- cn %>% mutate(region.age = paste(region,age,sep=".")) %>% 
  group_by(region.age) %>% 
  summarise(percent_N = mean(percent_N), C_N = mean(C_N))

cn_pop <- cn %>% mutate(pop_age = paste(pop,age,sep="_")) %>% 
  group_by(pop_age) %>% 
  summarise(percent_N_mean = mean(percent_N), 
            percent_N_sd = sd(percent_N), 
            percent_N_n = length(pop_age), 
            percent_N_se = percent_N_sd/sqrt(percent_N_n), 
            C_N_mean = mean(C_N), 
            C_N_sd = sd(C_N), 
            C_N_n = length(pop_age), 
            C_N_se = C_N_sd/sqrt(C_N_n)) %>% 
  select(-c(percent_N_n,percent_N_sd,C_N_n,C_N_sd))

#combine data
reg.means <- as.data.frame(left_join(left_join(chem_region,cn_region),tough_region))
reg.means <- reshape(reg.means, varying = colnames(reg.means)[2:ncol(reg.means)], v.names = "mean", timevar = "trait", times = colnames(reg.means)[2:ncol(reg.means)], direction = "long") %>% select(-id) %>% remove_rownames()


pop.means <- as.data.frame(left_join(chem_pop,left_join(tough_pop,cn_pop)))
#pop.means <- reshape(pop.means, varying = colnames(pop.means)[2:ncol(reg.means)], v.names = "mean", timevar = "trait", times = colnames(reg.means)[2:ncol(reg.means)], direction = "long") %>% select(-id) %>% remove_rownames()

pop.means$species = factor(pop.means$species, levels=c("PHRI","PHAM"))

rm(chem_region,cn_region,tough_region)

################
# functions and results table setup
###############
#shows results in the console, check assumptions
lme_output <- function(lme_model)
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

#creates a summary table with results (F-statistics, p-values) for all planned models
all.stats <- data.frame(trait = sort(rep(trait_list,6)), 
                        species=rep(sort(rep(c("PHAM","PHRI"),3)),6),
                        predictor="",stringsAsFactors=FALSE) %>% 
  add_column(numDF=NA,denDF=NA, "F-value"=NA, "p-value"=NA)

#adds ANOVA results of a model to the summary table which has already been created. Need to write 
#all.stats <- ano_table(a,b,c)
extract_stats <- function(response_string, species_string, lme_model,placeholder)
{
  ano <- as.data.frame(anova(lme_model))[2:4,] %>% rownames_to_column(var = "predictor") %>% add_column(trait = response_string, species = species_string, .before = "predictor") 
  all.stats[all.stats$trait == response_string & all.stats$species == species_string,] <- ano
  return(all.stats)
  rm(ano)
}

#creates summary table with all PHAM slopes
all.slopes <- data.frame(trait = sort(rep(trait_list,4))) %>% add_column(species = "", age = "", m = NA, b = NA)
all.slopes$species = factor(all.slopes$species, levels=c("PHRI","PHAM"))

#creates summary table with all the PHRI tukey comparisons. it also includes the regional means for each trait, which are helpful for calculating effect sizes and necessary for figuring out where to place letters on plots
all.tukey <- data.frame(trait = sort(rep(trait_list,6)), 
                        region = rep(rep(c("north temperate","subtropical","tropical"),2),6), 
                        age = rep(sort(rep(c("mature","young"),3) ),6)
) %>%
  mutate(species = ifelse(region == "tropical","PHRI","PHAM")) %>%
  left_join(lat_regions, by="region") %>%
  mutate(region.age = paste(region,age,sep=".")) %>%
  left_join(reg.means) %>%
  add_column(letter = "", letter_y = NA) 
all.tukey$species = factor(all.tukey$species, levels=c("PHRI","PHAM"))

#adds slopes to summary table which has already been created. For now I will keep this simple and assume a significant interaction; if not, it can be changed in all.slopes
#For PHAM, need to call as all.slopes <- extract_param(a,b,c,d) where d = 0 or NULL, whatever
#For PHRI, need to call as all.tukey <- extract_param(a,b,c,d)
extract_param <- function(response_string,species,lme_model,letter_y_constant)
{
  if(species == "PHAM"){
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
  
  slopes = add_column(slopes,trait = response_string,.before = "species")
  all.slopes[all.slopes$trait == response_string,] <- slopes
  return(all.slopes)
  rm(slopes)
  }
  else if (species == "PHRI"){
    interaction_model <- update(lme_model, fixed = paste(response_string,"~ reg_age"))
    tuk <- glht(interaction_model, linfct=mcp(reg_age="Tukey"))
    (tuk.cld <- cld(tuk))   # letters
    #convert the Tukey letter results to a dataframe for plotting letters. this seems incredibly clunky but I like being sure that I'm not just inserting letters in the larger summary data without merging somehow; they could be in the wrong order
    tukey <- as.data.frame(tuk.cld$mcletters$Letters) %>%
      `colnames<-`("letter") %>%
      rownames_to_column(var = "region.age")
    tukey <- left_join(all.tukey[all.tukey$trait == response_string,1:9] ,tukey,by="region.age")  %>%
      mutate(letter_y = mean + letter_y_constant)
    
    all.tukey[all.tukey$trait == response_string,] = tukey
    return(all.tukey)
    rm(tukey)
  }
  else{
    return("error. species is misspelled?")
  }
}


#################
# CHEMISTRY: ABUNDANCE, DIVERSITY, AND RICHNESS
#################
#MODELS.
#latitudinal gradient with PHAM
abund.pham <- lme(log.abund ~ lat*age , random=~1|pop, method = "REML", data=chem[chem$lat>20,], weights=varIdent(form=~1|pop))
#weights argument necessary to deal with fan-shaped residuals
lme_output(abund.pham) #lat and age are significant but not lat*age
(all.stats <- extract_stats("log.abund", "PHAM", abund.pham))
(all.slopes <- extract_param("log.abund", "PHAM", abund.pham))
#lat and age are significant but not their interaction, so use the mature slope for both leaf ages, with different intercepts
# all.slopes[all.slopes$trait == "log.abund" & all.slopes$species == "PHAM" & all.slopes$age == "young" ,4] = all.slopes[all.slopes$trait == "log.abund" & all.slopes$species == "PHAM" & all.slopes$age == "mature" ,4]
# all.slopes
rm(abund.pham)

#regional analysis with PHRI
abund.phri <- lme(log.abund ~ region*age , random=~1|pop, method = "REML", data=droplevels(chem[chem$region!="temperate",]), weights=varIdent(form=~1|pop))
lme_output(abund.phri) #significant region*age
#I wouldn't say residuals are fan-shaped without the pop weights, but the assumption-fitting does look a little nicer with those weights
(all.stats <- extract_stats("log.abund","PHRI", abund.phri))
#tukey test for interaction
(all.tukey <- extract_param("log.abund","PHRI",abund.phri,0.4))

#PLOT.


#--------------could make a function(dataframe, yvar, ylabel, plot_segment_constant, final export name)
plot_interactions <- function(pop_df, yvar, ){
fac = ggplot(pop_df,aes(x=lat, y=log.abund_mean,group=region))+  
  geom_point(size=3,stroke=1,aes(shape = age, color=region)) + 
  geom_errorbar(aes(ymin=log.abund_mean-log.abund_se, ymax=log.abund_mean+log.abund_se),width=.2,size=.4) +
  scale_shape_manual(values=c(19,1),name="Leaf age") + 
  scale_colour_manual(values=c("steelblue1","navyblue","grey","maroon2")) +
  geom_text(data = all.tukey[all.tukey$trait == "log.abund",], aes(x = lat.middle, y = letter_y, label = letter), size = 6) +
  geom_segment(data = all.tukey[all.tukey$trait == "log.abund",], aes(x = lat.min, xend = lat.max, y = letter_y-.1, yend = letter_y-.1)) +
  facet_grid(.~species,scales="free_x",switch="x") +
  force_panelsizes(cols = c(0.3, 1), rows = c(1.3), respect = TRUE) +
  facetted_pos_scales(x = list(scale_x_continuous(breaks = c(9, 10), limits = c(8.4,10.7))),NULL) +
  geom_abline(data=all.slopes[all.slopes$trait == "log.abund",], aes(slope=m,intercept=b,linetype = age)) +
  scale_linetype_manual(values=c("solid", "dashed"),name="Leaf age") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.spacing = unit(0.2, "lines"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 20), strip.background = element_blank(),strip.placement="outside",legend.position = "none") + 
  ylab("LC/MS peak log-abundance") + 
  xlab(expression("Latitude (\u00B0N)"))
fac

# setEPS()
# postscript("FiguresTables/ChemAbund_Lat.eps", height = 5, width = 5)
# fac
# dev.off()
}




fac = ggplot(chem_pop,aes(x=lat, y=yval,group=region))+  
  geom_point(size=3,stroke=1,aes(shape = age, color=region)) + 
  geom_errorbar(aes(ymin=log.abund_mean-log.abund_se, ymax=log.abund_mean+log.abund_se),width=.2,size=.4) +
  scale_shape_manual(values=c(19,1),name="Leaf age") + 
  scale_colour_manual(values=c("steelblue1","navyblue","grey","maroon2")) +
  geom_text(data = all.tukey[all.tukey$trait == "log.abund",], aes(x = lat.middle, y = letter_y, label = letter), size = 6) +
  geom_segment(data = all.tukey[all.tukey$trait == "log.abund",], aes(x = lat.min, xend = lat.max, y = letter_y-.1, yend = letter_y-.1)) +
  facet_grid(.~species,scales="free_x",switch="x") +
  force_panelsizes(cols = c(0.3, 1), rows = c(1.3), respect = TRUE) +
  facetted_pos_scales(x = list(scale_x_continuous(breaks = c(9, 10), limits = c(8.4,10.7))),NULL) +
  geom_abline(data=all.slopes[all.slopes$trait == "log.abund",], aes(slope=m,intercept=b,linetype = age)) +
  scale_linetype_manual(values=c("solid", "dashed"),name="Leaf age") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.spacing = unit(0.2, "lines"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 20), strip.background = element_blank(),strip.placement="outside",legend.position = "none") + 
  ylab("LC/MS peak log-abundance") + 
  xlab(expression("Latitude (\u00B0N)"))
fac

# setEPS()
# postscript("FiguresTables/ChemAbund_Lat.eps", height = 5, width = 5)
# fac
# dev.off()
