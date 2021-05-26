#analysis and figures for analyzing latitudinal distributions of leaf traits
#started 18 May 2021, modeled after palat_lat and 2018 scripts traits_lat_pham and traits_lat_phri


#libraries
###############
library(nlme)
library(multcomp)
library(ggplot2)
library(tidyverse)
library(ggh4x)
library(cowplot)

##############
# Notes and structure of analysis
##############
#three datasets for the six traits of interest: toughness, %N, C:N, and chemical abund, richness, and diversity. Chemistry was measured on one individual per line, the others on multiple. For each trait, for a latitude and regional model, we need to make sure the nesting is correct, check assumptions to decide if weights of pop or pop_age are needed, extract model results to a big table, and then extract slopes and/or Tukey comparisons for plotting if significant. Then calculate population and regional means for plotting.
#I'll try and write functions for everything to minimize pasting big blocks.

###############
#import data, check out boxplots
###############
chem <- read.csv("Processing/1b_out_ChemSummaries_Indiv.csv", header=T) %>% select(-c(NMDS1,NMDS2))
length(unique(chem$line_age)) == nrow(chem)
chem <- reshape(chem, varying = c("log.abund","richness","diversity"), times = c("log.abund","richness","diversity"), v.names = "value", timevar = "trait", direction = "long") %>% select(-id) %>% remove_rownames() %>% select(trait,pos,age,pos_age,line,line_age,pop,lat,region,value)

tough <- read.csv("Processing/1c_out_Toughness_Indiv.csv", header=T)  %>% rename(value = tough) %>% add_column(trait = "tough") %>% select(trait,pos,age,pos_age,line,line_age,pop,lat,region,value)
length(unique(tough$line_age)) == nrow(tough) #multiple indiv per line
length(unique(tough$pos_age)) == nrow(tough)

cn <- read.csv("Processing/1d_out_CarbonNitrogen_Indiv.csv", header=T) %>% select(-c(percent_C))
length(unique(cn$line_age)) == nrow(cn) #multiple indiv per line
length(unique(cn$pos_age)) == nrow(cn)
cn <- reshape(cn, varying = c("percent_N","C_N"), times = c("percent_N","C_N"), v.names = "value", timevar = "trait", direction = "long") %>% select(-id) %>% remove_rownames() %>% select(trait,pos,age,pos_age,line,line_age,pop,lat,region,value)

#combine trait datasets and tweak some column formats
names(chem) == names(tough)
names(chem) == names(cn)
all.traits <- na.omit(bind_rows(chem,bind_rows(tough,cn)) %>% mutate(species = ifelse(region == "tropical","PHRI","PHAM"), .before="value"))
all.traits$pop_age <- as.factor(paste(all.traits$pop, all.traits$age, sep = "_"))
all.traits$reg_age <- interaction(all.traits$region, all.traits$age)
#switches the order of the species, which is needed way down the line for plotting PHRI on the left of PHAM
all.traits$species = factor(all.traits$species, levels=c("PHRI","PHAM"))
all.traits$age <- as.factor(all.traits$age)
all.traits$region <- as.factor(all.traits$region)
rm(chem,cn,tough)

#get trait list vector
trait_list <- unique(all.traits$trait)

#check for heterogeneity of variance in each trait
par(mfrow=c(2,3))
#young
for (i in 1:length(trait_list)){
  plotdata.y <- subset(all.traits, all.traits$age == "young" & all.traits$trait == trait_list[i])
  boxplot(value~pop,
          data=plotdata.y,
          main=paste(trait_list[i],"by pop_young"),
          xlab="population" ,
          ylab = trait_list[i])
  rm(plotdata.y,i)
}
for (i in 1:length(trait_list)){
  plotdata.y <- subset(all.traits, all.traits$age == "mature" & all.traits$trait == trait_list[i])
  boxplot(value~pop,
          data=plotdata.y,
          main=paste(trait_list[i],"by pop_mature"),
          xlab="population" ,
          ylab = trait_list[i])
  rm(plotdata.y,i)
}
#percent_N, C_N, diversity, log.abund look like there will be problems with heterogeneity of variance. tough, richness might be ok.
par(mfrow=c(1,1))


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

#adds ANOVA results of a model to the summary table which has already been created. Need to call as 
#all.stats <- extract_stats(a,b,c)
extract_stats <- function(response_string, species_string, lme_model, placeholder)
{
  ano <- as.data.frame(anova(lme_model, type = "marginal", test= "F"))[2:4,] %>% rownames_to_column(var = "predictor") %>% add_column(trait = response_string, species = species_string, .before = "predictor") 
  all.stats[all.stats$trait == response_string & all.stats$species == species_string,] <- ano
  return(all.stats)
  rm(ano)
}

#creates summary table with all PHAM slopes
all.slopes <- data.frame(trait = sort(rep(trait_list,4))) %>% add_column(species = rep(c("PHRI","PHRI","PHAM","PHAM"),6), age = "", m = NA, b = NA)
all.slopes$species = factor(all.slopes$species, levels=c("PHRI","PHAM"))

#creates summary table with all the PHRI tukey comparisons. it also includes the regional means for each trait, which are helpful for calculating effect sizes and necessary for figuring out where to place letters on plots. Latitude is also needed for where to put letters and lines.
all.tukey <- as.data.frame(all.traits %>% mutate(region.age = paste(region,age,sep=".")) %>%
  group_by(trait, region.age,region,age,species) %>%
  filter(region != "temperate") %>%
  summarise(lat.min = min(lat), lat.max = max(lat), lat.middle = mean(c(lat.min,lat.max)), mean = mean(value)) %>%
  add_column(letter = "", letter_y = NA, segment_y = NA) )

#adds slopes to summary table which has already been created. For now I will keep this simple and assume a significant interaction; if not, it can be changed in all.slopes
#For PHAM, need to call as all.slopes <- extract_param(a,b,c,d) where d = 0 or NULL
#For PHRI, need to call as all.tukey <- extract_param(a,b,c,d)
extract_param <- function(response_string,species,lme_model,letter_y_constant,segment_y_constant)
{
  if(species == "PHAM"){
  slopes = data.frame(species = c("PHRI","PHRI","PHAM","PHAM"),age = c("mature","young"), m = c(0,0,NA,NA), b = c(0,0,NA,NA))
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
    interaction_model <- update(lme_model, fixed = "value ~ reg_age")
    tuk <- glht(interaction_model, linfct=mcp(reg_age="Tukey"))
    (tuk.cld <- cld(tuk))   # letters
    #convert the Tukey letter results to a dataframe for plotting letters. this seems clunky but I like being sure that I'm not just inserting letters in the larger summary data without merging somehow; they could be in the wrong order
    tukey <- as.data.frame(tuk.cld$mcletters$Letters) %>%
      `colnames<-`("letter") %>%
      rownames_to_column(var = "region.age")
    tukey <- left_join(all.tukey[all.tukey$trait == response_string,1:9] ,tukey,by="region.age")  %>%
      mutate(letter_y = mean + letter_y_constant, segment_y = letter_y - segment_y_constant)
    
    all.tukey[all.tukey$trait == response_string,] = tukey
    return(all.tukey)
    rm(tukey)
  }
  else{
    return("error. species is misspelled?")
  }
}

#calculate population means, jigger lat by 0.3 for young leaves so error bars of young and mature don't overlap
pop.means <- as.data.frame(all.traits %>% group_by(trait,pop_age,region,age,species) %>%
                             summarise(lat = mean(lat),
                                       trait.mean = mean(value),
                                       trait.sd = sd(value),
                                       trait.n = length(pop_age),
                                       trait.se = trait.sd/sqrt(trait.n)) %>%
                             mutate(lat = ifelse(age == "young",lat-0.3,lat)) %>%
                             select(-c(trait.sd,trait.n)))

#function to make standardized plot for each response
facet_plot <- function(response_string, ylabel_string){
  vv = ggplot(pop.means[pop.means$trait == response_string,],aes(x=lat, y=trait.mean,group=region))+  
    geom_point(size=3,stroke=1,aes(shape = age, color=region)) + 
    geom_errorbar(aes(ymin=trait.mean-trait.se, ymax=trait.mean+trait.se),width=.2,size=.4) +
    scale_shape_manual(values=c(19,1),name="Leaf age") + 
    scale_colour_manual(values=c("steelblue1","navyblue","grey","maroon2")) +
    geom_text(data = all.tukey[all.tukey$trait == response_string,], aes(x = lat.middle, y = letter_y, label = letter), size = 6) +
    geom_segment(data = all.tukey[all.tukey$trait == response_string,], aes(x = lat.min-0.3, xend = lat.max, y = segment_y, yend = segment_y)) +
    facet_grid(.~species,scales="free_x",switch="x") +
    force_panelsizes(cols = c(0.3, 1), rows = c(1.3), respect = TRUE) +
    facetted_pos_scales(x = list(scale_x_continuous(breaks = c(9, 10), limits = c(8.1,10.7))),NULL) +
    geom_abline(data=all.slopes[all.slopes$trait == response_string,], aes(slope=m,intercept=b,linetype = age)) +
    scale_linetype_manual(values=c("solid", "dashed"),name="Leaf age") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.spacing = unit(0.2, "lines"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 20), strip.background = element_blank(),strip.placement="outside",legend.position = "none", panel.background = element_rect(colour=NA, fill = "transparent"), plot.background = element_rect(colour=NA, fill = "transparent"), plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in")) +
    ylab(ylabel_string) + 
    xlab(expression("Latitude (\u00B0N)"))
  vv
}


#################
#################
# MODELS AND PLOTS
#################
#################

#################
#CHEMICAL ABUNDANCE
#################
#latitudinal gradient with PHAM
abund.pham <- lme(value ~ lat*age , random=~1|pop, method = "REML", data=droplevels(all.traits[all.traits$lat>20 & all.traits$trait == "log.abund",]), weights=varIdent(form=~1|pop))
#weights argument necessary to deal with fan-shaped residuals
lme_output(abund.pham) #lat and age are significant but not lat*age
(all.stats <- extract_stats("log.abund", "PHAM", abund.pham))
(all.slopes <- extract_param("log.abund", "PHAM", abund.pham))
#lat and age are significant but not their interaction. But if I use the mature slope for both, the fit looks weird because the slopes and intercepts are fit with the lat*age term in the model. Not sure yet how to deal with this properly
# all.slopes[all.slopes$trait == "log.abund" & all.slopes$species == "PHAM" & all.slopes$age == "young" ,4] = all.slopes[all.slopes$trait == "log.abund" & all.slopes$species == "PHAM" & all.slopes$age == "mature" ,4]
# all.slopes
rm(abund.pham)

#regional analysis with PHRI
abund.phri <- lme(value ~ region*age , random=~1|pop, method = "REML", data=droplevels(all.traits[all.traits$region != "temperate" & all.traits$trait == "log.abund",]), weights=varIdent(form=~1|pop))
lme_output(abund.phri) #significant region*age
#I wouldn't say residuals are fan-shaped without the pop weights, but the assumption-fitting does look a little nicer with those weights
(all.stats <- extract_stats("log.abund","PHRI", abund.phri))
#tukey test for interaction
(all.tukey <- extract_param("log.abund","PHRI",abund.phri,0.4,0.1))
rm(abund.phri)

(abund <- facet_plot("log.abund","Chemical log-abundance"))

setEPS()
postscript("FiguresTables/Fig_Geography_Abund.eps", height = 5, width = 5)
abund
dev.off()
rm(abund)

#################
#CHEMICAL DIVERSITY
#################
div.pham <- lme(value ~ lat*age , random=~1|pop, method = "REML", data=droplevels(all.traits[all.traits$lat>20 & all.traits$trait == "diversity",]), weights=varIdent(form=~1|pop))
#weights argument improves assumption-fitting
lme_output(div.pham) #significant lat*age
(all.stats <- extract_stats("diversity", "PHAM", div.pham))
(all.slopes <- extract_param("diversity", "PHAM", div.pham))
rm(div.pham)

div.phri <- lme(value ~ region + age + region:age , random=~1|pop, method = "REML", data=droplevels(all.traits[all.traits$region != "temperate" & all.traits$trait == "diversity",]), weights=varIdent(form=~1|pop))
lme_output(div.phri) #only region is significant
#I wouldn't say residuals are fan-shaped without the pop weights, but the assumption-fitting does look a little nicer with those weights. Long tail in the normal plot though
(all.stats <- extract_stats("diversity","PHRI", div.phri))
#tukey test for region only
tuk <- glht(div.phri, linfct=mcp(region="Tukey"))
(tuk.cld <- cld(tuk))
tukey <- as.data.frame(tuk.cld$mcletters$Letters) %>%
  `colnames<-`("letter") %>%
  rownames_to_column(var = "region")
tukey <- left_join(all.tukey[all.tukey$trait == "diversity",1:9] ,tukey,by="region")  %>%
  mutate(letter_y = mean + 0.25, segment_y = letter_y - 0.04) %>% filter(age == "mature")
all.tukey[all.tukey$trait == "diversity" & all.tukey$age == "mature",] = tukey
all.tukey

(diversity <- facet_plot("diversity","Chemical diversity"))
rm(tukey, tuk, tuk.cld,div.phri)

setEPS()
postscript("FiguresTables/Fig_Geography_Diversity.eps", height = 5, width = 5)
diversity
dev.off()
rm(diversity)

#################
#CHEMICAL RICHNESS
#################
richness.pham <- lme(value ~ lat*age , random=~1|pop, method = "REML", data=droplevels(all.traits[all.traits$lat>20 & all.traits$trait == "richness",]))
#weights argument does not improve assumption-fitting
lme_output(richness.pham) #nothing is significant
(all.stats <- extract_stats("richness", "PHAM", richness.pham))
rm(richness.pham)

richness.phri <- lme(value ~ region + age + region:age , random=~1|pop, method = "REML", data=droplevels(all.traits[all.traits$region != "temperate" & all.traits$trait == "richness",]))
#weights argument does not improve assumption-fitting
lme_output(richness.phri) #significant region*age
(all.stats <- extract_stats("richness","PHRI", richness.phri))
(all.tukey <- extract_param("richness","PHRI", richness.phri,150,25))
#manually move three of the lower tukey bars/letters to below instead of above the regions
all.tukey[all.tukey$trait == "richness" & (all.tukey$region.age == "tropical.mature" | all.tukey$region.age == "subtropical.young" | all.tukey$region.age == "north temperate.young"),11] = c(575,655,725)
all.tukey[all.tukey$trait == "richness" & (all.tukey$region.age == "tropical.mature" | all.tukey$region.age == "subtropical.young" | all.tukey$region.age == "north temperate.young"),12] = c(575,655,725)+25
rm(richness.phri)

(richness <- facet_plot("richness","Chemical richness"))

setEPS()
postscript("FiguresTables/Fig_Geography_Richness.eps", height = 5, width = 5)
richness
dev.off()
rm(richness)

#################
#LEAF TOUGHNESS: now we need to nest line in pop
#################
ctrl <- lmeControl(opt = 'optim')
tough.pham <- lme(value ~ lat*age , random=~1|pop/line, method = "REML", data=droplevels(all.traits[all.traits$lat>20 & all.traits$trait == "tough",]), weights=varIdent(form=~1|pop_age), control = ctrl)
#weights pop_age needed to deal with fan-shaped residuals
lme_output(tough.pham) #lat*age is sig
(all.stats <- extract_stats("tough", "PHAM", tough.pham))
(all.slopes <- extract_param("tough", "PHAM", tough.pham))
rm(tough.pham)

tough.phri <- lme(value ~ region + age + region:age , random=~1|pop/line, method = "REML", data=droplevels(all.traits[all.traits$region != "temperate" & all.traits$trait == "tough",]), weights=varIdent(form=~1|pop_age), control = ctrl)
#weights pop_age needed to deal with fan-shaped residuals
lme_output(tough.phri) #significant region*age
(all.stats <- extract_stats("tough","PHRI", tough.phri))
(all.tukey <- extract_param("tough","PHRI", tough.phri,50,10))
rm(tough.phri,ctrl)

(tough <- facet_plot("tough","Toughness (g force)"))

setEPS()
postscript("FiguresTables/Fig_Geography_Toughness.eps", height = 5, width = 5)
tough
dev.off()
rm(tough)

#################
#PERCENT NITROGEN: nest line in pop
#################
ctrl <- lmeControl(opt = 'optim')
percent_N.pham <- lme(value ~ lat*age , random=~1|pop/line, method = "REML", data=droplevels(all.traits[all.traits$lat>20 & all.traits$trait == "percent_N",]), weights=varIdent(form=~1|pop_age), control = ctrl)
#weights pop_age needed to deal with fan-shaped residuals. still a long tail in normal plot but we have some weirdly high values
lme_output(percent_N.pham) #lat*age is sig
(all.stats <- extract_stats("percent_N", "PHAM", percent_N.pham))
(all.slopes <- extract_param("percent_N", "PHAM", percent_N.pham))
rm(percent_N.pham)

percent_N.phri <- lme(value ~ region + age + region:age , random=~1|pop/line, method = "REML", data=droplevels(all.traits[all.traits$region != "temperate" & all.traits$trait == "percent_N",]), weights=varIdent(form=~1|pop))
#weights pop needed to deal with fan-shaped residuals (but not pop_age)
lme_output(percent_N.phri) #significant region*age
(all.stats <- extract_stats("percent_N","PHRI", percent_N.phri))
(all.tukey <- extract_param("percent_N","PHRI", percent_N.phri,.75,0.15))
rm(percent_N.phri,ctrl)

#manually move some of the tukey bars and letters
all.tukey[all.tukey$trait == "percent_N" & all.tukey$age == "mature" ,11] = c(1.1,0.9,1.6)
all.tukey[all.tukey$trait == "percent_N" & all.tukey$age == "mature" ,12] = c(1.1,0.9,1.6) + 0.15
all.tukey[all.tukey$trait == "percent_N" & all.tukey$region.age == "north temperate.mature",10] = "ab (M)"
all.tukey[all.tukey$trait == "percent_N" & all.tukey$region.age == "north temperate.young",10] = "ab (Y)"
all.tukey[all.tukey$trait == "percent_N" & all.tukey$region.age == "north temperate.young",11] = 3
all.tukey[all.tukey$trait == "percent_N" & all.tukey$region.age == "north temperate.young",12] = 3 - 0.15

(percent_N <- facet_plot("percent_N","% Nitrogen"))

setEPS()
postscript("FiguresTables/Fig_Geography_Nitrogen.eps", height = 5, width = 5)
percent_N
dev.off()
rm(percent_N)

#################
#CARBON:NITROGEN: nest line in pop
#################
ctrl <- lmeControl(opt = 'optim')
C_N.pham <- lme(value ~ lat*age , random=~1|pop/line, method = "REML", data=droplevels(all.traits[all.traits$lat>20 & all.traits$trait == "C_N",]))
#weights argument not needed
lme_output(C_N.pham) #lat*age is sig
(all.stats <- extract_stats("C_N", "PHAM", C_N.pham))
(all.slopes <- extract_param("C_N", "PHAM", C_N.pham))
#rm(C_N.pham)

C_N.phri <- lme(value ~ region + age + region:age , random=~1|pop/line, method = "REML", data=droplevels(all.traits[all.traits$region != "temperate" & all.traits$trait == "C_N",]), weights=varIdent(form=~1|pop))
#weights pop needed to deal with fan-shaped residuals (but not pop_age)
lme_output(C_N.phri) #significant region*age
(all.stats <- extract_stats("C_N","PHRI", C_N.phri))
(all.tukey <- extract_param("C_N","PHRI", C_N.phri,6,1))
#rm(C_N.phri)

#manually move some of the tukey bars and letters
all.tukey[all.tukey$trait == "C_N" & all.tukey$age == "young" ,11] = c(17,14,8)
all.tukey[all.tukey$trait == "C_N" & all.tukey$age == "young" ,12] = c(17,14,8) + 1
all.tukey[all.tukey$trait == "C_N" & all.tukey$region.age == "north temperate.mature",10] = "c (M)"
all.tukey[all.tukey$trait == "C_N" & all.tukey$region.age == "north temperate.young",10] = "c (Y)"

(C_N <- facet_plot("C_N","Carbon:Nitrogen"))

setEPS()
postscript("FiguresTables/Fig_Geography_CN.eps", height = 5, width = 5)
C_N
dev.off()
rm(C_N)


################
# a bit of table formatting
################
form.stats <- all.stats %>% mutate("d.f." = paste(numDF,denDF,sep=","),.before = "p-value") %>% 
  select(-c(numDF,denDF))
#formats the F and p to round and keep trailing zeros
form.stats$"F-value" <- sprintf("%.2f", round(form.stats$"F-value",2))
form.stats$"p-value" <- ifelse(form.stats$"p-value" < 0.001, "< 0.001", sprintf("%.3f", round(form.stats$"p-value",3) ) )
#format traits and predictors, sort by trait in the order I like and then by species
form.stats$predictor <- plyr::revalue(as.factor(form.stats$predictor), replace = c("lat" = "Latitude", "lat:age" = "Lat x Age", "region" = "Region", "age" = "Leaf age", "region:age" = "Reg x Age"))
form.stats$trait <- plyr::revalue(as.factor(form.stats$trait), replace = c("C_N" = "Carbon:Nitrogen", "diversity" = "Chemical diversity", "log.abund" = "Chemical abundance", "percent_N" = "% Nitrogen", "richness" = "Chemical richness", "tough" = "Leaf toughness"))
form.stats$trait = factor(form.stats$trait, levels=c("Chemical abundance","Chemical richness","Chemical diversity","Leaf toughness","% Nitrogen","Carbon:Nitrogen"))
form.stats <- form.stats %>% arrange(trait) %>% arrange(species) %>% rename(Trait = trait, Species = species, "Fixed effect" = predictor, "F" = "F-value", "p" = "p-value")

#export so I can make it wide with an extra header of species in excel (not sure how to do that here)
write.csv(form.stats, "FiguresTables/Table_Geography_Traits.csv",row.names=F)

################
# multipanel figure
###############
#remove x axis labels and just keep ticks for four
#adding a negative plot margin to the theme does get rid of some white space in individual plots (weirdly only top or bottom; if you change both to -0.5, it's like they are back to 0 again), but I tried adding it to these and building the multipanel and it doesn't change anything.
#+ theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.text.x = element_blank(), plot.margin = unit(c(-.5, 0, 0, 0), "in"))
multi_abund <- facet_plot("log.abund","Chemical log-abundance") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.text.x = element_blank()) 
multi_richness <- facet_plot("richness","Chemical richness") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.text.x = element_blank()) 
multi_diversity <- facet_plot("diversity","Chemical diversity") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.text.x = element_blank()) 
multi_tough <- facet_plot("tough","Toughness (g force)") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.text.x = element_blank()) 
multi_percent_N <- facet_plot("percent_N","% Nitrogen")
multi_C_N <- facet_plot("C_N","Carbon:Nitrogen")


multipanel <- plot_grid(multi_abund,multi_richness,multi_diversity,multi_tough,multi_percent_N,multi_C_N, nrow = 3, ncol = 2, align = "hv", labels = "AUTO")
save_plot("FiguresTables/Fig_Geography_All.pdf", multipanel, ncol = 2, nrow =3, base_height = 5, base_width = 5)

############### examples for multipanels
#https://github.com/wilkelab/cowplot/issues/31
# p1 <- qplot(1:10, 1:10) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
# p2 <- qplot(1:10, (1:10)^2) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
# plot_grid(p1, p2, p1, p2, p1, p2, p1, p2, p1, p2, p1, p2)
# 
# https://wilkelab.org/cowplot/reference/save_plot.html
# p1 <- ggplot(mpg, aes(x = cty, y = hwy, color = factor(cyl))) +
#   geom_point(size = 2) +
#   theme_half_open()
# 
# file1 <- tempfile("file1", fileext = ".png")
# file2 <- tempfile("file2", fileext = ".png")
# save_plot(file1, p1)
# # same as file1 but determine base_width given base_height
# save_plot(file2, p1, base_height = NULL, base_width = 6)
# 
# # save a single plot without legend, adjust aspect ratio
# x <- (1:100)/10
# p3 <- ggplot(data.frame(x = x, y = x*sin(x)), aes(x, y)) +
#   geom_line() +
#   theme_minimal_hgrid()
# save_plot("FiguresTables/example1.pdf", p3, base_asp = 1.1)
# 
# # now combine with a second plot and save
# p3b <- ggplot(data.frame(x = x, y = cos(x)+x), aes(x, y)) +
#   geom_line() +
#   theme_minimal_hgrid()
# p4 <- plot_grid(p3, p3b, labels = "AUTO")
# save_plot("FiguresTables/example2.pdf", p4, ncol = 2, base_asp = 1.1)


#this didn't quite work
# https://stackoverflow.com/questions/55151531/ggplot2-arrange-multiple-plots-all-the-same-size-no-gaps-in-between
# allplotslist <- align_plots(multi_abund, multi_richness, multi_diversity, multi_tough, multi_percent_N, multi_C_N, align = "hv")
# #(x co-ord, y co-ord, width, height)
# multipanel <- ggdraw() + 
#   draw_plot(allplotslist[[5]], 0,   0,   0.4, 0.4) + 
#   draw_plot(allplotslist[[6]], 0.375, 0,   0.4, 0.4) + 
#   draw_plot(allplotslist[[3]], 0,   0.31, 0.4, 0.4) +
#   draw_plot(allplotslist[[4]], 0.375, 0.31, 0.4, 0.4) + 
#   draw_plot(allplotslist[[1]], 0,   0.62, 0.4, 0.4) + 
#   draw_plot(allplotslist[[2]], 0.375, 0.62,   0.4, 0.4) 
# multipanel
# 
# save_plot("FiguresTables/Fig_Geography_All.pdf", multipanel, ncol = 2, nrow = 3, base_height = 5, base_width = 5)
# 
# setEPS()
# postscript("FiguresTables/Fig_Geography_All.eps", height = 15, width = 12)
# multipanel
# dev.off()




# #original plotting code
# fac = ggplot(pop.means[pop.means$trait == "diversity",],aes(x=lat, y=trait.mean,group=region))+
#   geom_point(size=3,stroke=1,aes(shape = age, color=region)) +
#   geom_errorbar(aes(ymin=trait.mean-trait.se, ymax=trait.mean+trait.se),width=.2,size=.4) +
#   scale_shape_manual(values=c(19,1),name="Leaf age") +
#   scale_colour_manual(values=c("steelblue1","navyblue","grey","maroon2")) +
#   geom_text(data = all.tukey[all.tukey$trait == "diversity",], aes(x = lat.middle, y = letter_y, label = letter), size = 6) +
#   geom_segment(data = all.tukey[all.tukey$trait == "diversity",], aes(x = lat.min, xend = lat.max, y = letter_y-.02, yend = letter_y-.02)) +
#   facet_grid(.~species,scales="free_x",switch="x") +
#   force_panelsizes(cols = c(0.3, 1), rows = c(1.3), respect = TRUE) +
#   facetted_pos_scales(x = list(scale_x_continuous(breaks = c(9, 10), limits = c(8.4,10.7))),NULL) +
#   geom_abline(data=all.slopes[all.slopes$trait == "diversity",], aes(slope=m,intercept=b,linetype = age)) +
#   scale_linetype_manual(values=c("solid", "dashed"),name="Leaf age") +
#   theme_bw() + theme(panel.grid.major = element_blank(), panel.spacing = unit(0.2, "lines"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 20), strip.background = element_blank(),strip.placement="outside",legend.position = "top") +
#   ylab("Chemical diversity (Shannon index)") +
#   xlab(expression("Latitude (\u00B0N)"))
# fac
# setEPS()
# postscript("FiguresTables/Fig_Geography_Key.eps", height = 8, width = 12)
# fac
# dev.off()
