library(ggplot2)
library(tidyverse)
library(varhandle)
library(cowplot)

###########################
#set up some functions to save and plot results
###########################
#will plot sorted importance scores and then output a vector with the number of features, model R^2, and SE of R^2
imp_plot <- function(age_oldnum)
{
  #example age_oldnum is "young_1447"
  imp_file = paste("RF_R/",paste(age_oldnum,"_imp",sep=""),sep="")
  results_file = paste("RF_R/",paste(age_oldnum,"_results.txt",sep=""),sep="")
  imp <- read.csv(imp_file, sep = "\t", header=T)
  plot(imp$mean_imp ~ rownames(imp))
  results <- read.csv(results_file, sep = "\t", header=T)
  (r <- c(nrow(imp),as.numeric(results[24,1]),as.numeric(results[26,1])))
}

#will make a new list of features to use in the next model, based on a cutoff of the importance file from last model
new_list <- function(age_oldnum,age_newnum)
{
  imp_file = paste("RF_R/",age_oldnum,"_imp",sep="")
  cutoff = as.numeric(strsplit(age_newnum,split="_")[[1]][2])
  tag = paste("RF_R/",paste(age_newnum,"_features",sep=""),sep="")
  imp <- read.csv(imp_file, sep = "\t", header=T)
  name <- as.character(paste(tag,".txt",sep=""))
  write.table(imp[1:cutoff,1], name, sep = "\t", row.names = FALSE, quote=FALSE, col.names=FALSE )
}

# # sends dataset and feature list to python to run random forest. See documentation here
# https://github.com/azodichr/ML-Pipeline
#df = data frame as tab-delimited with NAs removed
#alg = algorithm (RF = random forest)
#y_name = response
#?? optional drop_na = drop rows? with NA, T/F
#gs = perform grid search, T/F (?)
#cv = cross-validation (e.g. 5-fold learns on 80% and tests on 20% of the data)
#n = 100; how often model is trained
#tag = optional, will appear at end of output file name
#feat = feature list as text file

#here's the call. takes a long time!
#system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_max -feat ./RF_R/young_max_features.txt")

new_model <- function(age_newnum)
{
  call = paste("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name log.area -gs T -cv 5 -n 100 -save ./RF_R/",paste(age_newnum,paste(" -feat ./RF_R/",paste(age_newnum,"_features.txt",sep=""),sep=""),sep=""),sep="")
  system(call)
}

###########################
# import all traits and subset to young leaves, list factors in model
###########################
#read all traits data for young leaves, each maternal line is a row. leave out extra palatability data and NMDS
young <- read.csv("Processing/2_out_AllTraits.csv",header=T) %>% filter(age == "young") %>%  select(-c(line_age,age,surv,NMDS1,NMDS2,percent_C))
#change region and species to dummy variables. not sure why I had trouble doing this in the same line as initially creating young
young <- bind_cols(             bind_cols(young[1:12],      
                                bind_cols(as.data.frame(to.dummy(young$region, "reg")), 
                                          as.data.frame(to.dummy(young$species,"spp"))     )    ),     
                      young[13:ncol(young)]) %>% 
          select(-c(region,species))
#check that all lines are unique (this is the unique ID now)
length(unique(young$line)) == nrow(young)
#check where the NAs are and how many samples are lost
sum(is.na(young$tough))
sum(is.na(young$C_N))
sum(is.na(young$area))
sum(is.na(young$richness))
nrow(young) - nrow(na.omit(young))
#4 tough, 6 C:N, 9 area, 3 for LC/MS stuff. 13 rows are lost when NAs are removed. So be it!

#write the dataset as tab-delimited for python to use, removing NAs
write.table(na.omit(young), "RF_R/RF_young_tab.csv",row.names=F,sep="\t")

#create feature list, leaving out log.area (the response) and pop and line
write(colnames(young[c(3:6,8:ncol(young))]),"RF_R/young_max_features.txt")


#########################
# iterative model runs. could probably make a for loop but it takes many hours (maybe 18?) to run everything, so not great for testing loops!
#########################
# For each step, need to change 3 things: the larger model is the first in new_list, and the smaller model is the second in new_list and the only in new_model

#steps: reduce by one-third with each step
steps = c(1929,1286,857,571,381,254,169,113,75,50,33,22,15,10)

new_model("young_max")

#make a new list
new_list("young_max","young_1286")
#run next set
new_model("young_1286")

new_list("young_1286","young_857")
new_model("young_857")

new_list("young_857","young_571")
new_model("young_571")

new_list("young_571","young_381")
new_model("young_381")

new_list("young_381","young_254")
new_model("young_254")

new_list("young_254","young_169")
new_model("young_169")

new_list("young_169","young_113")
new_model("young_113")

new_list("young_113","young_75")
new_model("young_75")

new_list("young_75","young_50")
new_model("young_50")

new_list("young_50","young_33")
new_model("young_33")

new_list("young_33","young_22")
new_model("young_22")

new_list("young_22","young_15")
new_model("young_15")

new_list("young_15","young_10")
new_model("young_10")

#########################
# visualize model fit: how R^2 changes with the number of features
#########################

#this adds results to a dataframe called feat_plot and shows a plot of the importance scores
feat_plot <- data.frame(matrix(NA, nrow = 1, ncol = 3))
names(feat_plot) <- c("features","R^2","SE")
steps_course <- c("max",steps[2:length(steps)])
for(i in 1:length(steps_course)){
  feat_plot <- rbind(feat_plot,imp_plot(paste("young_",steps_course[i],sep="")))
}
feat_plot <- na.omit(feat_plot)
feat_plot
plot(`R^2` ~ features, data = feat_plot, xlab = "Number of features", ylab = expression("Model "~R^2), main = "A) Feature selection: 1928 to 10 features", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)

#export plot
setEPS()
postscript("FiguresTables/Fig_RF_young_A.eps")
par(mar=c(5,5,4,2))
plot(`R^2` ~ features, data = feat_plot, xlab = "Number of features", ylab = expression("Model "~R^2), main = "A) Feature selection: 1929 to 10 features", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
dev.off()

#########################
# iterative model runs exploring different numbers of features around the peak in R2
#########################
#steps: reduce by five features with each step
steps_fine = c(50,45,40,35,30,25,24,23,22,21,20)

new_list("young_50","youngfine_45")
new_model("youngfine_45")

for(i in 2:(length(steps_fine)-1)){
  new_list(paste("youngfine_",steps_fine[i],sep=""),paste("youngfine_",steps_fine[i+1],sep=""))
  new_model(paste("youngfine_",steps_fine[i+1],sep=""))
}


#########################
# Results and figures for fine-scale exploration
#########################
#this adds results to a dataframe called feat_plot_fine and shows a plot of the importance scores
feat_plot_fine <- data.frame(matrix(NA, nrow = 1, ncol = 3))
names(feat_plot_fine) <- c("features","R^2","SE")
for(i in 2:length(steps_fine)){
  feat_plot_fine <- rbind(feat_plot_fine,imp_plot(paste("youngfine_",steps_fine[i],sep="")))
}
feat_plot_fine <- na.omit(feat_plot_fine)
feat_plot_fine
plot(`R^2` ~ features, data = feat_plot_fine, xlab = "Number of features", ylab = expression("Model "~R^2), main = "B) Feature selection: 45 to 20 features", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)

#export plot
setEPS()
postscript("FiguresTables/Fig_RF_young_B.eps")
par(mar=c(5,5,4,2))
plot(`R^2` ~ features, data = feat_plot_fine, xlab = "Number of features", ylab = expression("Model "~R^2), main = "B) Feature selection: 45 to 20 features", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
dev.off()

# #actual v predicted plot for final model
setEPS()
postscript("FiguresTables/Fig_RF_young_C.eps")
par(mar=c(5,5,4,2))
plot(Y ~ Mean, data = read.csv("RF_R/youngfine_24_scores.txt", sep = "\t", header=TRUE), xlab = "Predicted values", ylab = "Actual values", main = "C) Fit of best model", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
dev.off()


#########################
# Calculate percent variance explained for each feature in final model using LOFO method (leave one feature out)
#########################


 ##################
 # analyze palatability as a function of pop alone and then pop + chemicals. I probably could have done these wrangling steps at the very beginning but I'm not going back to fix it now after running all those models!
 ##################
 
 #create new df from young with population as a dummy variable (0/1 for each population)
 popdummy <- bind_cols(bind_cols(young[2:3],as.data.frame(to.dummy(young$pop,"pop"))),young[4:ncol(young)])
 
 #write tab-delimited dataset for python, removing NA
 write.table(na.omit(popdummy), "RF_R/RF_young_popdummy_tab.csv",row.names=F,sep="\t")
 
 #feature list for model with population only
 poplist <- colnames(popdummy[c(3:18)])
 write(poplist,"RF_R/young_poponly_features.txt")
 system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_popdummy_tab.csv -alg RF -y_name log.area -gs T -cv 5 -n 100 -save ./RF_R/young_pop -feat ./RF_R/young_poponly_features.txt")
 
 #feature list for model with population and the features from the best model
 chemlist <- read.csv("RF_R/youngfine_24_features.txt",sep="\t",header=F)$V1
 write(c(poplist,chemlist),"RF_R/young_popchem_features.txt")
 system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_popdummy_tab.csv -alg RF -y_name log.area -gs T -cv 5 -n 100 -save ./RF_R/young_popchem -feat ./RF_R/young_popchem_features.txt")
 
 #################
 #combine three figures into one panel (to go to the right of young leaves, letters are B,D,F)
 
 plotA = ggplot(data=feat_plot,aes(x=features, y=`R^2`))+  
   geom_point(size=3, shape = 1) +
   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 16), panel.background = element_rect(colour=NA, fill = "transparent"), plot.background = element_rect(colour=NA, fill = "transparent"), plot.margin = unit(c(0, 0.3, 0, 0), "in")) +
   ylab(expression("Model "~R^2)) + 
   xlab("Number of predictors")
 
 plotC = ggplot(data=feat_plot_fine,aes(x=features, y=`R^2`))+  
   geom_point(size=3, shape = 1) +
   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 16), panel.background = element_rect(colour=NA, fill = "transparent"), plot.background = element_rect(colour=NA, fill = "transparent")) +
   ylab(expression("Model "~R^2)) + 
   xlab("Number of predictors")
 
 plotE = ggplot(data=read.csv("RF_R/youngfine_24_scores.txt", sep = "\t", header=TRUE),aes(x=Mean, y=Y))+  
   geom_point(size=3, shape = 1) +
   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 16), panel.background = element_rect(colour=NA, fill = "transparent"), plot.background = element_rect(colour=NA, fill = "transparent")) +
   ylab("Actual values") + 
   xlab("Predicted values")
 
 
 multipanel <- plot_grid(plotA,plotC,plotE, nrow = 3, ncol = 1, align = "hv", labels = c("A","C","E"))
 save_plot("FiguresTables/Fig_RF_young_combined.pdf", multipanel, ncol = 1, nrow =3, base_height = 4, base_width = 4)
