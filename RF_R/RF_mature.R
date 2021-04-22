library(ggplot2)
library(tidyverse)
library(varhandle)

###########################
#set up some functions to save and plot results
###########################
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
#system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_max -feat ./RF_R/mature_max_features.txt")

new_model <- function(age_newnum)
{
  call = paste("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/",paste(age_newnum,paste(" -feat ./RF_R/",paste(age_newnum,"_features.txt",sep=""),sep=""),sep=""),sep="")
  system(call)
}

#will plot sorted importance scores and then output a vector with the number of features, model R^2, and SE of R^2
imp_plot <- function(age_oldnum)
{
  #example age_oldnum is "mature_1447"
  imp_file = paste("RF_R/",paste(age_oldnum,"_imp",sep=""),sep="")
  results_file = paste("RF_R/",paste(age_oldnum,"_results.txt",sep=""),sep="")
  imp <- read.csv(imp_file, sep = "\t", header=T)
  plot(imp$mean_imp ~ rownames(imp))
  results <- read.csv(results_file, sep = "\t", header=T)
  (r <- c(nrow(imp),as.numeric(results[24,1]),as.numeric(results[26,1])))
}

###########################
# import all traits and subset to mature leaves, list factors in model
###########################
#read all traits data for mature leaves, each maternal line is a row. leave out extra palatability data and NMDS
mature <- read.csv("Processing/2_out_AllTraits.csv",header=T) %>% filter(age == "mature") %>%  select(-c(line_age,age,initial,surv,mass_surv,mass_cup,NMDS1,NMDS2,percent_C))
#change region to dummy variable. not sure why I had trouble doing this in the same line as initially creating mature
mature <- bind_cols(bind_cols(mature[1:11],as.data.frame(to.dummy(mature$region, "reg"))),mature[12:ncol(mature)]) %>% select(-region)
#check that all lines are unique (this is the unique ID now)
length(unique(mature$line)) == nrow(mature)
#check where the NAs are and how many samples are lost
sum(is.na(mature$tough))
sum(is.na(mature$C_N))
sum(is.na(mature$area))
sum(is.na(mature$richness))
nrow(mature) - nrow(na.omit(mature))
#4 tough, 6 C:N, 9 area, 3 for LC/MS stuff. 13 rows are lost when NAs are removed. So be it!

#write the dataset as tab-delimited for python to use, removing NAs
write.table(na.omit(mature), "RF_R/RF_mature_tab.csv",row.names=F,sep="\t")

#create feature list, leaving out area (the response) and pop and line
write(colnames(mature[c(3:6,8:ncol(mature))]),"RF_R/mature_max_features.txt")

new_model("mature_max")

#########################
# iterative model runs. could probably make a for loop but it takes many hours (maybe 18?) to run everything, so not great for testing loops!
#########################
#steps: reduce by one-third with each step
steps = c(1929,1286,857,571,381,254,169,113,75,50,33,22,15,10)

#make a new list
new_list("mature_max","mature_1286")
#run next set
new_model("mature_1286")

new_list("mature_1286","mature_857")
new_model("mature_857")

new_list("mature_857","mature_571")
new_model("mature_571")

new_list("mature_571","mature_381")
new_model("mature_381")

new_list("mature_381","mature_254")
new_model("mature_254")

new_list("mature_254","mature_169")
new_model("mature_169")

new_list("mature_169","mature_113")
new_model("mature_113")

new_list("mature_113","mature_75")
new_model("mature_75")

new_list("mature_75","mature_50")
new_model("mature_50")

new_list("mature_50","mature_33")
new_model("mature_33")

new_list("mature_33","mature_22")
new_model("mature_22")

new_list("mature_22","mature_15")
new_model("mature_15")

new_list("mature_15","mature_10")
new_model("mature_10")

#########################
# visualize model fit: how R^2 changes with the number of features
#########################

#this adds results to a dataframe called feat_plot and shows a plot of the importance scores
feat_plot <- data.frame(matrix(NA, nrow = 1, ncol = 3))
names(feat_plot) <- c("features","R^2","SE")
steps_course <- c("max",steps[2:length(steps)])
for(i in 1:length(steps_course)){
  feat_plot <- rbind(feat_plot,imp_plot(paste("mature_",steps_course[i],sep="")))
}
feat_plot <- na.omit(feat_plot)
feat_plot
plot(`R^2` ~ features, data = feat_plot, xlab = "Number of features", ylab = expression("Model "~R^2), main = "A) Feature selection: 1928 to 10 features", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)

#export plot
setEPS()
postscript("FiguresTables/Fig_RF_mature_A.eps")
par(mar=c(5,5,4,2))
plot(`R^2` ~ features, data = feat_plot, xlab = "Number of features", ylab = expression("Model "~R^2), main = "A) Feature selection: 1928 to 10 features", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
dev.off()

#########################
# iterative model runs exploring different numbers of features around the peak in R2
#########################
#steps: reduce by one-third with each step
steps_fine = c(22,21,20,19,18,17,16,15,14,13,12,11,10)

new_list("mature_50","maturefine_40")
new_model("maturefine_40")

for(i in 2:(length(steps_fine)-1)){
  new_list(paste("maturefine_",steps_fine[i],sep=""),paste("maturefine_",steps_fine[i+1],sep=""))
  new_model(paste("maturefine_",steps_fine[i+1],sep=""))
}


#########################
# Results and figures for fine-scale exploration
#########################
#this adds results to a dataframe called feat_plot_fine and shows a plot of the importance scores
feat_plot_fine <- data.frame(matrix(NA, nrow = 1, ncol = 3))
names(feat_plot_fine) <- c("features","R^2","SE")
for(i in 2:length(steps_fine)){
  feat_plot_fine <- rbind(feat_plot_fine,imp_plot(paste("maturefine_",steps_fine[i],sep="")))
}
feat_plot_fine <- na.omit(feat_plot_fine)
feat_plot_fine
plot(`R^2` ~ features, data = feat_plot_fine, xlab = "Number of features", ylab = expression("Model "~R^2), main = "B) Feature selection: 40 to 10 features", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)

#export plot
setEPS()
postscript("FiguresTables/Fig_RF_mature_B.eps")
par(mar=c(5,5,4,2))
plot(`R^2` ~ features, data = feat_plot_fine, xlab = "Number of features", ylab = expression("Model "~R^2), main = "B) Feature selection: 40 to 10 features", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
dev.off()

# #actual v predicted plot for final model
setEPS()
postscript("FiguresTables/Fig_RF_mature_C.eps")
par(mar=c(5,5,4,2))
plot(Y ~ Mean, data = read.csv("RF_R/maturefine_XX_scores.txt", sep = "\t", header=TRUE), xlab = "Predicted values", ylab = "Actual values", main = "C) Fit of best model", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
dev.off()








##################
# analyze palatability as a function of pop alone and then pop + chemicals. I probably could have done these wrangling steps at the very beginning but I'm not going back to fix it now after running all those models!
##################
#I checked, but it only increases the dataset by 1 line to try and use a larger dataset without toughness and carbon. Same for mature.

#create new df from mature with population as a dummy variable (0/1 for each population)
popdummy <- bind_cols(bind_cols(mature[2:3],as.data.frame(to.dummy(mature$pop,"pop"))),mature[4:ncol(mature)])

#write tab-delimited dataset for python, removing NA
write.table(na.omit(popdummy), "RF_R/RF_mature_popdummy_tab.csv",row.names=F,sep="\t")

#feature list for model with population only
poplist <- colnames(popdummy[c(3:18)])
write(poplist,"RF_R/mature_poponly_features.txt")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_popdummy_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_pop -feat ./RF_R/mature_poponly_features.txt")

#feature list for model with population and the features from the best model
chemlist <- read.csv("RF_R/mature_21_features.txt",sep="\t",header=F)$V1
write(c(poplist,chemlist),"RF_R/mature_popchem_features.txt")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_popdummy_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_popchem -feat ./RF_R/mature_popchem_features.txt")


