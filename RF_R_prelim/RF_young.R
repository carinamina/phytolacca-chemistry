library(ggplot2)
library(tidyverse)
library(varhandle)

###########################
#set up some functions to save and plot results
###########################
#will plot sorted importance scores and then output a vector with the number of features, model R^2, and SE of R^2
imp_plot <- function(imp_file,results_file)
{
  imp <- read.csv(imp_file, sep = "\t", header=T)
  plot(imp$mean_imp ~ rownames(imp))
  results <- read.csv(results_file, sep = "\t", header=T)
  (r <- c(nrow(imp),as.numeric(results[24,1]),as.numeric(results[26,1])))
}

#will make a new list of features to use in the next model, based on a cutoff of the importance file from last model
new_list <- function(imp_file,cutoff,tag)
{
  imp <- read.csv(imp_file, sep = "\t", header=T)
  name <- as.character(paste(tag,".txt",sep=""))
  write.table(imp[1:cutoff,1], name, sep = "\t", row.names = FALSE, quote=FALSE, col.names=FALSE )
}

feat_plot <- data.frame(matrix(NA, nrow = 1, ncol = 3))
names(feat_plot) <- c("features","R^2","SE")

###########################
# import all traits and subset to young leaves, list factors in model
###########################
#read all traits data for young leaves, each maternal line is a row. leave out extra palatability data and NMDS
young <- read.csv("Processing/2_out_AllTraits.csv",header=T) %>% filter(age == "young") %>%  select(-c(line_age,age,initial,surv,mass_surv,mass_cup,NMDS1,NMDS2))
#change region to dummy variable. not sure why I had trouble doing this in the same line as initially creating young
young <- bind_cols(bind_cols(young[1:12],as.data.frame(to.dummy(young$region, "reg"))),young[13:ncol(young)]) %>% select(-region)
#check that all lines are unique (this is the unique ID now)
length(unique(young$line)) == nrow(young)
#check where the NAs are and how many samples are lost
sum(is.na(young$tough))
nrow(young) - nrow(na.omit(young))
#4 tough, 6 C:N, 9 area, 3 for LC/MS stuff. 13 rows are lost when NAs are removed. So be it!

#write the dataset as tab-delimited for python to use, removing NAs
write.table(na.omit(young), "RF_R/RF_young_tab.csv",row.names=F,sep="\t")

#create feature list, leaving out area (the response) and pop and line
write(colnames(young[c(3:7,9:ncol(young))]),"RF_R/young_max_features.txt")

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
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_max -feat ./RF_R/young_max_features.txt")

#troubleshoot version
#quick check with a small set of features, just leaving out area
#system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 2 -n 3 -save ./RF_R/testsave -feat ./RF_R/test_features.txt")

#########################
# visualizing results and removing features for next iteration; running next iteration
#########################
# For each step, need to change 7 things: the inputs for feat_plot (2) and new_list (3) and system (2). The larger, previous number goes in the first three slots, and the smaller new number goes in the last 4.

#this adds results to a dataframe called feat_plot and shows a plot of the importance scores
(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_max_imp","RF_R/young_max_results.txt")))
feat_plot <- na.omit(feat_plot)
#how many should be kept in the list; I'll cut by 25% each time until getting close to the end
new_list("RF_R/young_max_imp",1447,"RF_R/young_1447_features")

#run next set; always need to modify the number in both "save" and "feat"
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_1447 -feat ./RF_R/young_1447_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_1447_imp","RF_R/young_1447_results.txt")))
new_list("RF_R/young_1447_imp",1086,"RF_R/young_1086_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_1086 -feat ./RF_R/young_1086_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_1086_imp","RF_R/young_1086_results.txt")))
new_list("RF_R/young_1086_imp",815,"RF_R/young_815_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_815 -feat ./RF_R/young_815_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_815_imp","RF_R/young_815_results.txt")))
new_list("RF_R/young_815_imp",612,"RF_R/young_612_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_612 -feat ./RF_R/young_612_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_612_imp","RF_R/young_612_results.txt")))
new_list("RF_R/young_612_imp",459,"RF_R/young_459_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_459 -feat ./RF_R/young_459_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_459_imp","RF_R/young_459_results.txt")))
new_list("RF_R/young_459_imp",345,"RF_R/young_345_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_345 -feat ./RF_R/young_345_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_345_imp","RF_R/young_345_results.txt")))
new_list("RF_R/young_345_imp",259,"RF_R/young_259_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_259 -feat ./RF_R/young_259_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_259_imp","RF_R/young_259_results.txt")))
new_list("RF_R/young_259_imp",195,"RF_R/young_195_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_195 -feat ./RF_R/young_195_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_195_imp","RF_R/young_195_results.txt")))
new_list("RF_R/young_195_imp",147,"RF_R/young_147_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_147 -feat ./RF_R/young_147_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_147_imp","RF_R/young_147_results.txt")))
new_list("RF_R/young_147_imp",111,"RF_R/young_111_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_111 -feat ./RF_R/young_111_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_111_imp","RF_R/young_111_results.txt")))
new_list("RF_R/young_111_imp",84,"RF_R/young_84_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_84 -feat ./RF_R/young_84_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_84_imp","RF_R/young_84_results.txt")))
new_list("RF_R/young_84_imp",63,"RF_R/young_63_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_63 -feat ./RF_R/young_63_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_63_imp","RF_R/young_63_results.txt")))
new_list("RF_R/young_63_imp",48,"RF_R/young_48_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_48 -feat ./RF_R/young_48_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_48_imp","RF_R/young_48_results.txt")))
new_list("RF_R/young_48_imp",36,"RF_R/young_36_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_36 -feat ./RF_R/young_36_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_36_imp","RF_R/young_36_results.txt")))
new_list("RF_R/young_36_imp",27,"RF_R/young_27_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_27 -feat ./RF_R/young_27_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_27_imp","RF_R/young_27_results.txt")))
new_list("RF_R/young_27_imp",21,"RF_R/young_21_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_21 -feat ./RF_R/young_21_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_21_imp","RF_R/young_21_results.txt")))
new_list("RF_R/young_21_imp",16,"RF_R/young_16_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_16 -feat ./RF_R/young_16_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_16_imp","RF_R/young_16_results.txt")))
new_list("RF_R/young_16_imp",12,"RF_R/young_12_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_12 -feat ./RF_R/young_12_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_12_imp","RF_R/young_12_results.txt")))
new_list("RF_R/young_12_imp",9,"RF_R/young_9_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_9 -feat ./RF_R/young_9_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_9_imp","RF_R/young_9_results.txt")))
new_list("RF_R/young_9_imp",7,"RF_R/young_7_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_7 -feat ./RF_R/young_7_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_7_imp","RF_R/young_7_results.txt")))
new_list("RF_R/young_7_imp",6,"RF_R/young_6_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_6 -feat ./RF_R/young_6_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_6_imp","RF_R/young_6_results.txt")))
new_list("RF_R/young_6_imp",5,"RF_R/young_5_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_5 -feat ./RF_R/young_5_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_5_imp","RF_R/young_5_results.txt")))
new_list("RF_R/young_5_imp",4,"RF_R/young_4_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_4 -feat ./RF_R/young_4_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_4_imp","RF_R/young_4_results.txt")))
new_list("RF_R/young_4_imp",3,"RF_R/young_3_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_3 -feat ./RF_R/young_3_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/young_3_imp","RF_R/young_3_results.txt")))

# # export an importance plot of interest
imp <- read.csv("RF_R/young_63_imp", sep = "\t", header=T)
setEPS()
postscript("FiguresTables/FigS2A_young_imp.eps")
par(mar=c(5,5,4,2))
plot(mean_imp~rownames(imp), data = imp, xlab = "Feature rank", ylab = "Feature importance score (Gini index)", main = "A) young-leaf palatability ~ 63 features", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
dev.off()
rm(imp)


#########################
# visualize model fit: how R^2 changes with the number of features, actual vs predicted plot
#########################

feat_plot
plot(`R^2` ~ features, data = feat_plot, xlab = "Number of features", ylab = expression("Model "~R^2), main = "B) Feature selection", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)

#export plot
setEPS()
postscript("FiguresTables/FigS2B_young_R2.eps")
par(mar=c(5,5,4,2))
plot(`R^2` ~ features, data = feat_plot, xlab = "Number of features", ylab = expression("Model "~R^2), main = "B) Feature selection", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
dev.off()


#I checked, but it only increases the dataset by 1 line to try and use a larger dataset without toughness and carbon. Same for mature.

# #actual v predicted plot for final model
 setEPS()
 postscript("FiguresTables/FigS2C_young_fit.eps")
 par(mar=c(5,5,4,2))
plot(Y ~ Mean, data = read.csv("RF_R/young_48_scores.txt", sep = "\t", header=TRUE), xlab = "Predicted values", ylab = "Actual values", main = "C) Fit of best model", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
 dev.off()

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
 system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_popdummy_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_pop -feat ./RF_R/young_poponly_features.txt")
 
 #feature list for model with population and the features from the best model
 chemlist <- read.csv("RF_R/young_21_features.txt",sep="\t",header=F)$V1
 write(c(poplist,chemlist),"RF_R/young_popchem_features.txt")
 system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_young_popdummy_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/young_popchem -feat ./RF_R/young_popchem_features.txt")
 
 
