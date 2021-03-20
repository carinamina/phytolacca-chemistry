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
# import all traits and subset to mature leaves, list factors in model
###########################
#read all traits data for mature leaves, each maternal line is a row. leave out extra palatability data and NMDS
mature <- read.csv("Processing/2_out_AllTraits.csv",header=T) %>% filter(age == "mature") %>%  select(-c(line_age,age,initial,surv,mass_surv,mass_cup,NMDS1,NMDS2,NMDS3,NMDS4))
#change region to dummy variable. not sure why I had trouble doing this in the same line as initially creating mature
mature <- bind_cols(bind_cols(mature[1:12],as.data.frame(to.dummy(mature$region, "reg"))),mature[13:ncol(mature)]) %>% select(-region)
#check that all lines are unique (this is the unique ID now)
length(unique(mature$line)) == nrow(mature)
#check where the NAs are and how many samples are lost
sum(is.na(mature$richness))
nrow(mature) - nrow(na.omit(mature))
#4 tough, 6 C:N, 6 area, 3 for LC/MS stuff. 10 rows are lost when NAs are removed. So be it!

#write the dataset as tab-delimited for python to use, removing NAs
write.table(na.omit(mature), "RF_R/RF_mature_tab.csv",row.names=F,sep="\t")

#create feature list, leaving out area (the response) and pop and line
write(colnames(mature[c(3:7,9:ncol(mature))]),"RF_R/mature_max_features.txt")

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
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_max -feat ./RF_R/mature_max_features.txt")

#troubleshoot version
#quick check with a small set of features, just leaving out area
#write(colnames(mature[c(3:7,9:12)]),"RF_R/RF_mature_MaxFeatures.txt")
#system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 2 -n 3 -save ./RF_R/testsave -feat ./RF_R/test_features.txt")

#########################
# visualizing results and removing features for next round
# For each step, need to change 7 things: the inputs for feat_plot (2) and new_list (3) and system (2). The larger, previous number goes in the first three slots, and the smaller new number goes in the last 4.

#this adds results to a dataframe called feat_plot and shows a plot of the importance scores
(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_max_imp","RF_R/mature_max_results.txt")))
feat_plot <- na.omit(feat_plot)
#how many should be kept in the list; I'll cut by 25% each time until getting close to the end
new_list("RF_R/mature_max_imp",1447,"RF_R/mature_1447_features")

#run next set; always need to modify the number in both "save" and "feat"
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_1447 -feat ./RF_R/mature_1447_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_1447_imp","RF_R/mature_1447_results.txt")))
new_list("RF_R/mature_1447_imp",1086,"RF_R/mature_1086_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_1086 -feat ./RF_R/mature_1086_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_1086_imp","RF_R/mature_1086_results.txt")))
new_list("RF_R/mature_1086_imp",815,"RF_R/mature_815_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_815 -feat ./RF_R/mature_815_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_815_imp","RF_R/mature_815_results.txt")))
new_list("RF_R/mature_815_imp",612,"RF_R/mature_612_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_612 -feat ./RF_R/mature_612_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_612_imp","RF_R/mature_612_results.txt")))
new_list("RF_R/mature_612_imp",459,"RF_R/mature_459_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_459 -feat ./RF_R/mature_459_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_459_imp","RF_R/mature_459_results.txt")))
new_list("RF_R/mature_459_imp",345,"RF_R/mature_345_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_345 -feat ./RF_R/mature_345_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_345_imp","RF_R/mature_345_results.txt")))
new_list("RF_R/mature_345_imp",259,"RF_R/mature_259_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_259 -feat ./RF_R/mature_259_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_259_imp","RF_R/mature_259_results.txt")))
new_list("RF_R/mature_259_imp",195,"RF_R/mature_195_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_195 -feat ./RF_R/mature_195_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_195_imp","RF_R/mature_195_results.txt")))
new_list("RF_R/mature_195_imp",147,"RF_R/mature_147_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_147 -feat ./RF_R/mature_147_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_147_imp","RF_R/mature_147_results.txt")))
new_list("RF_R/mature_147_imp",111,"RF_R/mature_111_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_111 -feat ./RF_R/mature_111_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_111_imp","RF_R/mature_111_results.txt")))
new_list("RF_R/mature_111_imp",84,"RF_R/mature_84_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_84 -feat ./RF_R/mature_84_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_84_imp","RF_R/mature_84_results.txt")))
new_list("RF_R/mature_84_imp",63,"RF_R/mature_63_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_63 -feat ./RF_R/mature_63_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_63_imp","RF_R/mature_63_results.txt")))
new_list("RF_R/mature_63_imp",48,"RF_R/mature_48_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_48 -feat ./RF_R/mature_48_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_48_imp","RF_R/mature_48_results.txt")))
new_list("RF_R/mature_48_imp",36,"RF_R/mature_36_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_36 -feat ./RF_R/mature_36_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_36_imp","RF_R/mature_36_results.txt")))
new_list("RF_R/mature_36_imp",27,"RF_R/mature_27_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_27 -feat ./RF_R/mature_27_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_27_imp","RF_R/mature_27_results.txt")))
new_list("RF_R/mature_27_imp",21,"RF_R/mature_21_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_21 -feat ./RF_R/mature_21_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_21_imp","RF_R/mature_21_results.txt")))
new_list("RF_R/mature_21_imp",16,"RF_R/mature_16_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_16 -feat ./RF_R/mature_16_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_16_imp","RF_R/mature_16_results.txt")))
new_list("RF_R/mature_16_imp",12,"RF_R/mature_12_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_12 -feat ./RF_R/mature_12_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_12_imp","RF_R/mature_12_results.txt")))
new_list("RF_R/mature_12_imp",9,"RF_R/mature_9_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_9 -feat ./RF_R/mature_9_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_9_imp","RF_R/mature_9_results.txt")))
new_list("RF_R/mature_9_imp",7,"RF_R/mature_7_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_7 -feat ./RF_R/mature_7_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_7_imp","RF_R/mature_7_results.txt")))
new_list("RF_R/mature_7_imp",6,"RF_R/mature_6_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_6 -feat ./RF_R/mature_6_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_6_imp","RF_R/mature_6_results.txt")))
new_list("RF_R/mature_6_imp",5,"RF_R/mature_5_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_5 -feat ./RF_R/mature_5_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_5_imp","RF_R/mature_5_results.txt")))
new_list("RF_R/mature_5_imp",4,"RF_R/mature_4_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_4 -feat ./RF_R/mature_4_features.txt")

(feat_plot <- rbind(feat_plot,imp_plot("RF_R/mature_4_imp","RF_R/mature_4_results.txt")))
new_list("RF_R/mature_4_imp",3,"RF_R/mature_3_features")
system("python3.9 ./RF_python_scripts/ML_regression.py -df ./RF_R/RF_mature_tab.csv -alg RF -y_name area -gs T -cv 5 -n 100 -save ./RF_R/mature_3 -feat ./RF_R/mature_3_features.txt")

# #special export of the first importance plot
# # imp <- read.csv("mature_20180507_RF_chem_imp", sep = "\t", header=FALSE)
# # setEPS()
# # postscript("20180509_mature_imp.eps")
# # par(mar=c(5,5,4,2)) 
# # plot(V2~rownames(imp), data = imp, xlab = "Feature rank", ylab = "Feature importance score (Gini index)", main = "A) Mature-leaf palatability ~ 110 peaks", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
# # dev.off()
# # rm(imp)
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem2_imp","mature_20180507_RF_chem2_results.txt"))
# #remove 22 again (keep 66)
# new_list("mature_20180507_RF_chem2_imp",66,"mature_chem3")
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem3_imp","mature_20180507_RF_chem3_results.txt"))
# #remove 22 again (keep 44)
# new_list("mature_20180507_RF_chem3_imp",44,"mature_chem4")
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem4_imp","mature_20180507_RF_chem4_results.txt"))
# #remove 22 again (keep 22)
# new_list("mature_20180507_RF_chem4_imp",22,"mature_chem5")
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem5_imp","mature_20180507_RF_chem5_results.txt"))
# #keep half, 11. I tried 6 the first time and it did improve R^2 a bit, but it declined with 5 and 4 so I wondered whether I missed a peak.
# new_list("mature_20180507_RF_chem5_imp",11,"mature_chem6")
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem6_imp","mature_20180507_RF_chem6_results.txt"))
# #keep 9
# new_list("mature_20180507_RF_chem6_imp",9,"mature_chem7")
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem7_imp","mature_20180507_RF_chem7_results.txt"))
# #keep 8
# new_list("mature_20180507_RF_chem7_imp",8,"mature_chem8")
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem8_imp","mature_20180507_RF_chem8_results.txt"))
# #keep 7
# new_list("mature_20180507_RF_chem8_imp",7,"mature_chem9")
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem9_imp","mature_20180507_RF_chem9_results.txt"))
# #keep 6
# new_list("mature_20180507_RF_chem9_imp",6,"mature_chem10")
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem10_imp","mature_20180507_RF_chem10_results.txt"))
# #keep 5
# new_list("mature_20180507_RF_chem10_imp",5,"mature_chem11")
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem11_imp","mature_20180507_RF_chem11_results.txt"))
# #keep 4
# new_list("mature_20180507_RF_chem11_imp",4,"mature_chem12")
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem12_imp","mature_20180507_RF_chem12_results.txt"))
# #keep 3
# new_list("mature_20180507_RF_chem12_imp",3,"mature_chem13")
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem13_imp","mature_20180507_RF_chem13_results.txt"))
# #keep 2
# new_list("mature_20180507_RF_chem13_imp",2,"mature_chem14")
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem14_imp","mature_20180507_RF_chem14_results.txt"))
# #keep 2
# new_list("mature_20180507_RF_chem14_imp",1,"mature_chem15")
# 
# feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem15_imp","mature_20180507_RF_chem15_results.txt"))
# 
# feat_plot
# # plot(`R^2` ~ features, data = feat_plot, xlab = "Number of features", ylab = expression("Model "~R^2), main = "B) Feature selection", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
# # setEPS()
# # postscript("20180509_mature_chem_R.eps")
# # par(mar=c(5,5,4,2))
# # plot(`R^2` ~ features, data = feat_plot, xlab = "Number of features", ylab = expression("Model "~R^2), main = "B) Feature selection", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
# # dev.off()
# 
# #the winner is 3 compounds! wow! Model #13
# 
# #actual v predicted plot for final model
# # setEPS()
# # postscript("20180509_mature_chem_scores.eps")
# # par(mar=c(5,5,4,2))
# # plot(Y ~ Mean, data = read.csv("mature_20180507_RF_chem13_scores.txt", sep = "\t", header=TRUE), xlab = "Predicted values", ylab = "Actual values", main = "C) Fit of best model", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
# # dev.off()
# 
# rm(feat_plot)
# ```
# 
# 
# Feature lists for "pop + chem" model that uses both population and chemistry 
# ```{r}
# x1 <- read.csv("pop_list.txt", sep = "\t", header=FALSE)
# x2 <- read.csv("mature_chem13.txt", sep = "\t", header=FALSE)
# write.table(rbind(x1,x2), "mature_chem_pop.txt" , sep = "\t", row.names = FALSE, quote=FALSE, col.names=FALSE )
# 
# #python /mnt/home/azodichr/GitHub/ML-Pipeline/ML_regression.py -df mature_20180507 -alg RF -y_name conv -gs T -cv 5 -n 100 -tag chempop -feat mature_chem_pop.txt
# ```