##---
##title: "RF model selection"
##author: "Carina Baskett"
##date: "5/9/2018"
##output: html_document
##---
##model selection for RF analysis of re-aligned chemistry data (the best!)
##started May 8, 2018

##```{r setup, include=FALSE}
##knitr::opts_chunk$set(echo = FALSE)
library(ggplot2)
#setwd("/Volumes/baskettc/Documents/")
#will plot sorted importance scores and then output a vector with the number of features, model R^2, and SE of R^2
imp_plot <- function(imp_file,results_file)
{
  imp <- read.csv(imp_file, sep = "\t", header=FALSE)
  plot(V2~rownames(imp), data = imp)
  results <- read.csv(results_file, sep = "\t", header=FALSE)
  (r <- c(nrow(imp),as.numeric(as.character(results[24,1])),as.numeric(as.character(results[26,1]))))
}

#will make a new list of features to use in the next model, based on a cutoff of the importance file from last model
new_list <- function(imp_file,cutoff,tag)
{
  imp <- read.csv(imp_file, sep = "\t", header=FALSE)
  name <- as.character(paste(tag,".txt",sep=""))
  write.table(imp[1:cutoff,1], name, sep = "\t", row.names = FALSE, quote=FALSE, col.names=FALSE )
}

extract_r2 <- function(results_file)
{
  results <- read.csv(results_file, sep = "\t", header=FALSE)
  (as.numeric(as.character(results[24,1])))
}
extract_se <- function(results_file)
{
  results <- read.csv(results_file, sep = "\t", header=FALSE)
  (as.numeric(as.character(results[26,1])))
}

##```

##Model selection: mature leaf chemistry
##```{r}

feat_plot <- data.frame(matrix(NA, nrow = 1, ncol = 3))
names(feat_plot) <- c("features","R^2","SE")

#python /mnt/home/azodichr/GitHub/ML-Pipeline/ML_regression.py -df mature_20180507 -alg RF -y_name conv -gs T -cv 5 -n 100 -tag chem -feat chem_list.txt
#python3.8 ...ML_regression.py -df mature_20180507 -alg RF -y_name conv -gs T -cv 5 -n 100 -tag chem -feat chem_list.txt

chem_palat <- read.csv("Processing/2_out_AllTraits.csv",header=T)
write(colnames(chem_palat[24:ncol(chem_palat)]),"Analysis/chem_list.txt")
other_traits <- c("pop","age","lat","region","tough","percent_N","percent_C","C_N","log.abund","richness","diversity")
write(c(other_traits,colnames(chem_palat[24:ncol(chem_palat)])),"Analysis/all_traits_list.txt")

area_chem_only <- read.table("Processing/2_out_AllTraitsTab.csv",sep="\t")

write.table(drop_na(), "Processing/2_out_AllTraitsTab.csv",row.names=F,sep="\t")

system("python3.8 ./RF_python_scripts/ML_regression.py -df ./Processing/2_out_AllTraitsTab.csv -alg RF -y_name area -gs T -cv 2 -n 3 -tag chem -feat ./Analysis/chem_list_a.txt")
#system("python3.8 ./chri/ML_regression.py -df ./data_out/mature_20180531 -alg RF -y_name conv -gs T -cv 2 -n 3 -tag chem -feat ./data_out/chem_list.txt")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem_imp","mature_20180507_RF_chem_results.txt"))
#remove 20% (22 compounds; so keep 88) for first cut
new_list("mature_20180507_RF_chem_imp",20,"mature_chem2") #88
system("python3.8 ./chri/ML_regression.py -df ./chri/mature_20180507 -alg RF -y_name conv -gs T -cv 2 -n 3 -tag chem -feat ./mature_chem2.txt", ignore.stdout = FALSE)
feat_plot <- na.omit(feat_plot)

#special export of the first importance plot
# imp <- read.csv("mature_20180507_RF_chem_imp", sep = "\t", header=FALSE)
# setEPS()
# postscript("20180509_mature_imp.eps")
# par(mar=c(5,5,4,2)) 
# plot(V2~rownames(imp), data = imp, xlab = "Feature rank", ylab = "Feature importance score (Gini index)", main = "A) Mature-leaf palatability ~ 110 peaks", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
# dev.off()
# rm(imp)

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem2_imp","mature_20180507_RF_chem2_results.txt"))
#remove 22 again (keep 66)
new_list("mature_20180507_RF_chem2_imp",66,"mature_chem3")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem3_imp","mature_20180507_RF_chem3_results.txt"))
#remove 22 again (keep 44)
new_list("mature_20180507_RF_chem3_imp",44,"mature_chem4")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem4_imp","mature_20180507_RF_chem4_results.txt"))
#remove 22 again (keep 22)
new_list("mature_20180507_RF_chem4_imp",22,"mature_chem5")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem5_imp","mature_20180507_RF_chem5_results.txt"))
#keep half, 11. I tried 6 the first time and it did improve R^2 a bit, but it declined with 5 and 4 so I wondered whether I missed a peak.
new_list("mature_20180507_RF_chem5_imp",11,"mature_chem6")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem6_imp","mature_20180507_RF_chem6_results.txt"))
#keep 9
new_list("mature_20180507_RF_chem6_imp",9,"mature_chem7")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem7_imp","mature_20180507_RF_chem7_results.txt"))
#keep 8
new_list("mature_20180507_RF_chem7_imp",8,"mature_chem8")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem8_imp","mature_20180507_RF_chem8_results.txt"))
#keep 7
new_list("mature_20180507_RF_chem8_imp",7,"mature_chem9")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem9_imp","mature_20180507_RF_chem9_results.txt"))
#keep 6
new_list("mature_20180507_RF_chem9_imp",6,"mature_chem10")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem10_imp","mature_20180507_RF_chem10_results.txt"))
#keep 5
new_list("mature_20180507_RF_chem10_imp",5,"mature_chem11")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem11_imp","mature_20180507_RF_chem11_results.txt"))
#keep 4
new_list("mature_20180507_RF_chem11_imp",4,"mature_chem12")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem12_imp","mature_20180507_RF_chem12_results.txt"))
#keep 3
new_list("mature_20180507_RF_chem12_imp",3,"mature_chem13")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem13_imp","mature_20180507_RF_chem13_results.txt"))
#keep 2
new_list("mature_20180507_RF_chem13_imp",2,"mature_chem14")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem14_imp","mature_20180507_RF_chem14_results.txt"))
#keep 2
new_list("mature_20180507_RF_chem14_imp",1,"mature_chem15")

feat_plot <- rbind(feat_plot,imp_plot("mature_20180507_RF_chem15_imp","mature_20180507_RF_chem15_results.txt"))

feat_plot
# plot(`R^2` ~ features, data = feat_plot, xlab = "Number of features", ylab = expression("Model "~R^2), main = "B) Feature selection", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
# setEPS()
# postscript("20180509_mature_chem_R.eps")
# par(mar=c(5,5,4,2))
# plot(`R^2` ~ features, data = feat_plot, xlab = "Number of features", ylab = expression("Model "~R^2), main = "B) Feature selection", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
# dev.off()

#the winner is 3 compounds! wow! Model #13

#actual v predicted plot for final model
# setEPS()
# postscript("20180509_mature_chem_scores.eps")
# par(mar=c(5,5,4,2))
# plot(Y ~ Mean, data = read.csv("mature_20180507_RF_chem13_scores.txt", sep = "\t", header=TRUE), xlab = "Predicted values", ylab = "Actual values", main = "C) Fit of best model", cex.lab=1.75, cex.axis=1.75, cex.main=1.75, cex.sub=1.75)
# dev.off()

rm(feat_plot)
```


Feature lists for "pop + chem" model that uses both population and chemistry 
```{r}
x1 <- read.csv("pop_list.txt", sep = "\t", header=FALSE)
x2 <- read.csv("mature_chem13.txt", sep = "\t", header=FALSE)
write.table(rbind(x1,x2), "mature_chem_pop.txt" , sep = "\t", row.names = FALSE, quote=FALSE, col.names=FALSE )

#python /mnt/home/azodichr/GitHub/ML-Pipeline/ML_regression.py -df mature_20180507 -alg RF -y_name conv -gs T -cv 5 -n 100 -tag chempop -feat mature_chem_pop.txt
```