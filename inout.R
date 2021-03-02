### file 01_chem_wrangling.R

  #ln 43
  raw = read.csv("20180504_export_all.csv", header=TRUE)

  #ln 62
  feb_key <- read.csv("2018_feb_chem_key.csv", header=TRUE)
  #ln77
  dec_key <- read.csv("2017_dec_chem_key.csv", header=TRUE)
  
  #ln 92
  lat = read.csv("lats_long_names.csv", header = TRUE)
  
  #ln 149
  leaf <- read.csv("2017_chem_leaf_mass.csv", header=TRUE)
  
  #ln 177
  write.csv(abundant1, file = "20180531_chem_abundant1percent.csv")
  write.csv(abundant.5, file = "20180531_chem_abundant.5percent.csv")
  write.csv(abundant.1, file = "20180531_chem_abundant.1percent.csv")
  
### file 02_chem_NMDS
  
  #ln 32
  chem <- read.csv("20180507_chem_abundant1percent.csv", header = TRUE) #THIS IS output FROM 01_chem_wrangling
  
  #ln 41
  imp.y <- read.csv("20180517_imp_young.csv", header = TRUE)
  imp.m <- read.csv("20180517_imp_mature.csv", header = TRUE)
  
  #ln 137
  write.csv(scores, file = "20180531_nmds_scores.csv", row.names = FALSE )

### file 03_RF_model_select
  
  #python3.8 ML_regression.py -df mature_20180507 -alg RF -y_name conv -gs T -cv 5 -n 100 -tag chem -feat chem_list.txt
  ML_regression.py
  mature_20180507
  chem_list.txt
  
### file BB_all_traits_wrangling_RF_gri

 source("A1_cn_wrangling.R")
  # ln
  cn = read.csv("all_CN_20180422.csv", header=TRUE)
  
 source("A2_toughness_wrangling.R")
  raw = read.csv("20171005_toughness_age.csv", header=TRUE)
  palat = read.csv("line_level_palatability.csv", header = TRUE)

### --- --- --- 

  write.table(mature, "mature_20180507", sep="\t", row.names = FALSE)
  write.table(young, "young_20180507", sep="\t", row.names = FALSE)

### new dataset
  "LC-MS_redo all_jul 2018 -Sheet1.csv"
    # looks like
  "2018_feb_chem_key.csv"
    # but not cleaned up.
  
