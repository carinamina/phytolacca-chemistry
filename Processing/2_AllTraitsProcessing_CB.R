#output will be the following traits merged by line_age, with information on latitude and region:
# palatability (script 1e)
# chemical abundance, richness, diversity (1b)
# toughness (1c)
# carbon, nitrogen, C:N (1d)
# ~2000 chemical abundance peaks (1a)

palat <- read.csv("Processing/1e_out_Palat_Line.csv", header=T)
chemsums <- read.csv("Processing/1b_out_ChemSummaries_Indiv.csv",header=T)
tough <- read.csv("Processing/1c_out_Toughness_Line.csv",header=T)
cn <- read.csv("Processing/1d_out_CarbonNitrogen_Line.csv",header=T)
compounds <- read.csv("Processing/1a_out_LCMS_Indiv.csv",header=T)

#chemsums and compounds have latitude etc, the others all just line_age

all.traits <- full_join(tough,cn,by="line_age") %>% full_join(palat,by="line_age") %>% full_join(chemsums,by="line_age")

#why does chem have bulk but not toughness? probably we can't afford to keep that label bc we lose too much data


rm(list=ls())
