Data Processing Workflow

Each raw dataset goes through a numbered Processing script, which outputs two numbered datasets: Processed_Indiv and Processed_Line. Processed_Indiv has individual measurements, or means for any traits measured multiple times per individual (toughness). Processed_Line has line means for traits measured on multiple individuals. LCMS data and ChemSummaries only have Processed_Indiv file because only one individual was sampled for LC/MS per line.

To go from raw to processed, each dataset should end up with the following columns. Some datasets only have line_age instead of pop, line, and age because line_age is the important one for merging, and pop, line, and age can be extracted from line_age

pos_age       Pasted pos_age (only in Processed_Indiv)

line_age      Pasted pop_line_age

pos or CupID           Position in the greenhouse; unique identifier WITHIN a trait dataset but do not assume pos 23 in LCMS is the same as pos 23 in palatability because these traits were measured on different cohorts of plants. For palatability, this is CupID instead, as each cup is a replicate trial of leaves pooled from a maternal line

pop           Population of origin

line          Maternal line nested in population of origin (ID of individual from which seeds were collected in the field or greenhouse)

age           Leaf age; young is expanding and not toughened up, mature is fully expanded and toughened up

lat           Latitude for the population of origin (only in Processed_Indiv)

region        Region for the population of origin (only in Processed_Indiv)

TRAITS        unique columns for each trait dataset


In the script combining all traits, the Processed_Line datasets are combined based on line_age and then lat and region are added