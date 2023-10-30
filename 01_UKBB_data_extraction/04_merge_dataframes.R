#################################################################################################################
### Merge dataframes                                                                                          ###
### Author: Samuel Moix                                                                                       ###
### Date: 19.10.2022                                                                                          ###
#################################################################################################################

library(data.table)
library(dplyr)

data_folder = "/SET/PATH/TO/DIRECTORY"
export_folder = "/SET/PATH/TO/DIRECTORY"

traits <- as.data.frame(fread(file.path(data_folder,  "phenotype_wo_diseases.txt"), header = T)) 


disease <- as.data.frame(fread(file.path(data_folder,  "pheno_ICD10_All.txt"), header = T))

pcs_array <- as.data.frame(fread(file.path(data_folder,  "df_lm_samples.txt"), header = T, drop = c("age", "sex", "blood_cancer_TF")))
pcs_array <- rename(pcs_array, TL_check_array = Zadj_TS)

temp <- merge(traits, disease, by="eid", all.x = TRUE)
final_df <- merge(temp, pcs_array, by="eid", all.x = TRUE)

print(paste0("Dimensions of phenotype table: ", ncol(final_df), " pheno x ", nrow(final_df), " eids"))

fwrite(final_df, file.path(export_folder,  "complete_DF.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

