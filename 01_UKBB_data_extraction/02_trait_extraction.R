################################################################################
### Extract data from the UKB                                                ###
### Author: Samuel Moix                                                      ###
### Date: 19.10.2022                                                         ###
### Code adapted from Chiara Auwerx                                          ###
### https://github.com/cauwerx/CNV_GWAS_continuous_traits                    ###
################################################################################

################################################
### Libraries ##################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)

################################################
### Working directories ########################
data_folder <- "/SET/PATH/TO/DIRECTORY"
export_folder <- "/SET/PATH/TO/DIRECTORY"

################################################
### STEP 1: Load metadata ######################

### Load metadata file and filter rows with FieldIDs
metadata <- as.data.frame(fread(file.path(data_folder,  "TL_metadata.csv"), header = T))

### Adjust format 
pheno_list <- metadata %>% 
  select(FieldID, pheno, type) %>%
  filter(!is.na(FieldID)) %>%
  rename(PHENO = pheno) 

print("Metadata loaded")

#################################################
### STEP 2: Order raw phenotype files ###########

# As phenotypic data is spread throughout multiple files downloaded from the UKBB portal, we start by listing all phenotype files available and order them by date
raw_pheno_file <- file.info(list.files("/SET/PATH/TO/DIRECTORY", pattern = "^ukb.*csv", full.names = T, recursive = F))
raw_pheno_file <- data.frame(File = rownames(raw_pheno_file), Date = raw_pheno_file[,4])
raw_pheno_file <- raw_pheno_file[order(raw_pheno_file$Date, decreasing = T), ]
print(paste0("There are ", nrow(raw_pheno_file), " raw phenotype files"))

#################################################
### STEP 3: Identify most recent location #######
# As mentioned in STEP 2, there are multiple phenotype files, with some phenotypes being in several files
# For the purpose of the analysis, we use phenotype information originating from the most recent phenotype file

pheno_list$File <- NA
pheno_list$File <- as.character(pheno_list$File) 
pheno_list$Col <- NA
pheno_list$Col <- as.character(pheno_list$Col) 

# Loop over all phenotypes to detect the most recent location
for (p in 1:nrow(pheno_list)) {
  
  # Define the phenotype
  pheno <- pheno_list[p, "PHENO"] 
  ID <- pheno_list[p, "FieldID"] 
  print(paste0("Extracting most recent file for ", pheno))
  counter <- 1
  
  # Determine the most recent location --> loop through raw files
  while (is.na(pheno_list[p, "File"])) {
    
    header <- fread(as.character(raw_pheno_file[counter, "File"]), nrow = 0)
    col <- grep(paste0("^", ID, "-"), names(header))
    if (length(col) > 0) {
      pheno_list[p, "File"] <- as.character(raw_pheno_file[counter, "File"])
      print(paste0("Most recent location: ", sub(".*/", "", pheno_list[p, "File"])))
      pheno_list[p, "Col"] <- paste(col, collapse = "_")}
    counter <- counter +1
  }
} 
rm(p, pheno, ID, counter, header, col)
print("Writing phenotype file list to pheno_list.txt")
fwrite(pheno_list, file.path(export_folder,  "pheno_list.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

### Save full pheno_list
full_pheno_list <- pheno_list

#################################################
### STEP 4: Extract phenotypes  #################

### Get base File generated with 01_sample_filtering and extract main phenotypes
phenotypes <- as.data.frame(fread(file.path(data_folder, "df_lm_samples.txt") , header = T))
phenotypes <- select(phenotypes, c("eid","Zadj_TS","age","sex","blood_cancer_TF"))
phenotypes <- drop_na(phenotypes, c("eid","Zadj_TS","age","sex"))


fwrite(phenotypes, file.path(export_folder,  "e_tasc_pheno.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
print("Base file loaded (df_lm_samples.txt) and saved")

################################
### Select continuous traits ###
pheno_list <- filter(full_pheno_list, type == "continuous")

# Loop over all phenotypes to extract the data and calculate the average over different measurement instances
for(p in 1:nrow(pheno_list)) {
  
  # Define the phenotype
  pheno <- pheno_list[p, "PHENO"] 
  ID <- pheno_list[p, "FieldID"] 
  file <- pheno_list[p, "File"] 
  col <- pheno_list[p, "Col"] 
  print(paste0("Extracting ", pheno, " from ", sub(".*/", "", file)))
  
  # Extract columns corresponding to the defined phenotype
  temp_pheno <- as.data.frame(fread(file, header = T, select = c(1, as.numeric(str_split(col, pattern = "_")[[1]]))))
  
  # Put NA's in negative values
  traits_with_negatives <- c("TDI","TL","Zadj_TS") 
  if(!pheno %in% traits_with_negatives){
    for(i in 2:length(temp_pheno)){
      temp_pheno[,i][which(temp_pheno[,i] < 0)] <- NA
    }
  }
  
  
  # Average (if more than one column present)
  if (ncol(temp_pheno) > 2) {
    temp_pheno <- data.frame(eid = temp_pheno[, 1], pheno = rowMeans(temp_pheno[, 2:ncol(temp_pheno)], na.rm = T))}
  
  # Merge
  colnames(temp_pheno) <- c("eid", pheno)
  phenotypes <- merge(phenotypes, temp_pheno, by="eid", all.x = TRUE)
  
}
rm(pheno, ID, file, col, temp_pheno)
print(paste0("Dimensions of continuous phenotype table: ", ncol(phenotypes), " pheno x ", nrow(phenotypes), " eids"))

################################
### Select integer traits ######
pheno_list <- filter(full_pheno_list, type == "integer")

# Loop over all phenotypes to extract the data
for(p in 1:nrow(pheno_list)) {
  
  # Define the phenotype
  pheno <- pheno_list[p, "PHENO"] 
  ID <- pheno_list[p, "FieldID"] 
  file <- pheno_list[p, "File"] 
  col <- pheno_list[p, "Col"] 
  print(paste0("Extracting ", pheno, " from ", sub(".*/", "", file)))
  
  # Extract columns corresponding to the defined phenotype
  temp_pheno <- as.data.frame(fread(file, header = T, select = c(1, as.numeric(str_split(col, pattern = "_")[[1]]))))
  
  # Select (if more than one column present)
  if (ncol(temp_pheno) > 2) {
    temp_pheno <- data.frame(eid = temp_pheno[, 1], pheno = temp_pheno[, 2])}
  
  # Merge
  colnames(temp_pheno) <- c("eid", pheno)
  phenotypes <- merge(phenotypes, temp_pheno, by="eid", all.x = TRUE)
  
}
rm(pheno, ID, file, col, temp_pheno)
print(paste0("Dimensions of phenotype table (cont, int): ", ncol(phenotypes), " pheno x ", nrow(phenotypes), " eids"))

################################
### Select factors #############
pheno_list <- filter(full_pheno_list, type == "factor")

# Loop over all phenotypes to extract the data
for(p in 1:nrow(pheno_list)) {
  
  # Define the phenotype
  pheno <- pheno_list[p, "PHENO"] 
  ID <- pheno_list[p, "FieldID"] 
  file <- pheno_list[p, "File"] 
  col <- pheno_list[p, "Col"] 
  print(paste0("Extracting ", pheno, " from ", sub(".*/", "", file)))
  
  # Extract columns corresponding to the defined phenotype
  temp_pheno <- as.data.frame(fread(file, header = T, select = c(1, as.numeric(str_split(col, pattern = "_")[[1]]))))
  
  # Select (if more than one column present)
  if (ncol(temp_pheno) > 2) {
    temp_pheno <- data.frame(eid = temp_pheno[, 1], pheno = temp_pheno[, 2])}
  
  # Merge
  colnames(temp_pheno) <- c("eid", pheno)
  phenotypes <- merge(phenotypes, temp_pheno, by="eid", all.x = TRUE)
  
}
rm(pheno, ID, file, col, temp_pheno)
print(paste0("Dimensions of final phenotype table: ", ncol(phenotypes), " pheno x ", nrow(phenotypes), " eids"))

#################################################
### STEP 5: Export file  ########################
fwrite(phenotypes, file.path(export_folder,  "phenotype_wo_diseases.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

