#################################################################################################################
### Regression coefficients with multiple phenotypes                                                          ###
### Author: Samuel Moix (adapted from Chiara Auwerx)                                                          ###
### Date: 18.04.2022                                                                                          ###
#################################################################################################################


# Extract phenotype information for binary diseases and generatte a PLINK-compatible file

#################################################
### Libraries ###################################

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#################################################
### Load Main Data ##############################

data_folder = "/SET/PATH/TO/DIRECTORY"
export_folder = "/SET/PATH/TO/DIRECTORY"

# Phenotype list; Contains the definition of inclusion and exclusion lists for each disease
pheno_list <- as.data.frame(fread(file.path(data_folder,  "ICD10_binary_traits_final.txt"), header = T))


#################################################
### STEP 1: Determine columns with phenotypes ###

header <- as.data.frame(fread("ukb44073.csv", header = T, nrow = 0))
col_SRC <- grep("^20001-", names(header))
col_SRD <- grep("^20002-", names(header))
col_ICD <- grep("^41270-", names(header))


#################################################
### STEP 2: Generate index tables ###############

# Self-reported cancer (#20001)
SRC <- as.data.frame(fread("ukb44073.csv", header = T, select = c(1, col_SRC), colClass = "character"))
SRC[SRC==""] <- NA
SRC <- na.omit(gather(SRC, Instance, Code, 2:ncol(SRC)))
SRC <- unique(SRC[, c(1,3)]) # 49'904 self-reported cancers

# Self-reported disease (#20002)
SRD <- as.data.frame(fread("ukb44073.csv", header = T, select = c(1, col_SRD), colClass = "character"))
SRD[SRD==""] <- NA
SRD <- na.omit(gather(SRD, Instance, Code, 2:ncol(SRD)))
SRD <- unique(SRD[, c(1,3)]) # 1'069'671 self-reported diseases

# ICD10 diagnosis (#41270)
ICD <- as.data.frame(fread("ukb44073.csv", header = T, select = c(1, col_ICD), colClass = "character"))
ICD[ICD==""] <- NA
ICD <- na.omit(gather(ICD, Instance, Code, 2:ncol(ICD)))
ICD <- unique(ICD[, c(1,3)]) #  4'131'218 ICD10 diagnosis

#################################################
### STEP 3: Generate phenotype tables ###########

# Load selected eids
eid_all <- as.data.frame(fread(file.path(data_folder,  "e_tasc_pheno.txt"), header = T, select = c("eid")))

# Generate empty table
pheno_all <- as.data.frame(matrix(data = 1, nrow = nrow(eid_all), ncol = nrow(pheno_list), dimnames = list(c(1:nrow(eid_all)), c(pheno_list$phenotype))))
pheno_all <- cbind(eid_all, pheno_all)

### 3.1 Exclusion list --> NA ###################
for (p in 1:nrow(pheno_list)) {
  # Define phenotype
  pheno <- pheno_list[p, "phenotype"]
  print(paste0("Starting to exclude individuals for ", pheno))

  # Excluded SRC (#20001)
  excl_SRC <- str_split(pheno_list[p, "exclude_20001"], ", ")[[1]]
  excl_SRC_eid <- unique(SRC[SRC$Code %in% excl_SRC, "eid"])

  # Excluded SRD (#20002)
  excl_SRD <- str_split(pheno_list[p, "exclude_20002"], ", ")[[1]]
  excl_SRD_eid <- unique(SRD[SRD$Code %in% excl_SRD, "eid"])

  # Excluded ICD (#41270)
  excl_ICD <- str_split(pheno_list[p, "exclude_41270"], ", ")[[1]]
  excl_ICD <- paste(paste0(sub("\\.", "", excl_ICD), ".*"), collapse = "|")
  excl_ICD_eid <- ICD %>% filter(str_detect(Code, excl_ICD))
  excl_ICD_eid <- unique(excl_ICD_eid$eid)

  # Merged exclusion list
  excl <- as.numeric(unique(c(excl_SRC_eid, excl_SRD_eid, excl_ICD_eid)))
  print(paste0("Number of individuals excluded for ", pheno, ": ", length(excl)))

  # Correct the table
  pheno_all[pheno_all$eid %in% excl, pheno] <- NA

}

rm()

### 3.2 Case list --> 2 #########################

for (p in 1:nrow(pheno_list)) {
  # Define phenotype
  pheno <- pheno_list[p, "phenotype"]
  print(paste0("Starting to include individuals for ", pheno))

  # Include ICD (#41270)
  incl_ICD <- str_split(pheno_list[p, "include_41270"], ", ")[[1]]
  incl_ICD <- paste(paste0(sub("\\.", "", incl_ICD), ".*"), collapse = "|")
  incl_ICD_eid <- ICD %>% filter(str_detect(Code, incl_ICD))
  incl_ICD_eid <- as.numeric(unique(incl_ICD_eid$eid))
  print(paste0("Number of individuals included for ", pheno, ": ", length(incl_ICD_eid)))

  # Correct the table
  pheno_all[pheno_all$eid %in% incl_ICD_eid, pheno] <- 2
}


#################################################
### STEP 4: Add a composite traits ##############

# Composite traits are calculated based on 55 diagnosis (exclude scoliosis, sex-specific diseases)
exc_traits <- c("eid", "scoliosis", "BC", "OC", "PC", "endometriosis", "menstruation")

# Disease suceptibility (having at least one diagnosis)
pheno_all[rowSums(pheno_all[, !names(pheno_all) %in% exc_traits] == 2, na.rm = T) == 0, "AnyDisease"] <- 1
pheno_all[rowSums(pheno_all[, !names(pheno_all) %in% exc_traits] == 2, na.rm = T) > 0, "AnyDisease"] <- 2

# Disease burden (Number of diagnosis)
pheno_all$DiseaseBurden <- rowSums(pheno_all[, !names(pheno_all) %in% c(exc_traits, "AnyDisease")] == 2, na.rm = T)

# All - Without removing female-specific traits
print(paste0("Dimensions of disease table for selected individuals: ", ncol(pheno_all)-1, " pheno x ", nrow(pheno_all), " individuals"))

#################################################
### Save Phenotypes #############################

# All
fwrite(pheno_all, file.path(data_folder,  "pheno_ICD10_All.txt"), col.names = T, row.names = F, quote = F, sep = "\t", na = "NA")

# ATTENTION WHEN DOING THE ANALYSIS sex-specific traits should be considered just for their sex (PC males other females)

