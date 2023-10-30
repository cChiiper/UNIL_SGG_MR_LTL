#################################################################################################################
### Retrieve samples used for the linear regression analysis                                                  ###
### Author: Samuel Moix                                                                                       ###
### Date: 07.03.2022                                                                                          ###
### Code adapted from Chiara Auwerx: Retrieve samples used for the CNV-GWAS analysis                          ###
### https://github.com/cauwerx/CNV_GWAS_continuous_traits                                                     ###
#################################################################################################################

################################################
### Libraries ##################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

################################################################################################################
### Load Main Data #############################

# Sample QC; This corresponds to the "ukb_sqc_v2.txt" file available from the UKBB portal
sample_qc <- as.data.frame(fread("ukb_sqc_v2.txt", header = F, select = c(3:68)))
colnames(sample_qc) <- c("array", "batch", "plate", "well", "cluster_CR", "dQC", "dna_concentration", "submitted_gender", "inferred_gender", "X_intensity", "Y_intensity", "submitted_plate", "submitted_well", "missing_rate", "heterozygosity", "heterozygosity_pc_corrected", "heterozygosity_missing_outlier", "PSCA", "in_kinship", "excluded_kinship_inference", "excess_relatives", "white_british", "pca_calculation", paste0("PC", seq(1,40)), "phasing_autosome", "phasing_X", "phasing_Y")

# Sample eid; This corresponds to a list of all sample eids (with sex information), retrieved from a ".fam" file available from the UKBB portal
sample_eid <- as.data.frame(fread("_001_ukb_cal_chr1_v2.fam", header = F, select = c(1,5), col.names = c("eid", "sex")))

#################################################
### Merge Main Data #############################

df <- cbind(sample_eid, sample_qc)
print(paste0("Start: ", nrow(df), " individuals"))

#################################################
### Add phenotypes of interest ##################


### Add age ####################################
# Identify column containing age (field 21003) from "ukb28603.csv", a phenotype file available from the UKBB portal
header <- fread("pheno/ukb28603.csv", header = T, nrow = 0)
col_age <- grep("21003-", names(header))

# Extract age, under select 1 to take eid and col_age[1] for the first instance 
age <- as.data.frame(fread("ukb28603.csv", header = T, select = c(1, col_age[1]), col.names = c("eid", "age")))

# Add to main dataframe
df <- right_join(age, df, by = "eid")

### Add Telomere length ########################
header <- fread("ukb46831.csv", header = T, nrow = 0)
col_TS <- grep("22192-", names(header))

# Extract telomere length as Z adj. T/S
Zadj_TS <- as.data.frame(fread("ukb46831.csv", header = T, select = c(1, col_TS[1]), col.names = c("eid", "Zadj_TS")))

# Add to main dataframe
df <- right_join(Zadj_TS, df, by = "eid")

################################################################################################################
### STEP 1: Exclude related samples (pca_calculation = 1)

df <- df[which(df$pca_calculation == 1), ]
print(paste0("STEP 1: Exclude related samples: ", nrow(df), " individuals"))

#################################################
### STEP 2: Exclude non-white, non-British ancestry samples (white_british = 1)

df <- df[which(df$white_british == 1), ]  
print(paste0("STEP 2: Exclude non-white, non-British ancestry samples: ", nrow(df), " individuals"))


#################################################
### STEP 4: Exclude retracted samples; This file was available from the UKBB portal 
retracted <- as.data.frame(fread("w16389_20220222.csv", header = F, col.names = "eid"))
df <- df[!df$eid %in% retracted$eid, ]
print(paste0("STEP 4: Exclude retracted samples: ", nrow(df), " individuals"))

#################################################
### STEP 5: Exclude samples with non-matching submitted vs. inferred sex
df <- df[which(df$submitted_gender == df$inferred_gender), ]
print(paste0("STEP 5: Exclude sex mismatches: ", nrow(df), " individuals"))

#################################################
## STEP 6: Find samples with self-reported blood malignancies or ICD10 blood malignancy diagnosis
print("STEP 6: Mark samples with blood malignancies")

# UKBB field 20001 (self-reported cancers), codes: 1047 (lymphoma), 1048 (leukemia), 1050 (multiple myeloma), 1051 (myelofibrosis), 1052 (Hodgkin's lymphoma), 1053 (Non-Hodgkin's lymphoma), 1055 (chronic lymphocytic), 1056 (chronic myeloid), 1058 (other hematological malignancies)
SRBM_code <- "1047|1048|1050|1051|1052|1053|1055|1056|1058"


# All ICD10 diagnoses mapping to the PheCode "cancer of lymphatic and hematopoietic tissue"; "phecode_ICD10_091220.csv" was downoladed from https://phewascatalog.org/phecodes_icd10 on 11/03/2022 
# phecode <- as.data.frame(fread("phecode_icd10.csv ", header = T, select = c(1,6)))

phecode = read.csv("phecode_icd10.csv")

phecode = phecode[,c(1,6)]

DBM_code <- phecode[which(phecode$Excl_Phenotypes == "cancer of lymphatic and hematopoietic tissue"), "ICD10"]


print(paste0("Number of cancer of lymphatic and hematopoietic tissue: ", length(DBM_code)))


DBM_code<- paste(sub("\\.", "", DBM_code[grep("\\.", DBM_code)]), collapse = "|")


# Select appropriated columns; "ukb44073.csv" contains phenotypic data for field 20001 (self-reported cancer) and 41270 (ICD10 diagnoses) available from the UKBB portal
header_pheno <- fread("ukb44073.csv", nrow = 0)
col_SRBM <- grep("20001-", names(header_pheno))
col_DBM <- grep("41270-", names(header_pheno))


# Self-reported blood malignancy
SRBM <- as.data.frame(fread("ukb44073.csv", select = c(1, col_SRBM), header = T, colClass = "character"))
SRBM[SRBM == ""] <- NA
SRBM <- unite(SRBM, code, names(SRBM[2:ncol(SRBM)]), sep = "_", remove = T, na.rm = T)
SRBM <- SRBM %>% filter(str_detect(code, SRBM_code))
print(paste0("Number of samples with self-reported blood malignancy: ", nrow(SRBM)))


# ICD10 diagnosis
DBM <- as.data.frame(fread("ukb44073.csv", select = c(1, col_DBM), header = T))
DBM[DBM == ""] <- NA
DBM <- unite(DBM, code, names(DBM[2:ncol(DBM)]), sep = "_", remove = T, na.rm = T)
DBM <- DBM %>% filter(str_detect(code, DBM_code))
print(paste0("Number of samples ICD10 diagnosed blood malignancy: ", nrow(DBM)))

# Total number of blood malignancies
BM <- union(SRBM$eid, DBM$eid)
print(paste0("Total number of samples with blood malignancy: ", length(BM)))

# Add column with True for patient with blood malignancy 
df <- df %>% mutate(blood_cancer_TF=ifelse(eid %in% BM,T,F))


################################################################################################################
### Save ########################################
fwrite(df, "df_lm_samples.txt", col.names = T, row.names = F, quote = F, sep = "\t")

# For use with PLINK v2.0 - All
fwrite(df[, "eid", drop = F], "lm_samples_All.txt", col.names = F, row.names = F, quote = F, sep = "\t")
