################################################################################
### Regression coefficients with multiple phenotypes                         ###
### Author: Samuel Moix                                                      ###
### Date: 04.04.2022                                                         ###
################################################################################

################################################
### Libraries ##################################
library(data.table) # Read data 
library(ggplot2) # Plotting basics
library(ggExtra) # Extra plotting XP
library(ggpubr) # Plotting regression stuff
library(dplyr)
library(tibble)
library(forcats)
require(pheatmap)
library(ggforestplot) # Forest plot
library(conflicted)

################################################
### Working directories ########################
data_folder = "/SET/PATH/TO/DIRECTORY/"
export_folder = "/SET/PATH/TO/DIRECTORY/"

################################################
### Parameters  ################################

### Scale all data
par_scale <- TRUE
###  Use only continuous traits
par_only_cont <- FALSE
###  Use only factorial traits
par_only_fact <- FALSE
###  Calculate number of independent traits automatically if FALSE set (indep_trait <- XXX)
par_calc_indep_traits <- FALSE
indep_trait <- 141
### Whether to remove outliers from data
par_outlier_filter <- TRUE
### Whether to remove patients with blood cancer
par_rm_blood_cancer <- TRUE
### Correcting parameters
### Correct traits (cholesterol for medication and female traits for SES)
par_corr_traits <- FALSE
par_keep_not_cor <- FALSE # Whether to keep traits not corrected 

### Export file namer:
par_exp_name <- function(par_export_name){
  if(par_only_cont){par_export_name <- paste0(par_export_name, "_cont")}
  if(par_only_fact){par_export_name <- paste0(par_export_name, "_fact")}
  if(par_scale){par_export_name <- paste0(par_export_name, "_scaled")}
  if(par_corr_traits){par_export_name <- paste0(par_export_name, "_cor_traits")}
  if(par_keep_not_cor){par_export_name <- paste0(par_export_name, "_and_NC")}
  if(par_outlier_filter){par_export_name <- paste0(par_export_name, "_no_outliers")}
  if(par_rm_blood_cancer){par_export_name <- paste0(par_export_name, "_no_bc")}
  return(par_export_name)
}

################################################
### Load Data ##################################

df <- as.data.frame(fread(file.path(data_folder,  "filtered_DF.txt"), header = T)) 
metadata <- read.csv(file = "TL_metadata.csv", sep = ',', header = TRUE)

### Remove blood cancer patients 
if(par_rm_blood_cancer){
  # Remove patients with blood cancer
  df <- df[which(df$blood_cancer == 0 | is.na(df$blood_cancer)),]
  df <- select(df, -c("blood_cancer"))
}

### Remove missing values of blood traits

# define the columns to check for NAs
cols_to_check <- c("age", "sex", "age2", "array", "eosinophil_count", 
                   "lymphocyte_count", "monocyte_count", "neutrophil_count", 
                   "platelet_count", "RBC_count", "reticulocyte_count", "WBC_count")

# remove rows with NAs in the selected columns
df <- df[complete.cases(df[, cols_to_check]),]
rm(cols_to_check)

### Linear regression model TL corrected ALSO FOR BLOOD TRAITS
df$TLc <- residuals(lm(df$Zadj_TS ~ df$age + df$sex + df$age:df$sex + df$age2 + 
                         df$age2:df$sex + df$array + df$array:df$sex + 
                         df$eosinophil_count + df$lymphocyte_count + 
                         df$monocyte_count + df$neutrophil_count + 
                         df$platelet_count + df$RBC_count + 
                         df$reticulocyte_count + df$WBC_count))

################################################
### Correct Data ###############################

### Correct traits
if(par_corr_traits){
  # Correct cholesterol for cholesterol medication
  # Changes in mmol/L from https://pubmed.ncbi.nlm.nih.gov/22354153/
  # and from Simvastatin https://bmcprimcare.biomedcentral.com/articles/10.1186/1471-2296-4-18 
 
  df <- df %>% add_column(cholesterol_cor = ifelse(df$cholest_medication == 0, df$cholesterol,
                                                   df$cholesterol + 1.64) # Patients without CV
                    , .after = "cholesterol" )
  
  df <- df %>% add_column(HDL_cor = ifelse(df$cholest_medication == 0, df$HDL,
                                                   df$HDL - 0.1) # All doses
                          , .after = "HDL" )
  
  df <- df %>% add_column(LDL_cor = ifelse(df$cholest_medication == 0, df$LDL,
                                                   df$LDL + 1.4) # All doses
                          , .after = "LDL" )
  
  df <- df %>% add_column(TG_cor = ifelse(df$cholest_medication == 0, df$TG,
                                           df$TG + 0.4) # All doses
                          , .after = "TG" )
  
  if(!par_keep_not_cor){
    df$cholesterol <- df$cholesterol_cor
    df$HDL <- df$HDL_cor
    df$LDL <- df$LDL_cor
    df$TG <- df$TG_cor
    df <- select(df, -c("cholesterol_cor","HDL_cor","LDL_cor","TG_cor"))
  }
  
  # Correct female trait for SES
  trait_to_correct_SES <- c("age_first_birth", "age_last_birth",
                            "menopause", "reproductive_lp", "menstruation",
                            "birth_weight_first_child","nb_birth")
  
  if(par_keep_not_cor){
    df <- add_column(df, "age_first_birth_cor_SES" = NA, .after = "age_first_birth") 
    df <- add_column(df, "age_last_birth_cor_SES" = NA, .after = "age_last_birth") 
    df <- add_column(df, "menopause_cor_SES" = NA, .after = "menopause") 
    df <- add_column(df, "reproductive_lp_cor_SES" = NA, .after = "reproductive_lp") 
    df <- add_column(df, "menstruation_cor_SES" = NA, .after = "menstruation") 
    df <- add_column(df, "birth_weight_first_child_cor_SES" = NA, .after = "birth_weight_first_child") 
    df <- add_column(df, "nb_birth_cor_SES" = NA, .after = "nb_birth") 
  }
  
  temp_name <- ""
  if(par_keep_not_cor){
    temp_name <- "_cor_SES"
  }
  
  for (trait in trait_to_correct_SES) { # Correct for all traits in list
    df[paste0(trait, temp_name)] <- residuals(lm(df[,which(colnames(df) == trait)] ~ df$TDI + df$household_income +
                                                   df$age_end_education, na.action=na.exclude))
  }
  rm(trait)
  rm(temp_name)
  
  if(FALSE){
    # Correct traits for BMI
    trait_to_correct_BMI <- c("hip_cir") #Not WHR, WHRadjBMI
    
    if(par_keep_not_cor){
      df <- add_column(df, "hip_cir_cor_BMI" = NA, .after = "hip_cir") 
    }
    
    temp_name <- ""
    if(par_keep_not_cor){
      temp_name <- "_cor_BMI"
    }
    
    for (trait in trait_to_correct_BMI) { # Correct for all traits in list
      # Subset females for which both traits don't have any NA's
      df_F <- subset(df, sex == 2 & !is.na(df[,which(colnames(df) == trait)]) & !is.na(df$BMI))
      # Regression
      res_F <- residuals(lm(df_F[,which(colnames(df) == trait)] ~ df_F$BMI, na.action=na.exclude))
      # Add residuals to the main data frame
      df[which(df$sex == 2 & !is.na(df[,which(colnames(df) == trait)]) & !is.na(df$BMI)), paste0(trait, temp_name)] <- res_F
      rm(df_F);rm(res_F)
      # Repeat with males
      df_M <- subset(df, sex == 1 & !is.na(df[,which(colnames(df) == trait)]) & !is.na(df$BMI))
      res_M <- residuals(lm(df_M[,which(colnames(df) == trait)] ~ df_M$BMI, na.action=na.exclude))
      df[which(df$sex == 1 & !is.na(df[,which(colnames(df) == trait)]) & !is.na(df$BMI)), paste0(trait, temp_name)] <- res_M
      rm(df_M);rm(res_M)
      # Put NA's on trait's that contained NA's (needed because of the loop syntaxe)
      df[is.na(df[,which(colnames(df) == trait)]) | is.na(df$BMI),which(colnames(df) == paste0(trait, temp_name))] <- NA
    }
    rm(trait)
    rm(temp_name)
  }
  
}


################################################
### Segregation by type ########################

### Create df only containing continuous traits ###
continuous_cols <- metadata$type %in% c("continuous", "integer")
continuous_cols_names <- names(df)[match(metadata$pheno[continuous_cols], names(df))]
continuous_cols_names <- c(continuous_cols_names,"TLc")
contvec <- sort(match(continuous_cols_names, names(df)))
rm(continuous_cols, continuous_cols_names)

if (par_only_cont){
  df_continuous <- df[,contvec]

  if (par_scale){
    ### Scale continuous data and replace in data frame ###
    df_continuous <- as.data.frame(scale(df_continuous)) ###
  }
  
  #### Change input df
  df <- df_continuous
  rm(df_continuous)
}

if (par_only_fact){
  ### Create df only containing factorial traits
  factor_cols <- metadata$type == "factor"
  factor_cols_names <- names(df)[match(metadata$pheno[factor_cols], names(df))]
  factor_cols_names <- c(factor_cols_names,"TLc")
  factvec <- sort(match(factor_cols_names, names(df)))
  rm(factor_cols, factor_cols_names)
  df <- df[,factvec]
}

################################################
### Calculating number of independent traits ###

if (par_calc_indep_traits){
  ### Create a correlation matrix, do not consider NA
  pheno_cor <- cor(df, use="pairwise.complete.obs")
  ### Set missing values to 0 to allow svd() ###
  pheno_cor[which(is.na(pheno_cor))] <- 0
  ### Calculate the eigenvalues ###
  # svd() computes the singular-value decomposition of a rectangular matrix, with $d being a vector containing the singular values of the decomposition
  pheno_EV <- svd(pheno_cor)$d
  ### Calculate Neff ###
  # Neff is defined as the number of eigenvalues required to explain 99.5% of the variation from the CNV data 
  sum_EV <- 0
  count <- 0
  while(sum_EV/sum(pheno_EV) < 0.995) {
    count <- count + 1
    sum_EV <- sum_EV + pheno_EV[count]
  }
  ### Set number of independent traits
  indep_trait <- count
}

################################################
### Remove continuous traits outliers from data #

# Outliers defined as +/-5SDs
if(par_outlier_filter & !par_only_fact){
  print("Outliers filtered")
  for (i in contvec) {
    upper <- mean(df[,i], na.rm = T) + sd(df[,i], na.rm = T)*5
    lower <- mean(df[,i], na.rm = T) - sd(df[,i], na.rm = T)*5
    df[,i][which(df[,i] > upper | df[,i] < lower)] <- NA
    rm(upper)
    rm(lower)
  }
}
rm(contvec)

################################################
### Total scaling ##############################

if (par_scale){
  #### To put all traits as numeric
  df[,1:length(df)] <- lapply(df[,1:length(df)], as.numeric) ###
  ### Scale whole data
  df <- as.data.frame(scale(df)) ###
}

################################################
### Calculating regression coef ################

### Function to calculate p-value
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


### Check if it was already run
# Set file name
lm_res_file_name <- "blood_cor_lm_res"
lm_res_file_name <- par_exp_name(lm_res_file_name)
lm_res_file_name <- paste0(lm_res_file_name,".tsv")
print(paste("File name:", lm_res_file_name))

#Check weather file exists
res_file_path <- file.path(export_folder, lm_res_file_name)
if (file.exists(res_file_path)){
  print("Results file exists and loaded")
  df_res <- as.data.frame(fread(res_file_path, header = T))
}else{
  ### Regression coefficients on df
  print("Results file does not exist and is being generated")
  x <- c("pheno", "mean", "sd", "rcoef", "CI1", "CI2", "p_value", "se")
  
  df_res <- data.frame(matrix(ncol = length(x), nrow = ncol(df)))
  colnames(df_res) <- x
  
  for(i in 1:ncol(df)) {       # for-loop over columns
    reg_TLc_pheno <- lm(df$TLc ~ df[, i])
    df_res[i,] <- c(colnames(df)[i], 
                    if(class(df[,i]) == 'integer' | class(df[,i]) == 'numeric'){
                      mean(df[, i], na.rm = T)
                    } else NA,
                    if(class(df[,i]) == 'integer' | class(df[,i]) == 'numeric'){
                      sd(df[, i], na.rm = T)
                    } else NA,
                    reg_TLc_pheno$coefficients[[2]],
                    confint(reg_TLc_pheno)[2,1],
                    confint(reg_TLc_pheno)[2,2],
                    lmp(reg_TLc_pheno),
                    coef(summary(reg_TLc_pheno))[2, 2])
    print(colnames(df)[i])
    print(reg_TLc_pheno)
  }
  
  
  
  # Put data in df_res as numeric and factor
  df_res$pheno <- as.factor(df_res$pheno)
  df_res[,2:ncol(df_res)] <- lapply(df_res[,2:ncol(df_res)], as.numeric)
  
  # Remove TL data (too correlated)
  df_res <- dplyr::filter(df_res, !pheno %in% c("TLc", "Zadj_TS"))
  
  # Add metadata
  df_res <- merge(df_res, metadata[,c("pheno","description","Category")], 
                  by = "pheno", all.x = TRUE)
  
  
  # Arrange in descending order of rcoef
  df_res <- arrange(df_res, (rcoef))
  
  write.table(df_res, res_file_path, row.names = F, quote = F, sep = '\t')
}