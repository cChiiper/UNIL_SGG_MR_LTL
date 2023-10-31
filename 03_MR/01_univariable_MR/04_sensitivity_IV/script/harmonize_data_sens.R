#######################################
### R-script to harmonize SENS data ###
### Samuel Moix                     ###
### 08.03.2023                      ###
#######################################

################################################
### Load packages ##############################
library(data.table)
library(dplyr)
library(stringr)

################################################
### Load data and harmonize data ############### 
 
# Read snps
dt <- fread(snakemake@input[["snps"]], header = FALSE, col.names = "rsid")

### Load exposure data
exposure_df <- read.table(snakemake@input[["gwas_exp"]], header = T)
exposure_df <- exposure_df[exposure_df$SNP %in% dt$rsid,]
print(head(exposure_df))
print("Outcome df is in memory..")

### Load trait data
TRAITS <- snakemake@input[["gwas_trait"]]
counter_trait <- length(TRAITS)
for (TRAIT in TRAITS) {
  # Read each trait one after another  
  trait_df <- read.table(TRAIT, header = T)
  print(paste(TRAIT, "as trait df is in memory."))
    
  # Merge first trait with exposure df
  if(counter_trait == length(TRAITS)){
    df <- merge(exposure_df, trait_df, by = 'SNP', suffixes = c('_OUT', '_EXPO'), all = FALSE)
    rm(exposure_df)
  # Merge following traits
  }else{
    df <- merge(df, trait_df, by = 'SNP', suffixes = c(paste0('_',prev_EXP), '_EXPO'), all = FALSE)
    df <- df %>%
      rename(A1_EXPO = A1) %>%
      rename(A2_EXPO = A2)
    rm(trait_df)
  }
  print(paste0("Number of SNPs before allele harmonizing: ", nrow(df)))
  
  
  ### Harmonize alleles
  # Change columns data type 
  df$A1_EXPO <- as.character(df$A1_EXPO)
  df$A2_EXPO <- as.character(df$A2_EXPO)
  df$A1_OUT <- as.character(df$A1_OUT)
  df$A2_OUT <- as.character(df$A2_OUT)
  df$Freq_OUT <- as.numeric(df$Freq_OUT)
  df$b_OUT <- as.numeric(df$b_OUT)
  
  # Flip flipped alleles
  mask_flip <- ((df$A1_EXPO == df$A2_OUT) & (df$A2_EXPO == df$A1_OUT))
  df[mask_flip, "A1_OUT"] <- df[mask_flip, "A1_EXPO"]
  df[mask_flip, "A2_OUT"] <- df[mask_flip, "A2_EXPO"]
  df[mask_flip, "Freq_OUT"] <- 1-df[mask_flip, "Freq_OUT"] 
  df[mask_flip, "b_OUT"] <- -df[mask_flip, "b_OUT"]
  
  ### Check if remaining alleles are correct
  mask_correct <- ((as.character(df$A1_EXPO) == as.character(df$A1_OUT)) & (as.character(df$A2_EXPO) == as.character(df$A2_OUT)))
  df <- df[mask_correct,]
  
  print(paste0("Number of SNPs of allele harmonizing: ", nrow(df)))
  # Remove added allele reference double
  df <- select(df, -c("A1_EXPO","A2_EXPO"))
  # Rename suffixes
  trait_name <- stringr::str_extract(TRAIT, "(?<=_to_).+(?=_gwas_subset\\.ma)") # Could be better !!!!!!
  if(counter_trait > 1){
    colnames(df) = gsub("_EXPO", "", colnames(df))
    prev_EXP <- trait_name
    counter_trait <- counter_trait-1
  # Rename last suffixe
  }else{
    colnames(df) = gsub("_EXPO", paste0("_",trait_name), colnames(df))
  }
}

### Rename allele columns
df <- df %>% 
  rename(A1 = A1_OUT) %>%
  rename(A2 = A2_OUT) 

################################################
### Export results ############################# 

write.table(df, snakemake@output[[1]], row.names = F, quote = F)
write.table(df$SNP, snakemake@output[[2]], row.names = F, quote = F, col.names = F)

print("Harmonize data done!")