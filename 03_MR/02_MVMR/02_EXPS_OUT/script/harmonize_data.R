#####################################
### R-script to harmonize MVMR data #
### Samuel Moix                     #
### 08.03.2023                      #
#####################################

### Load packages
library(data.table)
library(dplyr)

### Merge clumped SNP files
# Create an empty data.table to store the unique IDs
unique_rsids <- data.table(rsid = character())

# Loop through all the files in the directory
for (filename in snakemake@input[["snps"]]) {
  # Read the file into a data.table
  dt <- fread(filename, header = FALSE, col.names = "rsid")
  
  # Add the IDs from this file to the unique_ids data.table
  unique_rsids <- unique(rbind(unique_rsids, dt))
}

# Print the number of unique rsids
cat("There are", nrow(unique_rsids), "unique SNPs.\n")


################################################
### Prepare MVMR data ########################## 
 
### Load outcome data
outcome_df <- read.table(snakemake@input[["gwas_out"]], header = T)
outcome_df <- outcome_df[outcome_df$SNP %in% unique_rsids$rsid,]
print(head(outcome_df))
print("Outcome df is in memory..")

### Load exposure data
EXPOSURES <- snakemake@input[["gwas_exp"]]
counter_exposure <- length(EXPOSURES)
for (EXPOSURE in EXPOSURES) {
  # Read each exposure one after another  
  expo_df <- read.table(EXPOSURE, header = T)
  print(paste(EXPOSURE, "as exposure df is in memory."))
  
  # Subset rsids
  expo_df <- expo_df[expo_df$SNP %in% unique_rsids$rsid,]
  print("Exposure df is subsetted to clumped SNPs.")
  
  # Merge first exposure with outcome df
  if(counter_exposure == length(EXPOSURES)){
    df <- merge(outcome_df, expo_df, by = 'SNP', suffixes = c('_OUT', '_EXPO'))
    rm(outcome_df)
  # Merge following exposures
  }else{
    df <- merge(df, expo_df, by = 'SNP', suffixes = c(paste0('_',prev_EXP), '_EXPO'))
    df <- df %>%
      rename(A1_EXPO = A1) %>%
      rename(A2_EXPO = A2)
    rm(expo_df)
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
  exp_name <- dplyr::last(as.vector(strsplit(strsplit(EXPOSURE, "_gwas_summary_uk10kck.ma")[[1]][1], "/")[[1]])) # Could be better
  if(counter_exposure > 1){
    colnames(df) = gsub("_EXPO", "", colnames(df))
    prev_EXP <- exp_name
    counter_exposure <- counter_exposure-1
  # Rename last suffixe
  }else{
    colnames(df) = gsub("_EXPO", paste0("_",exp_name), colnames(df))
  }
}

### Rename allele columns
df <- df %>% 
  rename(A1 = A1_OUT) %>%
  rename(A2 = A2_OUT) 

################################################
### Prepare clumping file based on ranked method 

df_rank <- df %>%
  select(c("SNP",contains("p_"))) %>% # Select columns
  select(-c("p_OUT")) %>% # Remove p-value outcome column
  mutate_if(is.numeric, list(~dense_rank(.))) # Replace p-value by rank

df_rank <- df_rank %>%
  mutate(ranksum = rowSums(.[2:ncol(df_rank)])) %>% # Sum ranks, new dplyr across
  mutate(rankmin = apply(.[2:ncol(df_rank)], 1, min, na.rm = TRUE)) %>% # Minimum rank
  select(c("SNP","ranksum","rankmin")) %>%
  mutate(ranksum = ranksum/nrow(df_rank)/100) %>% # Divide by 100 for plink
  mutate(rankmin = rankmin/nrow(df_rank)/10) # Divide by 10 for plink

################################################
### Export results ############################# 

write.table(df, snakemake@output[[1]], row.names = F, quote = F)
write.table(df_rank, snakemake@output[[2]], row.names = F, quote = F)

print("Harmonize data done!")