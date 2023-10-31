################################################
### R-script for steiger filtering SENS      ### 
### Samuel Moix                              ###
### 19.04.2023                               ###
################################################

################################################
### Load packages ##############################
library(data.table)
library(dplyr)

################################################
### Load data  ################################# 

# Load harmonized data
df <- read.table(snakemake@input[["data"]], header = T)
# Load clumped snps
df_snps <- fread(snakemake@input[["snps"]], header = FALSE, col.names = "rsid")
# Subset clumped snps
df <- df[df$SNP %in% df_snps$rsid,]
print(paste0("Number of clumped SNPs: ", nrow(df)))

# Extract exposures from column names [-1] to remove output
EXPOSURES <- sub("^b_","",grep("^b_", names(df), value = T))[-1]
print("List of exposures")
print(EXPOSURES)

HLA_rsids <- read.table(snakemake@input[["HLA_snps"]], header = F)
names(HLA_rsids) = c('SNP')

HBB_rsids <- read.table(snakemake@input[["HBB_snps"]], header = F)
names(HBB_rsids) = c('SNP')

### Functions
 '%ni%' <- Negate("%in%")

################################################
### Filtering SNPS ############################# 

### Remove SNPs from HLA locus (chr6 25MB-37MB)
df <- df[df$SNP %ni% HLA_rsids$SNP, ]

### Remove SNPs from HBB locus
df <- df[df$SNP %ni% HBB_rsids$SNP, ]

print(paste0("Number of SNPs without HLA and HBB locus SNPs: ", nrow(df)))

### Run steiger filter

### Steiger filter: remove SNPs with larger trait than exposure effects with mask
for (exp_i in EXPOSURES) {
  df$zval_steiger <- (abs(df$b_OUT)-abs(df[, paste0("b_",exp_i)]))/sqrt(df[, paste0("se_",exp_i)]**2 + df$se_OUT**2)
  df <- df[df$zval_steiger > -1.96, ]
}

print(paste0("Number of SNPs after exposure --> trait steiger filtering: ", nrow(df)))

################################################
### Export results ############################# 

write.table(df$SNP, snakemake@output[[1]], row.names = F, quote = F, col.names = F)
