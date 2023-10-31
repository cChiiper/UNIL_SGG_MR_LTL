### Prepare MR data
expo_df <- read.table(snakemake@input[["expo_df"]], header = T)
print("Exposure df is in memory.")
outcome_df <- read.table(snakemake@input[["outcome_df"]], header = T)
print("Outcome df is in memory..")

snps <- read.table(snakemake@input[["snps"]])
names(snps) <- c('SNP')

expo_df <- expo_df[expo_df$SNP %in% snps$SNP,]
print("Exposure df is subsetted to clumped SNPs.")

df <- merge(expo_df, outcome_df, by = 'SNP', suffixes = c('_EXPO', '_OUT'))
print(paste0("Number of SNPs before allele harmonizing: ", nrow(df)))

## harmonize alleles
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

## check if remaining alleles are correct
mask_correct <- ((as.character(df$A1_EXPO) == as.character(df$A1_OUT)) & (as.character(df$A2_EXPO) == as.character(df$A2_OUT)))
df <- df[mask_correct,]

print(paste0("Number of SNPs of allele harmonizing: ", nrow(df)))

df <- df[,c("SNP", "A1_EXPO", "A2_EXPO", "Freq_EXPO", "Freq_OUT", "b_EXPO", "se_EXPO", "p_EXPO", "N_EXPO", "b_OUT", "se_OUT", "p_OUT", "N_OUT")]
names(df) <- c("SNP", "A1", "A2", "Freq_EXPO", "Freq_OUT", "b_EXPO", "se_EXPO", "p_EXPO", "N_EXPO", "b_OUT", "se_OUT", "p_OUT", "N_OUT")

write.table(df, snakemake@output[[1]], row.names = F, quote = F)