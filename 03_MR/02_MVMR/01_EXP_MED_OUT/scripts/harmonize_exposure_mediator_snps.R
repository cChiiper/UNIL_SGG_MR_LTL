
expo_df <- read.table(snakemake@input[["expo_df"]], header = T)
print("Exposure df is in memory.")
med_df <- read.table(snakemake@input[["med_df"]], header = T)
print("Mediator df is in memory.")
outcome_df <- read.table(snakemake@input[["outcome_df"]], header = T)
print("Outcome df is in memory..")

clumped_expo_snps <- read.table(snakemake@input[["clumped_expo_snps"]])
names(clumped_expo_snps) <- c('SNP')

med_snps <- read.table(snakemake@input[["med_snps"]])
names(med_snps) <- c('SNP')

## subset datasets to relevant SNPs
df <- expo_df[(expo_df$SNP %in% clumped_expo_snps$SNP) | (expo_df$SNP %in% med_snps$SNP), ]
df = merge(df, med_df, by = 'SNP', suffixes = c('_EXPO', '_MED'))
df = merge(df, outcome_df, by = 'SNP')

print("Dataframes are subsetted to relevant SNPs.")

# prepare dataset for mediator SNP clumping
clump_df = df[ , c("SNP", "p")]
clump_df[clump_df$SNP %in% clumped_expo_snps$SNP, "p"] = 1e-300 # prioritize exposure SNPs over mediator SNPs

write.table(clump_df, snakemake@output[["for_clumping"]], row.names = F, quote = F, sep = '\t')

## harmonize outcome alleles -> exposure SNPs will be reference SNPs
print(head(df))
mask_flip = ((as.character(df$A1_EXPO) == as.character(df$A2)) & (as.character(df$A2_EXPO) == as.character(df$A1)))
df[mask_flip, "A1"] = df[mask_flip, "A1_EXPO"]
df[mask_flip, "A2"] = df[mask_flip, "A2_EXPO"]
df[mask_flip, "Freq"] = 1-df[mask_flip, "Freq"]
df[mask_flip, "b"] = -df[mask_flip, "b"]

## check if remaining outcome alleles are correct
mask_correct = ((as.character(df$A1_EXPO) == as.character(df$A1)) & (as.character(df$A2_EXPO) == as.character(df$A2)))
df = df[mask_correct,]

## harmonize mediator alleles -> exposure SNPs will be reference SNPs
mask_flip = ((as.character(df$A1_EXPO) == as.character(df$A2_MED)) & (as.character(df$A2_EXPO) == as.character(df$A1_MED)))
df[mask_flip, "A1_MED"] = df[mask_flip, "A1_EXPO"]
df[mask_flip, "A2_MED"] = df[mask_flip, "A2_EXPO"]
df[mask_flip, "Freq_MED"] = 1-df[mask_flip, "Freq_MED"]
df[mask_flip, "b_MED"] = -df[mask_flip, "b_MED"]

## check if remaining outcome alleles are correct
mask_correct = ((as.character(df$A1_EXPO) == as.character(df$A1_MED)) & (as.character(df$A2_EXPO) == as.character(df$A2_MED)))
df = df[mask_correct,]

df <- df[,c("SNP", "A1_EXPO", "A2_EXPO", "Freq_EXPO", "Freq_MED", "Freq", "b_EXPO", "se_EXPO", "p_EXPO", "N_EXPO", "b_MED", "se_MED", "p_MED", "N_MED", "b", "se", "p", "N")]
names(df) <- c("SNP", "A1", "A2", "Freq_EXPO", "Freq_MED", "Freq_OUT", "b_EXPO", "se_EXPO", "p_EXPO", "N_EXPO", "b_MED", "se_MED", "p_MED", "N_MED", "b_OUT", "se_OUT", "p_OUT", "N_OUT")

write.table(df, snakemake@output[["harmonized"]], row.names = F, quote = F)