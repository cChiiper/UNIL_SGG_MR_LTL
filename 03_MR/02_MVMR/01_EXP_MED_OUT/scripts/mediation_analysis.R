
df = read.table(snakemake@input[["data"]], header = T)
snps = read.table(snakemake@input[["snps"]])
names(snps) = c('SNP')

df = df[df$SNP %in% snps$SNP, ]

print(paste0("Number of clumped SNPs: ", nrow(df)))

## Removing specific snps
HLA_snps = read.table(snakemake@input[["HLA_snps"]])
HBB_snps = read.table(snakemake@input[["HBB_snps"]])
names(HLA_snps) = c('SNP')
names(HBB_snps) = c('SNP')

df = df[!df$SNP %in% HLA_snps$SNP, ]
df = df[!df$SNP %in% HBB_snps$SNP, ]

print(paste0("Number of SNPs after removing HLA/HBB snps: ", nrow(df)))

## frequency difference check
df = df[abs(df$Freq_OUT - df$Freq_EXPO) < 0.05,]
df = df[abs(df$Freq_OUT - df$Freq_MED) < 0.05,]
df = df[abs(df$Freq_EXPO - df$Freq_MED) < 0.05,]

print(paste0("Number of SNPs after frequency check: ", nrow(df)))

## Steiger filtering

# Exposure to outcome
mask_expo = (df$p_EXPO < 5e-8)
df$zval_expo_steiger = (abs(df$b_EXPO)-abs(df$b_OUT))/sqrt(df$se_EXPO**2 + df$se_OUT**2)
df = df[(mask_expo & (df$zval_expo_steiger > -1.96)) | !mask_expo, ]

df_expo = df[df$p_EXPO < 5e-8, ] # for total effect calculation

# Exposure to mediator
mask_expo = (df$p_EXPO < 5e-8)
df$zval_expo_med_steiger = (abs(df$b_EXPO)-abs(df$b_MED))/sqrt(df$se_EXPO**2 + df$se_MED**2)
df = df[(mask_expo & (df$zval_expo_med_steiger > -1.96)) | !mask_expo, ]

# Mediator to outcome
mask_med = (df$p_MED < 5e-8)
mask_expo = (df$p_EXPO < 5e-8)
df$zval_med_steiger = (abs(df$b_MED)-abs(df$b_OUT))/sqrt(df$se_OUT**2 + df$se_MED**2)
df = df[(mask_med & (df$zval_med_steiger > -1.96) & !mask_expo) | !mask_med, ]

print(paste0("Number of SNPs after Steiger filtering: ", nrow(df)))

# calculate total effect
ivw.res <- summary(lm(b_OUT ~ -1 + b_EXPO, weights = 1/se_OUT^2, data = df_expo))
b_tot <- ivw.res$coef["b_EXPO","Estimate"]
se_tot <- ivw.res$coef["b_EXPO","Std. Error"]/min(1,ivw.res$sigma) #sigma is the residual standard error
p_tot <- 2 * pnorm(abs(b_tot/se_tot), lower.tail=FALSE)
n_tot = nrow(df_expo)

# calculate exposure to mediator effect
ivw.res <- summary(lm(b_MED ~ -1 + b_EXPO, weights = 1/se_MED^2, data = df_expo))
b_xm <- ivw.res$coef["b_EXPO","Estimate"]
se_xm <- ivw.res$coef["b_EXPO","Std. Error"]/min(1,ivw.res$sigma) #sigma is the residual standard error

# calculate direct effects
multi.res = summary(lm(b_OUT ~ b_EXPO + b_MED -1, weights = 1/se_OUT^2, data = df))

b_d = multi.res$coef["b_EXPO", "Estimate"]
se_d = multi.res$coef["b_EXPO", "Std. Error"] 
p_d = multi.res$coef["b_EXPO", "Pr(>|t|)"]

b_med_d = multi.res$coef["b_MED", "Estimate"]
se_med_d = multi.res$coef["b_MED", "Std. Error"]

n_d = nrow(df)

# calculate indirect effect - method A
b_ind_A = b_tot - b_d
se_ind_A = sqrt(se_tot^2 + se_d^2)
p_ind_A = 2*pnorm(abs(b_ind_A/se_ind_A), lower.tail = FALSE)

# calculate indirect effect - method B
b_ind = b_xm*b_med_d
se_ind = sqrt(se_med_d^2*se_xm^2 + b_xm^2*se_med_d^2 + b_med_d^2*se_xm^2)
p_ind = 2*pnorm(abs(b_ind/se_ind), lower.tail = FALSE)

### Account for regression dilution bias

#1 set non-significant effects to zero
df_shrinked = df
df_shrinked[df_shrinked$p_EXPO > 5e-8, "b_EXPO"] = 0
multi.res = summary(lm(b_OUT ~ b_EXPO + b_MED -1, weights = 1/se_OUT^2, data = df_shrinked))

b_d_shrink = multi.res$coef["b_EXPO", "Estimate"]
se_d_shrink = multi.res$coef["b_EXPO", "Std. Error"] 
p_d_shrink = multi.res$coef["b_EXPO", "Pr(>|t|)"]

b_med_d_shrink = multi.res$coef["b_MED", "Estimate"]
se_med_d_shrink = multi.res$coef["b_MED", "Std. Error"]

# calculate indirect effect - method A
b_ind_A_shrink = b_tot - b_d_shrink
se_ind_A_shrink = sqrt(se_tot^2 + se_d_shrink^2)
p_ind_A_shrink = 2*pnorm(abs(b_ind_A_shrink/se_ind_A_shrink), lower.tail = FALSE)

# calculate indirect effect - method B
b_ind_shrink = b_xm*b_med_d_shrink
se_ind_shrink = sqrt(se_med_d_shrink^2*se_xm^2 + b_xm^2*se_med_d_shrink^2 + b_med_d_shrink^2*se_xm^2)
p_ind_shrink = 2*pnorm(abs(b_ind_shrink/se_ind_shrink), lower.tail = FALSE)

# create dataframe with output

result_df = data.frame(matrix(NA, 7, 5))
names(result_df) = c('Relation', 'beta', 'se', 'pval', 'nsnps')
result_df[1,1] = "Total effect"
result_df[2,1] = "Direct effect"
result_df[3,1] = "Indirect effect - A"
result_df[4,1] = "Indirect effect - B"
result_df[5,1] = "Direct effect - Shrinkage"
result_df[6,1] = "Indirect effect - Shrinkage - A"
result_df[7,1] = "Indirect effect - Shrinkage - B"

result_df[1,2] = b_tot
result_df[1,3] = se_tot
result_df[1,4] = p_tot
result_df[1,5] = n_tot

result_df[2,2] = b_d
result_df[2,3] = se_d
result_df[2,4] = p_d
result_df[2,5] = n_d

result_df[3,2] = b_ind_A
result_df[3,3] = se_ind_A
result_df[3,4] = p_ind_A
result_df[3,5] = n_d

result_df[4,2] = b_ind
result_df[4,3] = se_ind
result_df[4,4] = p_ind
result_df[4,5] = n_d

result_df[5,2] = b_d_shrink
result_df[5,3] = se_d_shrink
result_df[5,4] = p_d_shrink
result_df[5,5] = n_d

result_df[6,2] = b_ind_A_shrink
result_df[6,3] = se_ind_A_shrink
result_df[6,4] = p_ind_A_shrink
result_df[6,5] = n_d

result_df[7,2] = b_ind_shrink
result_df[7,3] = se_ind_shrink
result_df[7,4] = p_ind_shrink
result_df[7,5] = n_d

### Adding exposure, outcome, and mediator information in each result files
# 1. Extract the wildcard values
expo_value <- snakemake@wildcards[["expo_trait"]]
med_value <- snakemake@wildcards[["med_trait"]]
outcome_value <- snakemake@wildcards[["outcome_trait"]]

# 2. Expand the dataframe
result_df$exposure <- NA
result_df$mediator <- NA
result_df$outcome <- NA

# 3. Assign the values to the new columns
result_df$exposure <- expo_value
result_df$mediator <- med_value
result_df$outcome <- outcome_value

write.table(result_df, snakemake@output[[1]], row.names = F, quote = F, sep = '\t')



