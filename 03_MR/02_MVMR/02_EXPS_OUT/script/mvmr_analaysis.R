#####################################
### R-script to run MVMR            #
### Samuel Moix                     #
### 10.03.2023                      #
#####################################

################################################
### Setup working space ########################

### Loading data
df = read.table(snakemake@input[["data"]], header = T)
print("Input data:")
print(head(df))

# Extract exposures from column names [-1] to remove output
EXPOSURES <- sub("^b_","",grep("^b_", names(df), value = T))[-1]
print("List of exposures")
print(EXPOSURES)

snps = read.table(snakemake@input[["clumped_snps"]])
names(snps) = c('SNP')

HLA_rsids <- read.table(snakemake@input[["HLA_snps"]], header = F)
names(HLA_rsids) = c('SNP')

HBB_rsids <- read.table(snakemake@input[["HBB_snps"]], header = F)
names(HBB_rsids) = c('SNP')

### Functions
 '%ni%' <- Negate("%in%")

################################################
### Filtering SNPS ############################# 

# Keep only clumped SNPs
df <- df[df$SNP %in% snps$SNP, ]
print(paste0("Number of clumped SNPs: ", nrow(df)))

### Remove SNPs from HLA locus (chr6 25MB-37MB)
df <- df[df$SNP %ni% HLA_rsids$SNP, ]

### Remove SNPs from HBB locus
df <- df[df$SNP %ni% HBB_rsids$SNP, ]

print(paste0("Number of SNPs without HLA and HBB locus SNPs: ", nrow(df)))

### Allele frequency filter 
for (exp_i in EXPOSURES) {
  df <- df[abs(df[, paste0("Freq_",exp_i)] - df$Freq_OUT) < 0.05,]
}

print(paste0("Number of SNPs after frequency check: ", nrow(df)))

# ### Steiger filter: remove SNPs with larger outcome than exposure effects (no mask)
# for (exp_i in EXPOSURES) {
#   df$zval_steiger <- (abs(df[, paste0("b_",exp_i)])-abs(df$b_OUT))/sqrt(df[, paste0("se_",exp_i)]**2 + df$se_OUT**2)
#   df <- df[df$zval_steiger > -1.96, ]
# }

### Steiger filter: remove SNPs with larger outcome than exposure effects with mask
df$keep_steiger <- TRUE

# Subset snps for MVMR (no exposure to exposure filtering)
df_MVMR <- df

for (exp_i in EXPOSURES) {
  # Exposure to outcome
  mask_exp_i = (df_MVMR[, paste0("p_", exp_i)] < 5e-8)
  df_MVMR$zval_expo_out_steiger <- (abs(df_MVMR[, paste0("b_", exp_i)]) - abs(df_MVMR$b_OUT)) / sqrt(df_MVMR[, paste0("se_", exp_i)]**2 + df_MVMR$se_OUT**2)
  df_MVMR$keep_steiger <- df_MVMR$keep_steiger & ((mask_exp_i & (df_MVMR$zval_expo_out_steiger > -1.96)) | !mask_exp_i)
}

df_MVMR <- df_MVMR[df_MVMR$keep_steiger, ]
print(paste0("Number of SNPs after steiger filter for MVMR: ", nrow(df_MVMR)))

# Subset snps for IVW
for (exp_i in EXPOSURES) {
  # Exposure to outcome
  mask_exp_i = (df[, paste0("p_", exp_i)] < 5e-8)
  df$zval_expo_out_steiger <- (abs(df[, paste0("b_", exp_i)]) - abs(df$b_OUT)) / sqrt(df[, paste0("se_", exp_i)]**2 + df$se_OUT**2)
  df$keep_steiger <- df$keep_steiger & ((mask_exp_i & (df$zval_expo_out_steiger > -1.96)) | !mask_exp_i)
  
  # Exposure to exposure
  for (exp_j in EXPOSURES) {
    if (exp_i != exp_j) {
      mask_exp_j <- (df[, paste0("p_", exp_j)] < 5e-8)
      df$zval_steiger <- (abs(df[, paste0("b_", exp_i)]) - abs(df[, paste0("b_", exp_j)])) / sqrt(df[, paste0("se_", exp_i)]**2 + df[, paste0("se_", exp_j)]**2)
      df$keep_steiger <- df$keep_steiger & ((mask_exp_i & (df$zval_steiger > -1.96) & !mask_exp_j) | !mask_exp_i)
    }
  }
}

df <- df[df$keep_steiger, ]
print(paste0("Number of SNPs after steiger filter for IVW: ", nrow(df)))

################################################
### Analysis ################################### 

### MVMR regression
formula <- as.formula(paste("b_OUT ~ ", paste(paste0("b_", EXPOSURES), collapse = " + "),"-1"))
multi.res <- summary(lm(formula, weights = 1/se_OUT^2, data = df_MVMR))

# Get results for MVMR
results_df <- as.data.frame(multi.res$coefficients)
results_df$exposure <- sub("b_","",rownames(results_df))
colnames(results_df) <- c("b", "se", "t_value", "pval","exposure")
results_df$relation <- "MVMR"
results_df$nb_snp <- nrow(df_MVMR)

### Q-heterogeneity statistic

# Calculate the number of mediators (K) and the number of SNPs (L)
K <- length(EXPOSURES)
print(paste("Number of exposures (K):", K))
L <- nrow(df_MVMR)
print(paste("Number of SNPs (L):", L))

# Extract the relevant columns from the df_MVMR data frame
y <- as.vector(df_MVMR$b_OUT)
se_y <- as.vector(df_MVMR$se_OUT)
X <- as.matrix(df_MVMR[, paste0("b_", EXPOSURES)])
print("Head of exposure estimates matrix:")
print(head(X))

# Run IVW regression
ivw.res <- summary(lm(y ~ -1 + X, weights = 1/se_y^2))

# Calculate Q-heterogeneity statistics
Q_df <- L - K # Degrees of freedom
Q <- ivw.res$sigma^2 * Q_df #sigma is the residual standard error
Q_pval <- pchisq(Q, Q_df, lower.tail = FALSE)

# Add Q-heterogeneity statistics to the results data frame
results_df$Q <- Q
results_df$Q_pval <- Q_pval

# Clean
rm(Q_df,Q,Q_pval,y,se_y,ivw.res)

### F-conditional statistic
results_df$F_stat <- NA
for (exp_i in EXPOSURES) {
  # SE matrix
  SE <- as.matrix(df_MVMR[, paste0("se_", EXPOSURES)])
  # E-xposure 
  E <- as.vector(X[,paste0("b_",exp_i)])
  # Other exposures / M-ediators
  M <- as.matrix(X[,paste0("b_", EXPOSURES[-which(EXPOSURES %in% exp_i)])]) # without exp_i 
  # Do the magic  https://doi.org/10.1002/sim.9133
  delta.reg <- lm(E ~ -1 + M)
  delta <- delta.reg$coefficients
  res <- delta.reg$residuals
  delta <- as.vector(c(1,delta))
  var_IV <-  (SE**2) %*% (delta**2)
  Q <- sum(1/(var_IV) * res**2)
  Fstat <- Q/(L-(K-1))
  results_df$F_stat[results_df$exposure == exp_i] <- Fstat
  print(paste("Added F-stat for:", exp_i,":", Fstat))
}

# Clean
rm(Q,K,L,X,SE,delta.reg,delta,res,var_IV)

### MR regression
# Get ivw results for all exposures stringent steiger
for (exp_i in EXPOSURES) {
  # Select SNPs with p-value < 5e-8
  exp_i_df <- df[df[,paste0("p_", exp_i)] < 5e-8,]
  # Run IVW regression
  ivw.res <- summary(lm(as.formula(paste("b_OUT ~ -1 +", paste0("b_",exp_i))),
                        weights = 1/se_OUT^2, data = exp_i_df))
  # Save results in temporary dataframe
  temp_df <- as.data.frame(ivw.res$coefficients)
  temp_df$exposure <- sub("b_","",rownames(temp_df))
  colnames(temp_df) <- c("b", "se", "t_value", "pval","exposure")
  temp_df$relation <- "ivw_steiger"
  temp_df$nb_snp <- nrow(exp_i_df)
  # Calculate Q-heterogeneity statistics
  Q_df <- nrow(exp_i_df) - 1 # Degrees of freedom
  Q <- ivw.res$sigma^2 * Q_df
  Q_pval <- pchisq(Q, Q_df, lower.tail = FALSE)
  # Add Q-heterogeneity statistics to the temporary results data frame
  temp_df$Q <- Q
  temp_df$Q_pval <- Q_pval
  temp_df$F_stat <- NA
  # Add results to final dataframe
  results_df <- rbind(results_df,temp_df)
  rownames(results_df) <- NULL
}
rm(temp_df)

### IVW for all exposures normal
# Get ivw results for all exposures
for (exp_i in EXPOSURES) {
  # Select SNPs with p-value < 5e-8
  exp_i_df <- df_MVMR[df_MVMR[,paste0("p_", exp_i)] < 5e-8,]
  # Run IVW regression
  ivw.res <- summary(lm(as.formula(paste("b_OUT ~ -1 +", paste0("b_",exp_i))),
                        weights = 1/se_OUT^2, data = exp_i_df))
  # Save results in temporary dataframe
  temp_df <- as.data.frame(ivw.res$coefficients)
  temp_df$exposure <- sub("b_","",rownames(temp_df))
  colnames(temp_df) <- c("b", "se", "t_value", "pval","exposure")
  temp_df$relation <- "ivw"
  temp_df$nb_snp <- nrow(exp_i_df)
  # Calculate Q-heterogeneity statistics
  Q_df <- nrow(exp_i_df) - 1 # Degrees of freedom
  Q <- ivw.res$sigma^2 * Q_df
  Q_pval <- pchisq(Q, Q_df, lower.tail = FALSE)
  # Add Q-heterogeneity statistics to the temporary results data frame
  temp_df$Q <- Q
  temp_df$Q_pval <- Q_pval
  temp_df$F_stat <- NA
  # Add results to final dataframe
  results_df <- rbind(results_df,temp_df)
  rownames(results_df) <- NULL
}
rm(temp_df)

# Add outcome information in result
results_df$outcome <- snakemake@params[["out_name"]]

# Save results dataframe
results_df <- results_df[,c("exposure", "outcome", "b", "se", "pval","relation","nb_snp","Q","Q_pval","F_stat")]
print("Results:")
print(results_df)

write.table(results_df, snakemake@output[[1]], row.names = F, quote = F, sep = '\t')
