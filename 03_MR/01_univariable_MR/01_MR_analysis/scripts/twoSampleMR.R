source(paste0(getwd(), "/scripts/mr_functions.R"))

df <- read.table(snakemake@input[[1]], header = T)
print(paste0("Number of IVs: ", nrow(df)))

## allele frequency filter
df <- df[abs(df$Freq_EXPO - df$Freq_OUT) < 0.05,]
print(paste0("Number of SNPs after frequency check: ", nrow(df)))

## Steiger filter: remove SNPs with larger outcome than exposure effects

df$zval_steiger <- (abs(df$b_EXPO)-abs(df$b_OUT))/sqrt(df$se_EXPO**2 + df$se_OUT**2)
df <- df[df$zval_steiger > -1.96, ]
print(paste0("Number of SNPs after Steiger filter: ", nrow(df)))

# MR analyses
result_df <- data.frame(matrix(NA, 3, 7))
names(result_df) <- c('method', 'b_mr', 'se_mr', 'pval_mr', 'nsnps', 'Q_stat', 'Q_pval')

res <- mr_ivw(df$b_EXPO, df$b_OUT, df$se_EXPO, df$se_OUT)
result_df[1, 'method'] <- 'mr_ivw'
result_df[1, 'b_mr'] <- res$b[[1]]
result_df[1, 'se_mr'] <- res$se[[1]]
result_df[1, 'pval_mr'] <- res$pval[[1]]
result_df[1, 'nsnps'] <- res$nsnp[[1]]
result_df[1, 'Q_stat'] <- res$Q[[1]]
result_df[1, 'Q_pval'] <- res$Q_pval[[1]]

res <- mr_simple_median(df$b_EXPO, df$b_OUT, df$se_EXPO, df$se_OUT, 1000)
result_df[2, 'method'] <- 'mr_simple_median'
result_df[2, 'b_mr'] <- res$b[[1]]
result_df[2, 'se_mr'] <- res$se[[1]]
result_df[2, 'pval_mr'] <- res$pval[[1]]
result_df[2, 'nsnps'] <- res$nsnp[[1]]

res <- mr_simple_mode(df$b_EXPO, df$b_OUT, df$se_EXPO, df$se_OUT)
result_df[3, 'method'] <- 'mr_simple_mode'
result_df[3, 'b_mr'] <- res$b[[1]]
result_df[3, 'se_mr'] <- res$se[[1]]
result_df[3, 'pval_mr'] <- res$pval[[1]]
result_df[3, 'nsnps'] <- res$nsnp[[1]]

write.table(result_df, snakemake@output[[1]], row.names = F, quote = F)