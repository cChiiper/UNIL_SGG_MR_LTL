################################################################################
### MR Plots from MR prepared data                                           ###
### Author: Samuel Moix                                                      ###
### Date: 29.11.2023                                                         ###
################################################################################


################################################
### Libraries ##################################
library(TwoSampleMR) #TwoSampleMR version 0.5.6
library(ggplot2)
library(dplyr)
library(gridExtra)
library(data.table)
library(tidyr)
library(conflicted)

################################################
### Rsids to exclude files #####################
HLA_rsids <- read.table(snakemake@input[["HLA_snps"]], header = F)
names(HLA_rsids) = c('SNP')

HBB_rsids <- read.table(snakemake@input[["HBB_snps"]], header = F)
names(HBB_rsids) = c('SNP')

CHIPS_rsid <- read.table(snakemake@input[["CHIP_snps"]], header = F)
names(CHIPS_rsid) = c('SNP')

################################################
### TwoSampleMR ################################

EXPOSURE <- snakemake@params[["EXPNAME"]]
OUTCOME <- snakemake@params[["OUTNAME"]]
    
### Read exposure data
exposure_dat <- read_exposure_data(
    file.path(snakemake@input[["mr_data"]]),
                                   sep = " ",
                                   snp_col = "SNP", 
                                   effect_allele_col = "A1",
                                   other_allele_col = "A2",
                                   eaf_col = "Freq_EXPO",
                                   beta_col = "b_EXPO",
                                   se_col = "se_EXPO",
                                   pval_col = "p_EXPO",
                                   samplesize_col = "N_EXPO")

exposure_dat$exposure <- EXPOSURE

### Read outcome data
outcome_dat <- read_outcome_data(file.path(snakemake@input[["mr_data"]]),
                                 sep = " ",
                                 snp_col = "SNP", 
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2",
                                 eaf_col = "Freq_OUT",
                                 beta_col = "b_OUT",
                                 se_col = "se_OUT",
                                 pval_col = "p_OUT",
                                 samplesize_col = "N_OUT")

outcome_dat$outcome <- OUTCOME

### Harmonize data
dat <- harmonise_data(exposure_dat, outcome_dat)
# In cases where eaf = 0 NA's are introduced
if(sum(rowSums(is.na(dat)) > 0) != 0){
  print("Removing the following SNPs (eaf = 0)")
  print(dat[rowSums(is.na(dat)) > 0, "SNP"]) 
}

# Here we remove such cases
dat <- dat[rowSums(is.na(dat)) == 0,]

### allele frequency filter
print(paste0("Starting with ", nrow(dat), " SNPs"))
dat <- dat[abs(dat$eaf.exposure - dat$eaf.outcome) < 0.05,]
print(paste0("Number of SNPs after frequency check: ", nrow(dat)))

### Steiger filter: remove SNPs with larger outcome than exposure effects
dat$zval_steiger <- (abs(dat$beta.exposure)-abs(dat$beta.outcome))/sqrt(dat$se.exposure**2 + dat$se.outcome**2)
dat <- dat[dat$zval_steiger > -1.96, ]
print(paste0("Number of SNPs after steiger filter: ", nrow(dat)))

### With TwoSampleMR function
# dat <- steiger_filtering(dat)
# dat <- dat[(dat$steiger_dir) | (!dat$steiger_dir & dat$steiger_pval > 0.05), ]
# print(paste0("Number of SNPs after Steiger filter: ", nrow(dat)))

### Remove SNP from the HBB locus (chr11 5246696-5248301)
dat <- dplyr::filter(dat, !SNP %in% HBB_rsids$SNP)
print(paste0("Number of SNPs without HBB locus: ", nrow(dat)))
    
### Remove SNPs from HLA locus (chr6 25MB-37MB)
dat <- dplyr::filter(dat, !SNP %in% HLA_rsids$SNP)
print(paste0("Number of SNPs without HLA locus: ", nrow(dat)))

### Remove SNPs from CHIPS genes (23)
dat <- dplyr::filter(dat, !SNP %in% CHIPS_rsid$SNP)
print(paste0("Number of SNPs without CHIPS genes: ", nrow(dat)))

### MR analysis
mr_results <- mr(dat)

### MR results to export
mr_add_res_table <- mr_results %>% select(c("exposure", "outcome", "method", "nsnp", 
                                            "b", "se", "pval"))
if(nrow(dat) >= 2){ # If only one SNP can't compute heterogeneity
  q_stat <- mr_heterogeneity(dat)[,c("method","Q","Q_pval")]
} else{
  q_stat <- data.frame("method" = c("MR Egger","Inverse variance weighted"), 
                       "Q" = c(NA, NA),
                       "Q_pval" = c(NA,NA))
}
mr_add_res_table <- merge(mr_add_res_table, q_stat, by="method", all.x = T)
mr_add_res_table <- mr_add_res_table[,c("exposure", "outcome", "method", "nsnp", 
                                        "b", "se", "pval","Q","Q_pval")]

### MR-Egger additional statistics
mr_egger_res <- mr_egger_regression(b_exp = dat$beta.exposure, 
                                    b_out = dat$beta.outcome, 
                                    se_exp = dat$se.exposure,
                                    se_out = dat$se.outcome)

mr_add_res_table$Egger_intercept <- NA
mr_add_res_table$Egger_intercept[which(mr_add_res_table$method == "MR Egger")] <- mr_egger_res$b_i

mr_add_res_table$E_se_i <- NA
mr_add_res_table$E_se_i[which(mr_add_res_table$method == "MR Egger")] <- mr_egger_res$se_i

mr_add_res_table$E_pval_i <- NA
mr_add_res_table$E_pval_i[which(mr_add_res_table$method == "MR Egger")] <- mr_egger_res$pval_i


################################################
### Save results ###############################
fwrite(mr_add_res_table, snakemake@output[[1]], col.names = T, row.names = F, quote = F, sep = "\t")
