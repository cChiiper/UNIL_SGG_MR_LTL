################################################################################
### MR-APPS formatting                                                       ###
### Author: Samuel Moix                                                      ###
### Date: 28.02.2024                                                         ###
################################################################################

################################################
### Libraries ##################################
library(MRAPSS)
library(readr)

################################################
### Read data ##################################

exposure_gwas <- readr::read_delim(snakemake@input[["exp_gwas"]], "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
outcome_gwas <- readr::read_delim(snakemake@input[["out_gwas"]], "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

HLA_rsids <- read.table(snakemake@input[["HLA_snps"]], header = F)
names(HLA_rsids) = c('SNP')

HBB_rsids <- read.table(snakemake@input[["HBB_snps"]], header = F)
names(HBB_rsids) = c('SNP')


################################################
### Exclude rsids ##############################

print(paste0("Number of exposure SNPs: ", nrow(exposure_gwas)))
# Remove SNP from the HBB locus (chr11 5246696-5248301)
exposure_gwas <- subset(exposure_gwas, !SNP %in% HBB_rsids$SNP)
print(paste0("Number of exposure SNPs without HBB locus: ", nrow(exposure_gwas)))

# Remove SNPs from HLA locus (chr6 25MB-37MB)
exposure_gwas <- subset(exposure_gwas, !SNP %in% HLA_rsids$SNP)
print(paste0("Number of exposure SNPs without HLA locus: ", nrow(exposure_gwas)))

print(paste0("Number of outcome SNPs: ", nrow(outcome_gwas)))
# Remove SNP from the HBB locus (chr11 5246696-5248301)
outcome_gwas <- subset(outcome_gwas, !SNP %in% HBB_rsids$SNP)
print(paste0("Number of outcome SNPs without HBB locus: ", nrow(outcome_gwas)))
# Remove SNPs from HLA locus (chr6 25MB-37MB)
outcome_gwas <- subset(outcome_gwas, !SNP %in% HLA_rsids$SNP)
print(paste0("Number of outcome SNPs without HLA locus: ", nrow(outcome_gwas)))

################################################
### Format data ################################
exposure_format <- format_data(exposure_gwas,
                                snp_col = "SNP",
                                b_col = "b",
                                se_col = "se",
                                freq_col = "Freq",
                                A1_col = "A1",
                                A2_col = "A2",
                                p_col = "p",
                                n_col = "N")

outcome_format <- format_data(outcome_gwas,
                              snp_col = "SNP",
                              b_col = "b",
                              se_col = "se",
                              freq_col = "Freq",
                              A1_col = "A1",
                              A2_col = "A2",
                              p_col = "p",
                              n_col = "N")

### Harmonize data
paras <- est_paras(dat1 = exposure_format,
                   dat2 = outcome_format,
                   trait1.name = snakemake@wildcards[["expo_trait"]],
                   trait2.name = snakemake@wildcards[["outcome_trait"]],
                   ldscore.dir = snakemake@params[["ldsc_path"]])

### Print paras
print("dat:")
print(head(paras$dat))

print("ldsc_res:")
print(paras$ldsc_res)

print("C:")
print(paras$C)

print("Omega:")
print(paras$Omega)


################################################
### Clump ######################################

# @param dat  a data frame must have columns with information about SNPs and p values
# @param SNP_col column with SNP rsid. The default is `"SNP"`
# @param pval_col column with p value. The default is `"pval"`
# @param clump_kb clumping window in kb. Default is 1000.
# @param clump_r2 clumping r2 threshold. Default is 0.001.
# @param clump_p  clumping significance level for index variants. Default = 5e-05
# @param bfile     bfile as LD reference panel. If this is provided, then will use local PLINK. Default = NULL.
# @param plink_bin path to local plink binary. Default = NULL.

MRdat <- clump(paras$dat,
              IV.Threshold = 5e-05,
              SNP_col = "SNP",
              pval_col = "pval.exp",
              clump_kb = 1000,
              clump_r2 = 0.001,
              pop = "EUR",
              bfile = snakemake@params[["bfile_path"]],
              plink_bin = snakemake@params[["plink_path"]])


################################################
### Steiger filtering ##########################

print(paste0("Number of SNPs before Steiger filter: ", nrow(MRdat)))
# Calculate the Steiger z-value for each SNP
MRdat$zval_steiger <- (abs(MRdat$b.exp) - abs(MRdat$b.out)) / sqrt(MRdat$se.exp^2 + MRdat$se.out^2)
# Apply the Steiger filter to keep SNPs with zval_steiger > -1.96
MRdat <- MRdat[MRdat$zval_steiger > -1.96, ]
# Drop the zval_steiger column from the filtered dataset
MRdat <- within(MRdat, rm(zval_steiger))
# Display the number of SNPs after applying the Steiger filter
print(paste0("Number of SNPs after Steiger filter: ", nrow(MRdat)))


################################################
### Run MR-APPS ################################

MRres <- MRAPSS(MRdat,
               exposure = snakemake@wildcards[["expo_trait"]],
               outcome = snakemake@wildcards[["outcome_trait"]],
               C = paras$C,
               Omega = paras$Omega ,
               Cor.SelectionBias = T)

# names(MRres)
# "MRdat"        "exposure"     "outcome"      "beta"         "beta.se"     
# "pvalue"       "tau.sq"       "sigma.sq"     "pi0"          "post"        
# "IVsignal.sum" "likelihoods"  "Threshold"    "method"  

results_table <- data.frame(
  exposure = MRres$exposure,
  outcome = MRres$outcome,
  method = MRres$method,
  beta = MRres$beta,
  nsnp = nrow(MRdat),
  se = MRres$beta.se, 
  pvalue = MRres$pvalue,
  tau_sq = MRres$tau.sq, 
  sigma_sq = MRres$sigma.sq, 
  pi0 = MRres$pi0,
  IVsignal_sum = MRres$IVsignal.sum, 
  Threshold = MRres$Threshold
)

print("Results table:")
print(results_table)

################################################
### Save results ###############################

write.table(results_table, file = snakemake@output[["MRAPSS_res"]], col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
saveRDS(MRres, file = snakemake@output[["RDS_file"]])
saveRDS(paras, file = snakemake@output[["paras_RDS"]])

pdf(file = snakemake@output[["plot_file"]], width = 7, height = 5)
MRplot(MRres, exposure = snakemake@wildcards[["expo_trait"]], outcome = snakemake@wildcards[["outcome_trait"]])
dev.off()

################################################
### Sensitivity analysis #######################

# sensitivity(MRdat,
#             Omega=paras$Omega,
#             C = paras$C,
#             exposure = "TELOMERE",
#             outcome = "MCH")