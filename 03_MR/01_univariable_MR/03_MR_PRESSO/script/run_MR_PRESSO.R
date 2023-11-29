################################################################################
### MR Plots from MR prepared data                                           ###
### Author: Samuel Moix                                                      ###
### Date: 16.10.2023                                                         ###
################################################################################

################################################
### Libraries ##################################
library(TwoSampleMR) #TwoSampleMR version 0.5.6
library(dplyr)
library(data.table)
library(tidyr)
library(conflicted)

################################################
### Functions ##################################
round_up_1000 <- function(x) {
  return(ceiling(x / 1000) * 1000)
}

################################################
### Rsids to exclude files #####################
HLA_rsids <- read.table(snakemake@input[["HLA_snps"]], header = F)
names(HLA_rsids) = c('SNP')

HBB_rsids <- read.table(snakemake@input[["HBB_snps"]], header = F)
names(HBB_rsids) = c('SNP')

################################################
### TwoSampleMR-PRESSO #########################

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
    samplesize_col = "N_EXPO"
)
exposure_dat$exposure <- EXPOSURE
    
### Read outcome data
outcome_dat <- read_outcome_data(
    file.path(snakemake@input[["mr_data"]]),
    sep = " ",
    snp_col = "SNP", 
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "Freq_OUT",
    beta_col = "b_OUT",
    se_col = "se_OUT",
    pval_col = "p_OUT",
    samplesize_col = "N_OUT"
)
outcome_dat$outcome <- OUTCOME
    
### Harmonize data (redone to have it in the package format)
dat <- harmonise_data(exposure_dat, outcome_dat)
# In cases where eaf = 0 NA's are introduced
if (sum(rowSums(is.na(dat)) > 0) != 0) {
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
dat$zval_steiger <- (abs(dat$beta.exposure)-abs(dat$beta.outcome)) / sqrt(dat$se.exposure^2 + dat$se.outcome^2)
dat <- dat[dat$zval_steiger > -1.96, ]
print(paste0("Number of SNPs after steiger filter: ", nrow(dat)))
        
### Remove SNP from the HBB locus (chr11 5246696-5248301)
dat <- dplyr::filter(dat, !SNP %in% HBB_rsids$SNP)
print(paste0("Number of SNPs without HBB locus: ", nrow(dat)))
    
### Remove SNPs from HLA locus (chr6 25MB-37MB)
dat <- dplyr::filter(dat, !SNP %in% HLA_rsids$SNP)
print(paste0("Number of SNPs without HLA locus: ", nrow(dat)))
    
### MR analysis
if (nrow(dat) > 3) {
    ############################# RUN MR-PRESSO ##############################
    # Initialize variables
    nb_distribution <- round_up_1000(nrow(dat))  # starting value (at least upper 1000 of number of SNPs)
    max_nb_distribution <- 21000  # maximum allowed value
    warning_occurred <- TRUE  # to enter the while loop at least once
    nb_increase_dist <- 5000
      
    handle_warning <- function(w) {
        # Check if the warning message contains the expected string
        if (grepl("Outlier test unstable", conditionMessage(w))) {
            warning_occurred <<- TRUE  # set warning flag if a warning occurs
            nb_distribution <<- nb_distribution + nb_increase_dist  # increase nb_distribution
            if (nb_distribution <= max_nb_distribution) {
                message("Warning occurred: ", conditionMessage(w), 
                        "\nRetrying with NbDistribution = ", nb_distribution, ".")
            } else {
                message("Maximum NbDistribution reached without resolving the issue. Proceeding with NbDistribution = ", nb_distribution, ".")
            }
        } else {
            # Handle other warnings if needed, or simply print them
            message("Warning occurred: ", conditionMessage(w))
        }
    }
      
    # Loop to adjust nb_distribution and re-run MR-PRESSO
    while (warning_occurred & nb_distribution <= max_nb_distribution) {
        warning_occurred <- FALSE  # reset warning flag
        
        # Try running MR-PRESSO and catch any warnings
        tryCatch({
            mr_results <- TwoSampleMR::run_mr_presso(
                dat, 
                NbDistribution = nb_distribution
            )
        },
        warning = handle_warning,
        error = function(e) {
            stop("Error occurred: ", conditionMessage(e))  # stop on errors
        })
    }
      
    ##########################################################################
     
    # Extracting main MR results
    main_results <- mr_results[[1]]$`Main MR results`
    main_results_df <- as.data.frame(main_results)
      
    # Extracting MR-PRESSO results
    global_test_pvalue <- as.numeric(sub("<", "", mr_results[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue))
    global_test_rssobs <- mr_results[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs
    distortion_test_pvalue <- mr_results[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue
    distortion_coeff <- mr_results[[1]]$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`
    outliers <- length(mr_results[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
      
    # Adding MR-PRESSO results to the main results dataframe
    main_results_df$global_test_pvalue <- global_test_pvalue
    main_results_df$global_test_rssobs <- global_test_rssobs
    main_results_df$distortion_test_pvalue <- distortion_test_pvalue
    main_results_df$distortion_coeff <- distortion_coeff
    main_results_df$nb_distribution_set <- nb_distribution

    # Rename column
    main_results_df <- main_results_df %>% rename(exposure = "Exposure") %>%
        rename(method = "MR Analysis") %>% rename(b = "Causal Estimate") %>%
        rename(sd = "Sd") %>% rename(pval = "P-value") %>% rename(tstat = "T-stat")
      
    # Add number of snps (remove outliers)
    main_results_df$nsnp <- ifelse(main_results_df$method == "Raw", nrow(dat), nrow(dat) - outliers)
      
    # Add se
    main_results_df$se <- abs(main_results_df$b) / main_results_df$tstat
      
    # Add exposure and outcome
    main_results_df$exposure <- EXPOSURE
    main_results_df$outcome <- OUTCOME
      
    # Reordering columns
    main_results_df <- main_results_df %>% 
        select(exposure, outcome, method, nsnp, b, se, pval, everything())
      
    # Displaying the comprehensive dataframe
    print(main_results_df)
      
} else {
    methods_list <- c("Raw","Outlier-corrected")
    main_results_df <- data.frame(
        exposure = rep(EXPOSURE, length(methods_list)),
        outcome = rep(OUTCOME, length(methods_list)),
        method = methods_list,
        nsnp = rep(nrow(dat), length(methods_list)),
        b = rep(NA, length(methods_list)),
        se = rep(NA, length(methods_list)),
        pval = rep(NA, length(methods_list)),
        sd = rep(NA, length(methods_list)),
        tstat = rep(NA, length(methods_list)),
        global_test_pvalue = rep(NA, length(methods_list)),
        global_test_rssobs = rep(NA, length(methods_list)),
        distortion_test_pvalue = rep(NA, length(methods_list)),
        distortion_coeff = rep(NA, length(methods_list)),
        nb_distribution_set = rep(NA, length(methods_list))
    )
}

write.table(main_results_df, snakemake@output[[1]], row.names = F, quote = F, sep = '\t')
