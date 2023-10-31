################################################################################
### MR Plots from MR prepared data                                           ###
### Author: Samuel Moix                                                      ###
### Date: 11.07.2022                                                         ###
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
### Working directories ########################
data_folder <- "/SET/PATH/TO/DIRECTORY/"
export_folder <- "/SET/PATH/TO/DIRECTORY/"
tsv_file_name <- "_MR_data_sens.tsv"

################################################
### Parameters  ################################

### Whether to plot (existing plots won't be re-plotted)
par_plot <- TRUE
### Whether to recompute MR analysis
par_keep_old_results <- TRUE
### Whether to compute the reverse MR effects
par_reverse <- TRUE # traits --> st_TRAIT if TRUE
#                    # st_TRAIT --> traits if FALSE

################################################
### Set trait of interest ######################
st_TRAIT <- "TELOMERE"


################################################
### Rsids to exclude files #####################

HLA_rsids <- fread("HLA_rsids.txt", header = F)
HBB_rsids <- fread("HBB_rsids.txt", header = F)

################################################
### Load regression/already computed MR results 

### File containing FiledID, pheno and added .ma format phenotype names
lm_res <- read.csv(file = "TL_metadata.csv", sep = ',', header = TRUE)

# Load previous results to avoid recomputing them
### Load already computed MR results 
file_name_mr <- paste0(st_TRAIT, "_on_TRAIT_mr_results.txt")

if(par_reverse){
  file_name_mr <- paste0("TRAIT_on_",st_TRAIT, "_mr_results.txt")
}
if(par_keep_old_results){
  if(file.exists(file.path(export_folder,  file_name_mr))){
    mr_res_table <- read.table(file = file.path(export_folder,  file_name_mr), sep = '\t', header = TRUE)
  }
}

################################################
### TwoSampleMR ################################

### Studied traits causal effect on all traits
# Set bolean for first passage (result table)
first_trait <- TRUE

# Get list of traits to analyse and remove
list_of_traits <- lm_res$MR_name[which(lm_res$MR_name != "")]
list_of_traits <- list_of_traits[which(list_of_traits != st_TRAIT)]

if(par_keep_old_results){
  if(file.exists(file.path(export_folder,  file_name_mr))){
    if(!par_reverse){
      list_of_traits <- list_of_traits[which(!list_of_traits %in% mr_res_table$outcome)]
    }
    if(par_reverse){
      list_of_traits <- list_of_traits[which(!list_of_traits %in% mr_res_table$exposure)]
    }
    first_trait <- FALSE
  }
}

# Create vector to report non existing files
missing_traits <- c()
double_files <- c()

for (TRAIT in list_of_traits) {
  ################################################
  ### Selecting exposure and outcome #############
  EXPOSURE <- st_TRAIT
  OUTCOME <- TRAIT
  
  if(par_reverse){
    EXPOSURE <- TRAIT
    OUTCOME <- st_TRAIT
  }
  
  
  file_name <- paste0(EXPOSURE, "_to_", OUTCOME, tsv_file_name) 
  ### Adapt file name if in sub-directory
  file_name <- list.files(path = data_folder, pattern = file_name, recursive = TRUE)
  # Write if files exists more than twice
  if(length(file_name) > 1){
    double_files <- c(double_files, file_name[1])
  }
  
  ### Check weather file exists (was listed)
  if(!identical(file_name, character(0))){
    
    ### Read exposure data
    exposure_dat <- read_exposure_data(file.path(data_folder,  file_name[1]),
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
    outcome_dat <- read_outcome_data(file.path(data_folder,  file_name[1]),
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
    dat <- dplyr::filter(dat, !SNP %in% HBB_rsids$V1)
    print(paste0("Number of SNPs without HBB locus: ", nrow(dat)))
    
    ### Remove SNPs from HLA locus (chr6 25MB-37MB)
    dat <- dplyr::filter(dat, !SNP %in% HLA_rsids$V1)
    print(paste0("Number of SNPs without HLA locus: ", nrow(dat)))
    
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
    
    if(!first_trait){
      mr_res_table <- bind_rows(mr_res_table, mr_add_res_table)
    }
    
    if(first_trait){
      mr_res_table <- mr_add_res_table
      first_trait <- FALSE
    }
    
    ### Plotting
    if(par_plot){
      ### Plot export file name
      plot_file_name <- paste0(EXPOSURE, "_to_", OUTCOME, "_MR_plots.pdf")
      ### Avoid regenerating existing plots
      if(!file.exists(file.path(export_folder, plot_file_name))){
        # Scatter plot
        p1 <- mr_scatter_plot(mr_results, dat)
        
        # Forest plot
        res_single <- mr_singlesnp(dat)
        p2 <- mr_forest_plot(res_single)
        
        # Leave-one-out plot
        res_loo <- mr_leaveoneout(dat)
        p3 <- mr_leaveoneout_plot(res_loo)
        
        # Funnel plot
        p4 <- mr_funnel_plot(res_single)
        
        ### Plot all plots together
        mr_plots <- gridExtra::grid.arrange(grobs = c(p1,p2,p3,p4), ncol=2)
        
        export_file_name <- paste0(plot_file_name)
        ggsave(file.path(export_folder, export_file_name),
               mr_plots, width=3, height=3, units="in", scale=3)
      }
    }
  }
  else{
    missing_traits <- c(missing_traits,TRAIT)
  }
}

print("These traits are missing")
print(missing_traits)

### Save results
fwrite(mr_res_table, file.path(export_folder,  file_name_mr), col.names = T, row.names = F, quote = F, sep = "\t")

### Save info
fileConn<- file(file.path(export_folder,
                          paste0(strsplit(file_name_mr, split = "mr_results.txt"),"missing.txt")), "a")
write(paste(c(missing_traits, double_files)), file = fileConn, append = TRUE)
close(fileConn)

