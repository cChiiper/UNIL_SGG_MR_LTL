################################################################################
### Plot the different methods estimates                                     ###
### Author: Samuel Moix                                                      ###
### Date: 14.09.2023                                                         ###
################################################################################

################################################
### Libraries ##################################
library(dplyr)
library(ggplot2)
library(ggbreak)
library(conflicted)

library(ggforestplot)

################################################
### Working directories ########################

### Input files
metadata_file <- "/SET/PATH/TO/DIRECTORY/TL_metadata.csv"
MR_res_forward_file <- "/SET/PATH/TO/DIRECTORY/TELOMERE_on_TRAIT_mr_results.txt"
MR_res_reverse_file <-"/SET/PATH/TO/DIRECTORY/TRAIT_on_TELOMERE_mr_results.txt"
MR_PRESSO_file <- "/SET/PATH/TO/DIRECTORY/MR_PRESSO_results.tsv"
indep_trait <- 141

################################################
### Load data ##################################

### Load metadata
metadata <- read.csv(metadata_file, sep = ',', header = TRUE)

### Load bi-directional MR files
MR_F <- read.table(MR_res_forward_file, header = T, sep = "\t")
MR_R <- read.table(MR_res_reverse_file, header = T, sep = "\t")

### Add MR-PRESSO
MR_P <- read.table(MR_PRESSO_file, header = T, sep = "\t")

# Forward
MR_P_F <- MR_P %>% 
  select("exposure","outcome","method","nsnp","b","se","pval") %>%
  dplyr::filter(method == "Outlier-corrected" & exposure == "TELOMERE")
MR_P_F$method <- "MR-PRESSO"
MR_P_F$Q <- NA
MR_P_F$Q_pval <- NA
MR_F <- rbind(MR_F, MR_P_F)
rm(MR_P_F)

# Reverse
MR_P_R <- MR_P %>% 
  select("exposure","outcome","method","nsnp","b","se","pval") %>%
  dplyr::filter(method == "Outlier-corrected" & outcome == "TELOMERE")
MR_P_R$method <- "MR-PRESSO"
MR_P_R$Q <- NA
MR_P_R$Q_pval <- NA
MR_R <- rbind(MR_R, MR_P_R)
rm(MR_P_R)


################################################
### Main script forward ########################

### Select traits of interest
selected_outcomes <- MR_F %>%
  dplyr::filter(method == "Inverse variance weighted" & pval < (0.05/indep_trait)) %>%
  arrange(b) %>%
  pull(outcome) 

selected_outcomes <- selected_outcomes[!selected_outcomes %in% c("ENDOMET_UKBB")]

MR_F <- MR_F %>%
  dplyr::filter(outcome %in% selected_outcomes)

### Add trait description
MR_F <- dplyr::rename(MR_F, MR_name = outcome)
MR_F <- merge(MR_F, metadata[,c("MR_name","description","Category")], by = "MR_name")

MR_F <- MR_F %>%
  mutate(MR_name = factor(MR_name, levels = selected_outcomes)) %>%
  arrange(MR_name, method)

# Compute the confidence intervals
MR_F$lower_ci <- MR_F$b - 1.96 * MR_F$se
MR_F$upper_ci <- MR_F$b + 1.96 * MR_F$se

### Plotting

plot <- ggforestplot::forestplot(
  df = MR_F,
  name = description,
  estimate = b,
  pvalue = pval,
  psignif = 0.05/indep_trait,
  xlab = expression(""*alpha*" (95% CI)"), #"\u03b1(CI95%)"
  title = "",
  colour = method,
  se = se
) + 
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11)) # Text size


p <- plot  +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(
    breaks = c(-0.9, -0.4, -0.2, 0, 0.2, 0.4, 0.9),
    limits = c(round(min(MR_F$lower_ci),1) - 0.1, round(max(MR_F$upper_ci),1) + 0.1),
    position = "bottom"
  ) +
  ggbreak::scale_x_break((c(round(min(MR_F$lower_ci),1) + 0.01, -0.5))) + 
  ggbreak::scale_x_break(c(0.53, round(max(MR_F$upper_ci),1) + 0.01)) + 
    theme(legend.title=element_blank(),
          legend.position = c("top"),
          legend.margin = margin(1, 150, 1, 1),
          axis.title.x.top = element_blank(),
          axis.text.x.top = element_blank(),
          axis.title.x = element_text(hjust = 0.65, vjust = 1))

ggsave(file="/SET/PATH/TO/DIRECTORY/plot1.png", 
       plot = p, width = 8, height = 9.5, dpi = 400)
ggsave(file="/SET/PATH/TO/DIRECTORY/plot1.svg", 
       plot = p, width = 8, height = 9.5)

################################################
### Main script reverse ########################

### Select traits of interest
selected_outcomes <- MR_R %>%
  dplyr::filter(method == "Inverse variance weighted" & pval < (0.05/indep_trait)) %>%
  arrange(b) %>%
  pull(exposure) 

MR_R <- MR_R %>%
  dplyr::filter(exposure %in% selected_outcomes)

### Add trait description
MR_R <- dplyr::rename(MR_R, MR_name = exposure)
MR_R <- merge(MR_R, metadata[,c("MR_name","description","Category")], by = "MR_name")

MR_R <- MR_R %>%
  mutate(MR_name = factor(MR_name, levels = selected_outcomes)) %>%
  arrange(MR_name, method)

# Compute the confidence intervals
MR_R$lower_ci <- MR_R$b - 1.96 * MR_R$se
MR_R$upper_ci <- MR_R$b + 1.96 * MR_R$se

### Plotting

plot <- ggforestplot::forestplot(
  df = MR_R,
  name = description,
  estimate = b,
  pvalue = pval,
  psignif = 0.05/indep_trait,
  xlab = expression(""*alpha*" (95% CI)"), #"\u03b1(CI95%)"
  title = "",
  colour = method,
  se = se
) + 
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)) # Text size

p <- plot  +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(
    breaks = c(-0.5,-0.2,0, 0.2, 0.4,0.9, 1.1, 1.6),
    limits = c(round(min(MR_R$lower_ci),1) - 0.05, round(max(MR_R$upper_ci),1) + 0.1),
    position = "bottom"
  ) +
  ggbreak::scale_x_break(c(0.5,0.82)) + 
  ggbreak::scale_x_break(c(1.2, 1.6)) + 
  theme(legend.title=element_blank(),
        legend.position = c("top"),
        legend.margin = margin(1, 170, 1, 1),
        axis.title.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.title.x = element_text(hjust = 0.65, vjust = 1))

ggsave(file="/SET/PATH/TO/DIRECTORY/plot2.png", 
       plot = p, width = 8, height = 9.5, dpi = 500)
ggsave(file="/SET/PATH/TO/DIRECTORY/plot2.svg", 
       plot = p, width = 8, height = 9.5)
