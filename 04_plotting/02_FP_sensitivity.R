################################################################################
### Plotting alpha and beta estimates of TL relationship with complex traits ###
### Author: Samuel Moix                                                      ###
### Date: 17.03.2023                                                         ###
################################################################################

################################################
### Libraries ##################################
library(ggplot2)
library(dplyr)
library(cowplot)
library(conflicted)
library(gridExtra)

################################################
### Input files ################################
metadata_file <- "/SET/PATH/TO/DIRECTORY/TL_metadata.csv"
MR_res_reverse_sens_file <- "/SET/PATH/TO/DIRECTORY/TELOMERE_on_TRAIT_mr_results.txt"
# MR_res_reverse_sens_file <- "/SET/PATH/TO/DIRECTORY/TRAIT_on_TELOMERE_mr_results.txt"
MR_res_reverse_file <-"/SET/PATH/TO/DIRECTORY/TELOMERE_on_TRAIT_mr_results.txt"
# MR_res_reverse_file <- "/SET/PATH/TO/DIRECTORY/TRAIT_on_TELOMERE_mr_results.txt"
#
###############################################
### Parameters #################################
indep_trait <- 141

################################################
### Read files #################################
mr_res_table <- read.table(file = MR_res_reverse_file, sep = '\t', header = TRUE)
mr_res_sens_table <- read.table(file = MR_res_reverse_sens_file, sep = '\t', header = TRUE)
metadata <- read.csv(metadata_file, sep = ',', header = TRUE)

################################################
### Merge and filter dataframes ################
mr_res_table$analysis <- "all_snps"
mr_res_sens_table$analysis <- "sens_snps"
mr_res_all <- rbind(mr_res_table, mr_res_sens_table)

mr_res_all <- mr_res_all %>%
  dplyr::filter(method == "Inverse variance weighted") %>%
  rename(MR_name = outcome) # !!!!!!!!!!!!!!!!!!!!!! Has to changed all outcomes to exposure

mr_res_all <- merge(mr_res_all, metadata[,c("MR_name","description","Category")], by = "MR_name")

# Keep only traits with SENS data
traits_SENS <- mr_res_all$description[which(mr_res_all$analysis == "sens_snps")]
mr_res_all <- mr_res_all %>% 
  dplyr::filter(description %in% traits_SENS)

df_plot <- mr_res_all %>%
  group_by(MR_name) %>%
  dplyr::filter(any(pval < 0.05/indep_trait)) %>%
  ungroup() %>%
  arrange(ifelse(analysis == "all_snps", b, Inf))

require(tidyr)
df_integers <- df_plot %>%
  select(description, analysis, nsnp) %>%
  tidyr::spread(analysis, nsnp) %>%
  mutate(description_with_snps = paste(all_snps, sens_snps, description, sep = " | ")) %>%
  select(description, description_with_snps)

df_plot <- df_plot %>%
  left_join(df_integers, by = "description")

################################################
### Plot #######################################
p1 <- ggforestplot::forestplot(
  df = df_plot,
  name = description_with_snps,
  estimate = b,
  pvalue = pval,
  psignif = 0.05/indep_trait,
  xlab = expression(""*alpha*" (95% CI)"),
  title = "",
  colour = analysis,
  se = se
) + 
  scale_color_manual(
    values = c("sens_snps" = "#e69988", "all_snps" = "#CC3311"), #b2cce6 6699CC
    labels = c("sens_snps" = "Sensitvity SNPs", "all_snps" = "All SNPs")
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 3))  # Increase size of the dots in the legend
  ) +
  labs(colour = NULL) +  # Remove the legend title
  theme(
    axis.text.y = element_text(size = 10),
    legend.title=element_blank(),
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(1, 1, 1, 1),
    legend.background = element_rect(fill="#F1EDE6",
                                     size=0.5, linetype="solid", 
                                     colour ="black")
  )

print(p1)

# https://bggj.is/posts/forest-plot-table/


################################################
### Merge and filter dataframes for scatter ####

merged_df <- merge(mr_res_table[,c("exposure","outcome","method","b","se","pval")], 
                   mr_res_sens_table[,c("outcome","method","b","se","pval")], 
                   by=c("outcome","method"), suffixes = c("_all", "_sens"))

# Filter traits that are a least significant for one method
merged_df <- merged_df %>% 
  dplyr::filter(pval_all < 0.05/indep_trait | pval_sens < 0.05)

# Assuming your merged dataframe is named 'merged_df'
# Add a column indicating whether values are significant (pval <= 0.05)
merged_df$sig <- (merged_df$pval_all <= 0.05) & (merged_df$pval_sens <= 0.05)

# Calculate the 95% confidence intervals
merged_df$lower_ci_all <- merged_df$b_all - 1.96 * merged_df$se_all
merged_df$upper_ci_all <- merged_df$b_all + 1.96 * merged_df$se_all
merged_df$lower_ci_sens <- merged_df$b_sens - 1.96 * merged_df$se_sens
merged_df$upper_ci_sens <- merged_df$b_sens + 1.96 * merged_df$se_sens

# Add difference significance
merged_df <- dplyr::mutate(merged_df, 
                     pdiff = 2*pnorm(-abs((b_all-b_sens)/sqrt(se_all**2+se_sens**2)), mean = 0, sd = 1))

# Subset to Inverse variance weighted
merged_df <- dplyr::filter(merged_df, method == "Inverse variance weighted")

# Add metadata
merged_df <- dplyr::rename(merged_df, MR_name = outcome)
merged_df <- merge(merged_df, metadata[,c("MR_name","description","Category")], by = "MR_name")

# Remove non-wanted labels
merged_df$description[which(merged_df$pdiff > 0.05)] <- ""

# Compute correlation
correlation <- cor(merged_df$b_all, merged_df$b_sens)


################################################
### Plot #######################################

comp_plot <- ggplot(merged_df, aes(x=b_all, y=b_sens)) +
  geom_errorbar(aes(ymax = upper_ci_sens, ymin = lower_ci_sens, width = 0), color = "#c2c2c2",
                alpha = ifelse(merged_df$pdiff < 0.05,0.8,0)) +
  geom_errorbar(aes(xmax = upper_ci_all, xmin = lower_ci_all, width = 0), color = "#c2c2c2", 
                alpha = ifelse(merged_df$pdiff < 0.05,0.8,0)) +
  geom_point(pch=21, 
             fill =  ifelse(merged_df$pdiff < 0.05/indep_trait,"#8338ec","#CBBDDD"), 
             size = 3, 
             colour = "black",
             alpha = ifelse(merged_df$pdiff < 0.05,0.9,0.2)) + # Show dots
  geom_abline(slope = 1, color = "#685044", alpha = 0.4, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  ggrepel::geom_text_repel(
    label=merged_df$description, 
    size = 4,
    #fontface = "bold"
    color = ifelse(merged_df$pdiff < 0.05/indep_trait, "#E63946", "black"),
    max.overlaps = 250,
    box.padding = 2.5,
    seed = 1996,
    segment.color = "grey",
    segment.alpha = 0.9,
    segment.linetype = "dotted"
  )+
  xlab(expression(""*alpha*" all SNPs (95% CI)")) + 
  ylab(expression(""*alpha*" sensitivtiy SNPs (95% CI)"))

mylim <- 0.3
p2 <- comp_plot + 
  theme_light() +
  scale_x_continuous(breaks = seq(-0.4, 0.4, by = 0.1), limits = c(-mylim, mylim)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, by = 0.1), limits = c(-mylim, mylim)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title=element_text(size=15),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x = -mylim + 0.1, y = mylim - 0.01, 
           label = paste("rho =", round(correlation, 2)), size = 4)

p2 <- p2 + coord_fixed(ratio = 1)
print(p2)

### Plot together:
p3 <- grid.arrange(p1, p2, ncol=2)


ggsave(file="/SET/PATH/TO/DIRECTORY/sensitvtiy_steiger.png", 
       plot = p3, width = 14, height = 6, dpi = 400)

ggsave(file="/SET/PATH/TO/DIRECTORY/forward_sensitvtiy_steiger.svg", 
       plot = p3, width = 14, height = 6)
