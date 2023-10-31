################################################################################
### Comparing LTL corrected to LTL corrected also for blood counts           ###
### Author: Samuel Moix                                                      ###
### Date: 28.08.2023                                                         ###
################################################################################

################################################
### Libraries ##################################

library(data.table) # Read data 
library(ggplot2) # Plotting basics
library(dplyr)
library(forcats)
require(ggrepel)
require(tibble)
library(readr)


################################################
### Working directories ########################
data_folder <- "/SET/PATH/TO/DIRECTORY/"
export_folder <- "/SET/PATH/TO/DIRECTORY/"
metadata_file <- "/SET/PATH/TO/DIRECTORY/TL_metadata.csv"

### Parameters
indep_trait <- 141

################################################
### Load files #################################
metadata <- read.csv(metadata_file, sep = ',', header = TRUE)
lm_results <- read_delim(file.path(data_folder,"lm_res_scaled_no_outliers_no_bc.tsv"), 
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)
lm_results_blood <- read_delim(file.path(data_folder,"lm_res_scaled_no_outliers_no_bc_corbc.tsv"), 
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)


### Merge data
# Drop the 'mean' and 'sd' columns
lm_results <- lm_results[, !(names(lm_results) %in% c("mean", "sd"))]
lm_results_blood <- lm_results_blood[, !(names(lm_results_blood) %in% c("mean", "sd"))]

# Rename the specified columns in lm_results_blood
cols_to_rename <- c("rcoef", "CI1", "CI2", "p_value", "se")
new_names <- paste0(cols_to_rename, "_blood")
names(lm_results_blood)[names(lm_results_blood) %in% cols_to_rename] <- new_names

# Merge the two dataframes on 'pheno'
merged_df <- merge(lm_results, lm_results_blood, by = "pheno", all = TRUE)

# Drop the redundant 'description' and 'Category' columns from the merged dataframe
merged_df$description.y <- NULL
merged_df$Category.y <- NULL

# Rename the remaining 'description.x' and 'Category.x' columns to 'description' and 'Category'
names(merged_df)[names(merged_df) == "description.x"] <- "description"
names(merged_df)[names(merged_df) == "Category.x"] <- "Category"

# Clean
rm(cols_to_rename, new_names)

################################################
### Add pdiff and plot #########################

### Add pdiff
merged_df <- mutate(merged_df, 
                     pdiff = 2*pnorm(-abs((rcoef-rcoef_blood)/sqrt(se**2+se_blood**2)), mean = 0, sd = 1))

# Remove non-wanted labels
merged_df$description[which(merged_df$pdiff > 0.05)] <- ""

blood_count_traits <- c("lymphocyte_count", "WBC_count", "RBC_count", "platelet_count",
                        "monocyte_count", "eosinophil_count", "neutrophil_count",
                        "reticulocyte_count")

merged_df <- merged_df %>%
  dplyr::filter(!pheno %in% blood_count_traits)

### Plot
# Fit a linear model
model <- lm(rcoef_blood ~ rcoef, data = merged_df)
# Extract coefficients
intercept <- coef(model)[1]
slope <- coef(model)[2]
equation <- paste0("y = ", round(slope, 3), "x")

set.seed(1234)
comp_plot <- ggplot(merged_df, aes(x=rcoef, y=rcoef_blood)) + 
  # geom_abline(intercept = intercept, slope = slope, color = "#FF43FF", alpha = 0.8, linetype = "dashed") +
  # annotate("text", x = -0.1, y = 0.1, label = equation, hjust = "left", size = 5, color = "#FF43FF") + 
  geom_errorbar(aes(ymax = CI2_blood, ymin = CI1_blood, width = 0), color = "#c2c2c2",
                alpha = ifelse(merged_df$pdiff < 0.05,0.8,0)) +
  geom_errorbar(aes(xmax = CI2, xmin = CI1, width = 0), color = "#c2c2c2", 
                alpha = ifelse(merged_df$pdiff < 0.05,0.8,0)) +
  geom_point(pch=21, 
             fill =  ifelse(merged_df$pdiff < 0.05/indep_trait,"#FF0000","#FFC0C0"), 
             size =  ifelse(merged_df$pdiff < 0.05, 3.5, 2), 
             colour = "black",
             alpha = ifelse(merged_df$pdiff < 0.05,0.9,0.2)) + # Show dots
  geom_abline(slope = 1, color = "#685044", alpha = 0.4, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  ggrepel::geom_text_repel(
    label=merged_df$description, 
    segment.color = "#6c757d", 
    segment.alpha = 0.6,
    segment.linetype = 3,
    size = 4.5,
    #fontface = "bold"
    color = ifelse(merged_df$pdiff < 0.05/indep_trait, "#FF0000", "black"),
    alpha =  ifelse(merged_df$pdiff < 0.05/indep_trait, 1, 0.8),
    max.overlaps = 200,
    box.padding = 2.9,
    seed = 1996
  )+
  xlab(expression(""*beta*" (95% CI)")) + 
  ylab(expression(""*beta*" corrected for blood counts (95% CI)"))


mylim <- 0.1
blood_plot <- comp_plot + 
  theme_light() +
  scale_x_continuous(breaks = seq(-0.09, 0.09, by = 0.02), limits = c(-mylim, mylim)) +
  scale_y_continuous(breaks = seq(-0.09, 0.09, by = 0.02), limits = c(-mylim, mylim)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title=element_text(size=25),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 


ggsave(file="/SET/PATH/TO/DIRECTORY/plot1.svg", 
       plot = blood_plot, width=12, height=12)
ggsave(file="/SET/PATH/TO/DIRECTORY/plot1.png", 
       plot = blood_plot, width = 12, height = 12, dpi = 300)

