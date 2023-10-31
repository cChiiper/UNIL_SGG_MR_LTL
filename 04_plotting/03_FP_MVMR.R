################################################################################
### Plotting Forrest-plot of MVMR analysis (multiple exp/mediators)          ###
### Author: Samuel Moix                                                      ###
### Date: 15.03.2023                                                         ###
################################################################################

################################################
### Libraries ##################################
library(ggplot2)
library(dplyr)
require(ggforestplot)

################################################
### Working directories ########################
data_folder = "/SET/PATH/TO/DIRECTORY/"

################################################
### Load data ##################################

### Load metadata
metadata <- read.csv(file = "/SET/PATH/TO/DIRECTORY/TL_metadata.csv", sep = ',', header = TRUE)

### Load results
# FILES: HI_to_TELOMERE_MVMR_result.tsv and ACLT_to_TELOMERE_MVMR_result.tsv
plot_df <- read.table(file.path(data_folder,"ACLT_to_TELOMERE_MVMR_result.tsv"), header = T, sep = "\t")

################################################
### Prepare plotting dataframe #################

### Add description
plot_df$MR_name <- plot_df$exposure
plot_df <- merge(plot_df, metadata[,c("MR_name","description")], by="MR_name")
outcome_name <- metadata$description[which(metadata$MR_name == unique(plot_df$outcome))]
rm(metadata)

### Rename methods
plot_df$relation[which(plot_df$relation == "ivw_steiger")] <- "IVW MR (stringent steiger)"
plot_df$relation[which(plot_df$relation == "ivw")] <- "IVW MR"

################################################
### Plotting ###################################

ggforestplot::forestplot(
  df = plot_df,
  name = description,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = expression(""*alpha*" (CI95%)"),
  title = paste("Outcome:", outcome_name),
  colour = relation,
  se = se) +   
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(axis.text.y = element_text(size = 10),
        legend.position="top",
        legend.title=element_blank())


### Plotting Lipids
p <- forestplot(
  df = plot_df,
  name = description,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = expression(""*alpha*" (95% CI)"),
  title = "",
  colour = relation,
  se = se) +   
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values = c("#6699CC","#0C5197","#568259")) + # Use custom colors
  theme(axis.text.y = element_text(size = 12),
        legend.position="top",
        legend.title=element_blank())

ggsave(file="/SET/PATH/TO/DIRECTORY/plot1.png", 
       plot = p, width = 6.93, height = 3.35, dpi = 400)
ggsave(file="/SET/PATH/TO/DIRECTORY/plot1.svg", 
       plot = p, width = 6.93, height = 4)

ggsave(file="/SET/PATH/TO/DIRECTORY/plot2.png", 
       plot = p, width = 7, height = 5, dpi = 400)
