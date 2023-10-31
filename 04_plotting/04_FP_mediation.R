################################################################################
### Plotting all pairs of mediated traits by telomere length                 ###
### Author: Samuel Moix                                                      ###
### Date: 18.08.2023                                                         ###
################################################################################


################################################
### Libraries ##################################
library(ggplot2)
library(dplyr)
library(ggforestplot)
library(conflicted)

################################################
### Working directories ########################

med_file <- "/SET/PATH/TO/DIRECTORY/results.tsv"
par_two_methods <- TRUE
multiple_testing <- 359

################################################
### Load data ##################################

### Load metadata
metadata <- read.csv(file = "/SET/PATH/TO/DIRECTORY/TL_metadata.csv", sep = ',', header = TRUE)
metadata <- metadata[,c("MR_name","description")]
### Load mediation results
med_data <- read.table(med_file, sep ="\t", header = T)

################################################
### Filtering ##################################

### Filter methods
if(!par_two_methods){
  med_data <- med_data %>%
    dplyr::filter(Relation %in% c("Indirect effect - Shrinkage - B", "Direct effect - Shrinkage",
                                  "Total effect"))
  med_data$Relation[which(med_data$Relation == "Indirect effect - Shrinkage - B")] <- "Indirect effect"
  med_data$Relation[which(med_data$Relation == "Direct effect - Shrinkage")] <- "Direct effect"
}

if(par_two_methods){
  med_data <- med_data %>%
    dplyr::filter(Relation %in% c("Indirect effect - Shrinkage - B", "Direct effect - Shrinkage",
                                  "Total effect", "Indirect effect - Shrinkage - A"))
  med_data$Relation[which(med_data$Relation == "Indirect effect - Shrinkage - A")] <- "Indirect effect difference"
  med_data$Relation[which(med_data$Relation == "Indirect effect - Shrinkage - B")] <- "Indirect effect product"
  med_data$Relation[which(med_data$Relation == "Direct effect - Shrinkage")] <- "Direct effect"
}


### Remove trait to same trait
med_data <- med_data %>%
  dplyr::filter(exposure != outcome)

### Add description
# Make labels shorter
metadata[which(metadata$MR_name == "EDUAGE"),c("description")] <- "EA"
metadata[which(metadata$MR_name == "BMI"),c("description")] <- "BMI"
metadata[which(metadata$MR_name == "MCH"),c("description")] <- "MCH"


# Joining for exposure
med_data <- med_data %>%
  left_join(metadata %>% select(MR_name, description), by = c("exposure" = "MR_name")) %>%
  rename(exp_description = description)

# Joining for outcome
med_data <- med_data %>%
  left_join(metadata %>% select(MR_name, description), by = c("outcome" = "MR_name")) %>%
  rename(out_description = description)


### Add column for name
plot_df <- med_data %>%
  mutate(tag = paste(exposure, outcome, sep="_to_")) %>%
  mutate(lable = paste(exp_description, out_description, sep = " to "))

### Remove blood traits
blood_traits_MR <- c("MCH","EOSINO","LYMPH","MONO","NEUTRO","PLATELET","RBC","RETICULO","WBC")
plot_df <- plot_df %>%
  dplyr::filter(!exposure %in% blood_traits_MR) %>%
  dplyr::filter(!outcome %in% blood_traits_MR) %>%
  dplyr::filter(!exposure %in% c("TRI")) %>%
  dplyr::filter(!outcome %in% c("MS")) %>%
  dplyr::filter(!outcome %in% c("ENDOMET_UKBB")) 


### Filter by p-val
plot_df <- plot_df %>%
  dplyr::group_by(tag) %>%
  dplyr::filter(any(pval < 0.05/multiple_testing & Relation == "Total effect")) %>%
  dplyr::ungroup()

if(!par_two_methods){
  plot_df <- plot_df %>%
    dplyr::group_by(tag) %>%
    dplyr::filter(any(pval < 0.05/multiple_testing & Relation == "Indirect effect")) %>%
    dplyr::ungroup()
}

if(par_two_methods){
  plot_df <- plot_df %>%
    dplyr::group_by(tag) %>%
    dplyr::filter(any(pval < 0.05/multiple_testing & Relation == "Indirect effect product")) %>%
    dplyr::ungroup()
}




### Reorder by Total effect beta
order_df <- plot_df %>%
  dplyr::filter(Relation == "Total effect") %>%
  arrange(beta)

plot_df <- plot_df %>%
  group_by(tag) %>%  
  arrange(factor(tag, levels = order_df$tag)) %>%
  ungroup()


################################################
### Plotting ###################################

### Forest plot
p <- ggforestplot::forestplot(
  df = plot_df,
  name = lable,
  estimate = beta,
  pvalue = pval,
  psignif = 0.05/116,
  xlab = expression(""*alpha*" (CI95%)"),
  title = "",
  colour = Relation,
  se = se) 

if(!par_two_methods){
  p <- p + theme(axis.text.y = element_text(size = 13, hjust = 1),
            legend.position = "bottom", # Move legend to the bottom
            axis.text.x = element_text(size = 13),
            legend.text = element_text(size = 12),
            legend.title = element_blank(),
            axis.title.x = element_text(size = 15)) + 
    scale_color_manual(values=c("#FE9000","#DC63CB","#000000"),
                       breaks=c('Indirect effect', 'Direct effect','Total effect')) +
    guides(colour = guide_legend(override.aes = list(size=2)))
}


if(par_two_methods){
  p <- p + theme(axis.text.y = element_text(size = 13, hjust = 1),
             legend.position = "bottom", # Move legend to the bottom
             axis.text.x = element_text(size = 13),
             legend.text = element_text(size = 12),
             legend.title = element_blank(),
             axis.title.x = element_text(size = 15)) + 
    scale_color_manual(values=c("#000000","#DC63CB","#007BFF","#FE9000"),
                       breaks=c('Total effect', 'Direct effect','Indirect effect difference','Indirect effect product')) +
    guides(colour = guide_legend(override.aes = list(size=2)))
}

p

### Bar plot of frequencies

if(!par_two_methods){ 
  # Assuming you already have the tables
  count_df <- plot_df  %>% dplyr::filter(Relation == "Indirect effect")
  exp_table <- table(count_df$exp_description)
  out_table <- table(count_df$out_description)
  
  # Order the tables from highest to lowest
  exp_table <- exp_table[order(-exp_table)]
  out_table <- out_table[order(-out_table)]
  
  # Plot histograms side by side
  par(mfrow=c(1,2), mar=c(12, 4, 4, 2) + 0.1) # Increased bottom margin
  
  # Plot for exp_description without x-axis labels and capture the bar centers
  bar_centers_exp <- barplot(exp_table, 
                             main="Histogram of exposure", 
                             xlab="", 
                             ylab="Counts", 
                             col="lightblue", 
                             ylim=c(0, max(c(exp_table, out_table))),
                             names.arg=rep("", length(exp_table)))
  
  # Manually add tilted labels for exp_description
  text(x=bar_centers_exp, y=-0.3, srt=45, adj=1, labels=names(exp_table), xpd=TRUE)
  
  # Plot for out_description without x-axis labels and capture the bar centers
  bar_centers_out <- barplot(out_table, 
                             main="Histogram of outcome", 
                             xlab="", 
                             ylab="", 
                             col="lightcoral", 
                             ylim=c(0, max(c(exp_table, out_table))),
                             names.arg=rep("", length(out_table)))
  
  # Manually add tilted labels for out_description
  text(x=bar_centers_out, y=-0.3, srt=45, adj=1, labels=names(out_table), xpd=TRUE)
  
}

################################################
### Plotting mediation forestplot ##############

### Add information
df <- plot_df

# Compute L95 and U95 using se
df$L95 <- df$beta - 1.96 * df$se
df$U95 <- df$beta + 1.96 * df$se

# Determine shape based on pval
df$shape <- ifelse(df$pval < 0.05/multiple_testing, "full_circle", "circle")

df$exposure <- df$exp_description
df$outcome <- df$out_description


# Create a combined factor and reorder it
df <- df %>%
  mutate(combined = paste0(outcome, "_", exposure),
         combined = reorder(combined, -beta * (Relation == "Total effect")))

# Order outcome
df$outcome <- factor(df$outcome, levels = c("Ischemic heart disease","Proxied lifespan",
                                "Total protein","Forced vital capacity (FVC)",
                                "Albumin","IGF-1","Aspartate aminotransferase",
                                "Testosterone","SHBG","Disorders of menstruation"),
                     labels = c("Ischemic heart disease","Lifespan",
                                "Total protein","Forced vital capacity",
                                "Albumin","IGF-1","Aspartate aminotransferase",
                                "Testosterone","SHBG","Disorders of menstruation"))

### Plot
if(!par_two_methods){
  p <- ggplot(data = df) +
    # Background
    geom_stripes(aes(y = combined), odd = "#F4F4F4", even = "white") + 
    geom_vline(xintercept = 0, color = "#2B4050", linewidth = 0.8, lty = "dashed") +
    geom_vline(xintercept = c(-0.2,0.2), color = "#2B4050", linewidth = 0.2, lty = "dotted") +
    # Effect
    geom_effect(aes(x = beta, y = combined, 
                    xmin = L95, xmax = U95, colour = Relation, shape = shape),
                position = ggstance::position_dodgev(height = 0.5), size = 0.55) +
    scale_shape_manual(values = c("full_circle" = 19, "circle" = 1)) +
    scale_color_manual(values=c("#5e548e","#DC63CB","#FE9000"),
                       breaks=c('Total effect', 'Direct effect','Indirect effect')) +
    # Facet
    facet_grid(outcome~., scales = "free", space = "free") +
    # Adjust y-axis to show only exposure
    scale_y_discrete(labels = function(x) gsub(".*_", "", x)) +
    # Axis
    xlab(expression(""*alpha*" (95% CI)")) + ylab("") +
    coord_cartesian(xlim = c(-0.3, 0.3)) +
    # Edit legend
    guides(colour = guide_legend(title = NULL, override.aes = list(size = 2)),
           shape = "none") +
    # Layout
    theme_bw() +
    theme(text = element_text(family = "Arial"),
          axis.text.x = element_text(size = 11),
          axis.title.x = element_text(size = 12, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 9),
          strip.text.y = element_text(size = 9, color = "#2B4050", face = "bold",  angle = 0),
          strip.background = element_rect(fill = "#F4F4F4", color = "#2B4050"),
          # Legend stuff
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          #legend.justification = c(0, 0), # move to the left
          legend.margin = margin(t = -7, r = 0, b = 0, l = 0) # remove space above
    )
  
  ### Export plot
  ggsave(file="/SET/PATH/TO/DIRECTORY/plot1.png", 
         plot = p, width = 6.6, height = 6, dpi = 600)
  ggsave(file="/SET/PATH/TO/DIRECTORY/plot1.svg", 
         plot = p, width = 6.6, height = 6)
}


### Plot Supplemental
if(par_two_methods){
  p <- ggplot(data = df) +
    # Background
    geom_stripes(aes(y = combined), odd = "#F4F4F4", even = "white") + 
    geom_vline(xintercept = 0, color = "#2B4050", linewidth = 0.8, lty = "dashed") +
    geom_vline(xintercept = c(-0.2,0.2), color = "#2B4050", linewidth = 0.2, lty = "dotted") +
    # Effect
    geom_effect(aes(x = beta, y = combined, 
                    xmin = L95, xmax = U95, colour = Relation, shape = shape),
                position = ggstance::position_dodgev(height = 0.5), size = 0.55) +
    scale_shape_manual(values = c("full_circle" = 19, "circle" = 1)) + 
    scale_color_manual(values=c("#5e548e","#DC63CB","#007BFF","#FE9000"),
                       breaks=c('Total effect', 'Direct effect','Indirect effect difference','Indirect effect product')) +
    # Facet
    facet_grid(outcome~., scales = "free", space = "free") +
    # Adjust y-axis to show only exposure
    scale_y_discrete(labels = function(x) gsub(".*_", "", x)) +
    # Axis
    xlab(expression(""*alpha*" (95% CI)")) + ylab("") +
    coord_cartesian(xlim = c(-0.3, 0.3)) +
    # Edit legend
    guides(colour = guide_legend(title = NULL, override.aes = list(size = 2)),
           shape = "none") +
    # Layout
    theme_bw() +
    theme(text = element_text(family = "Arial"),
          axis.text.x = element_text(size = 11),
          axis.title.x = element_text(size = 12, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 9),
          strip.text.y = element_text(size = 9, color = "#2B4050", face = "bold",  angle = 0),
          strip.background = element_rect(fill = "#F4F4F4", color = "#2B4050"),
          # Legend stuff
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          #legend.justification = c(0, 0), # move to the left
          legend.margin = margin(t = -7, r = 0, b = 0, l = 0) # remove space above
    )
  
  ### Export plot
  ggsave(file="/SET/PATH/TO/DIRECTORY/plot2.png", 
         plot = p, width = 6.6, height = 6, dpi = 600)
  
  ggsave(file="/SET/PATH/TO/DIRECTORY/plot2.svg", 
         plot = p, width = 6.6, height = 6)
}


################################################
### Mediation proportion #######################
################################################

plot_df$description <- plot_df$tag

med_prop_df <- plot_df %>%
  mutate(sd = se * sqrt(nsnps)) %>%
  select(c("description","Relation","beta","sd","se")) %>%
  dplyr::filter(Relation != "Direct effect") %>%
  rename(mu = beta)


# Get unique descriptions
unique_descriptions <- unique(med_prop_df$description)
#plot_df$description <- rep(1:(nrow(plot_df)/3), each=3)

# Initialize empty dataframe to store results
results_df <- data.frame(description = character(),
                         mu_Z = numeric(),
                         CI_lower = numeric(),
                         CI_upper = numeric(),
                         stringsAsFactors = FALSE)

# Loop over each unique description
for(trait_interest in unique_descriptions){
  
  mu_X <- as.numeric(med_prop_df[which(med_prop_df$description == trait_interest & 
                                         med_prop_df$Relation == "Indirect effect"),c("mu")]) 
  sigma_X <- as.numeric(med_prop_df[which(med_prop_df$description == trait_interest & 
                                            med_prop_df$Relation == "Indirect effect"),c("se")]) 
  mu_Y <- as.numeric(med_prop_df[which(med_prop_df$description == trait_interest & 
                                         med_prop_df$Relation == "Total effect"),c("mu")])
  sigma_Y <- as.numeric(med_prop_df[which(med_prop_df$description == trait_interest & 
                                            med_prop_df$Relation == "Total effect"),c("se")]) 
  
  mu_Z <- mu_X/mu_Y
  
  # Simulate X and Y
  set.seed(1234)
  samples_X <- rnorm(10000, mean = mu_X, sd = sigma_X)
  samples_Y <- rnorm(10000, mean = mu_Y, sd = sigma_Y)
  
  # Compute Z = X/Y
  samples_Z <- samples_X / samples_Y
  
  # Compute 95% Confidence Interval for Z
  conf_level <- 0.95
  alpha <- 1 - conf_level
  lower_quantile <- alpha/2
  upper_quantile <- 1 - alpha/2
  
  CI_lower <- quantile(samples_Z, probs = lower_quantile)
  CI_upper <- quantile(samples_Z, probs = upper_quantile)
  
  # Append to results dataframe
  results_df <- rbind(results_df, data.frame(description = trait_interest, 
                                             mu_Z = mu_Z, 
                                             CI_lower = CI_lower, 
                                             CI_upper = CI_upper))
  
  # Print results
  print(trait_interest)
  print(paste0("$P_{\text{mediation}} = ",
               round(mu_Z*100,2),"% \text{}95% \text{ CI} [",
               round(CI_lower*100,2),"-",round(CI_upper*100,2),"%]$"))
}

# View results
print(results_df)

### Filter out lifespan and compute mean:
library(stringr)
results_df <- results_df %>%
  dplyr::filter(!str_detect(description, "PROX"))
print(mean(results_df$mu_Z)*100)
