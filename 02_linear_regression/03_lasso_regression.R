################################################################################
### Lasso regression to find the best predictors of TL                       ###
### Author: Samuel Moix                                                      ###
### Date: 13.04.2022                                                         ###
################################################################################

################################################
# Libraries ####################################
library(glmnet)
library(ggplot2)
library(data.table) # Read data 
library(dplyr)
library(cowplot)
library(ggplotify)
library(gridExtra)
library(forcats)
library(ggforestplot)
library(conflicted)


################################################
### Working directories ########################
data_folder <- "/SET/PATH/TO/DIRECTORY/"
export_folder <- "/SET/PATH/TO/DIRECTORY/"

################################################
### Parameters  ################################

# Weather to make the anylsis sex-stratified
par_sex_strat <- "B" # "B", "F", "M"

# Proportion of NA tolerance
prop_tresh <- 0.07

# Number of replications
par_nb_rep <- 50

# Force to keep traits
par_force_trait <- FALSE
par_traits_to_keep <- c("father_age_at_birth","mother_age_at_birth")

### Whether to remove outliers from data
par_outlier_filter <- TRUE
### Whether to remove patients with blood cancer
par_rm_blood_cancer <- TRUE
### Whether to scale data
par_scale <- TRUE

# Name file
par_exp_name <- function(par_export_name){
  if(par_scale){par_export_name <- paste0(par_export_name, "_scaled")}
  if(par_outlier_filter){par_export_name <- paste0(par_export_name, "_no_outl")}
  if(par_rm_blood_cancer){par_export_name <- paste0(par_export_name, "_no_bc")}
  if(par_force_trait){par_export_name <- paste0(par_export_name, "_t_forced")}
  par_export_name <- paste0(par_export_name, "_THna",prop_tresh*100)
  par_export_name <- paste0(par_export_name, "_", par_sex_strat)
  return(par_export_name)
}


################################################
### Load Data ##################################

df <- as.data.frame(fread(file.path(data_folder,  "filtered_DF.txt"), header = T)) 
metadata <- read.csv(file = "TL_metadata.csv", sep = ',', header = TRUE)

# Select only females
if(par_sex_strat == "F"){
  df <- df %>%
    dplyr::filter(sex == 2) %>%
    dplyr::select(-sex)
}

# Select only males
if(par_sex_strat == "M"){
  df <- df %>%
    dplyr::filter(sex == 1) %>%
    dplyr::select(-sex)
}  

### Adapt and add traits:
# Replace age2 by residuals of age2 by age
df$age2 <- residuals(lm(df$age2 ~ df$age))

# Replace WHR by residuals of WHR by sex
if(FALSE){
  df <- df %>%
    dplyr::filter(!is.na(WHR))
  df$WHR <- residuals(lm(df$WHR ~ df$sex))
}

if(par_sex_strat == "F"){
  df$years_since_menopause <- df$age - df$menopause
  df$years_since_menopause[which(df$had_menopause == 1)] <- 0
}

### Remove blood cancer patients 
if(par_rm_blood_cancer){
  # Remove patients with blood cancer
  df <- df[which(df$blood_cancer == 0 | is.na(df$blood_cancer)),]
  df <- select(df, -c("blood_cancer"))
}

### Remove outliers
if(par_outlier_filter){
  ### Create df only containing continuous traits ###
  continuous_cols <- metadata$type %in% c("continuous", "integer")
  continuous_cols_names <- names(df)[match(metadata$pheno[continuous_cols], names(df))]
  continuous_cols_names <- c(continuous_cols_names,"TLc")
  contvec <- sort(match(continuous_cols_names, names(df)))
  rm(continuous_cols, continuous_cols_names)
  
  # Outliers defined as +/-5SDs
  print("Outliers filtered")
  for (i in contvec) {
    upper <- mean(df[,i], na.rm = T) + sd(df[,i], na.rm = T)*5
    lower <- mean(df[,i], na.rm = T) - sd(df[,i], na.rm = T)*5
    df[,i][which(df[,i] > upper | df[,i] < lower)] <- NA
    rm(upper)
    rm(lower)
  }
  rm(contvec)
}


# Scale data
if(par_scale){
  df <- as.data.frame(scale(df))
  #df <- df[, colSums(is.na(df)) == 0]
}


# Calculate proportion of missing data in each column
missing_prop <- colMeans(is.na(df))

# Identify columns with more than X-% missing data
removed_cols <- names(missing_prop[missing_prop > prop_tresh])

if(par_force_trait){
  removed_cols <- removed_cols[!removed_cols %in% par_traits_to_keep]
}

print("The following traits were removed:")
print(removed_cols)

# Select only columns with less than or equal to proportion threshold of NAs
df <- df[, -which(names(df) %in% removed_cols)]

print("___________________________________")
print("The following traits were kept:")
print(colnames(df))

# Remove rows with NAs
df <- na.omit(df)


# Prepare data
outcome_var <- as.matrix(df$Zadj_TS)
predictor_vars <- as.matrix(df[, !(names(df) %in% "Zadj_TS")])

################################################
### Run regression anylsis #####################

r_squared_vector <- vector()
lasso_coefs_list <- list()
for (iteration in 1:par_nb_rep) {
  
  # Tutorial found at:
  # https://glmnet.stanford.edu/articles/glmnet.html
  
  ### Lasso regression
  fit <- glmnet(predictor_vars, outcome_var, alpha = 1)
  plot1 <- as.grob(function() plot(fit, label = TRUE))
  cv.lasso <- cv.glmnet(predictor_vars, outcome_var, alpha = 1) # alpha = 1 for lasso
  # Lasso path plot
  plot2 <- as.grob(function() plot(cv.lasso))
  # Best lambda 1se
  best_lambda <- cv.lasso$lambda.1se
  print(paste("Best lambda for Lasso:", best_lambda))
  lasso_coefs <- coef(cv.lasso, s = "lambda.1se")
  print("Lasso coefficients:")
  print(lasso_coefs)
  
  ################################################
  ### Compute r-squared ##########################
  
  # Note that Lasso prioritizes model simplicity and feature selection
  predicted_values <- predict(cv.lasso, predictor_vars, s = "lambda.1se")
  SSE <- sum((outcome_var - predicted_values)^2) # Sum of Squared Errors
  SST <- sum((outcome_var - mean(outcome_var))^2) # Total Sum of Squares
  r_squared <- 1 - (SSE / SST)
  print(paste("R-squared value:", r_squared))
  
  ################################################
  ### Plotting ###################################
  
  # Lasso coefficients plot
  lasso_coefs <- as.data.frame(as.matrix(lasso_coefs))
  lasso_coefs$Variable <- rownames(lasso_coefs)
  colnames(lasso_coefs) <- c("Coefficient", "Variable")
  lasso_coefs <- lasso_coefs[-1,] # Remove the intercept
  
  ### Add description
  
  colnames(metadata)[which(colnames(metadata) == "pheno")] <- "Variable"
  if(par_sex_strat == "F"){
    # Create a new row
    new_row <- data.frame(Variable = "years_since_menopause",
                          FieldID = NA,
                          MR_name = NA,
                          description = "Years since menopause",
                          mean = NA,
                          sd = NA,
                          N = NA,
                          Disease_cases = NA,
                          type = NA,
                          GWAS_source = NA,
                          Codd = NA,
                          Comment = NA,
                          Category = NA,
                          Group = NA)
    # Add the new row to the metadata dataframe
    metadata <- rbind(metadata, new_row)
    
  }
  
  lasso_coefs <- merge(lasso_coefs, metadata[,c("Variable","description")], by="Variable")
  
  lasso_coefs <- lasso_coefs[which(lasso_coefs$Coefficient != 0),]
  
  plot3 <- ggplot(lasso_coefs, aes(x = reorder(description, Coefficient), y = Coefficient)) +
    geom_col() +
    coord_flip() +
    theme_minimal(base_size = 17) +
    labs(title = "Lasso Coefficients", x = "", y = "Coefficient") +
    theme(axis.text.y = element_text(size = 19)) +
    annotate("text", x = Inf, y = -0.05, label = paste("R-squared:", round(r_squared, 4)), 
             hjust = 1, vjust = 1, size = 4)
  
  export_name <- par_exp_name("_lasso_plot")
  
  ggsave(filename = paste0(export_folder,iteration,export_name,".pdf"), plot = grid.arrange(
    plot_grid(plot1, plot2, ncol = 1),
    plot3,
    ncol = 2,
    widths = c(0.5, 0.5),
    layout_matrix = rbind(c(1, 2), c(1, 2))
  ), width = 18, height = 9, units = "in")
  
  ### Storing results of each iteration
  r_squared_vector <- c(r_squared_vector, r_squared)
  lasso_coefs_list[[iteration]] <- lasso_coefs
}

################################################
### Merge information of each iteration ########

# Create a list to store the counts of each variable
variable_counts <- list()

# Create a list to store the sum of coefficients for each variable
coefficient_sums <- list()

# Loop over each dataframe in the list
for(df_coefs_i in lasso_coefs_list){
  # Loop over each row in the dataframe
  for(i in 1:nrow(df_coefs_i)){
    # Get the variable name
    variable <- df_coefs_i$Variable[i]
    
    # If the variable is not in the counts list, add it
    if(is.null(variable_counts[[variable]])){
      variable_counts[[variable]] <- 1
      coefficient_sums[[variable]] <- df_coefs_i$Coefficient[i]
    } else {
      # If the variable is already in the counts list, increment the count
      variable_counts[[variable]] <- variable_counts[[variable]] + 1
      coefficient_sums[[variable]] <- coefficient_sums[[variable]] + df_coefs_i$Coefficient[i]
    }
  }
}

# Calculate the threshold for being in 95% of dataframes
threshold <- length(lasso_coefs_list) * 0.95

# Create a new dataframe to store the final results
result_df <- data.frame(Variable = character(), Coefficient = numeric(), description = character())

# Loop over each variable in the counts list
for(variable in names(variable_counts)){
  # If the variable is in at least 95% of dataframes, add it to the result dataframe
  if(variable_counts[[variable]] >= threshold){
    # Find the description of the variable from the first dataframe in the list
    description <- lasso_coefs_list[[1]]$description[lasso_coefs_list[[1]]$Variable == variable]
    
    # Calculate the average coefficient
    average_coefficient <- coefficient_sums[[variable]] / variable_counts[[variable]]
    
    # Add the variable to the result dataframe
    result_df <- rbind(result_df, data.frame(Variable = variable, Coefficient = average_coefficient, description = description))
  }
}

# Print the result dataframe
print(result_df)
lasso_coefs<- result_df

plot3 <- ggplot(lasso_coefs, aes(x = reorder(description, Coefficient), y = Coefficient)) +
  geom_col() +
  coord_flip() +
  theme_minimal(base_size = 17) +
  labs(title = "Average lasso Coefficients", x = "", y = "Coefficient") +
  theme(axis.text.y = element_text(size = 19)) +
  annotate("text", x = Inf, y = -0.05, label = paste("Avg. R-squared:", round(mean(r_squared_vector), 4)), 
           hjust = 1, vjust = 1, size = 4)

export_name <- par_exp_name("it_avg_lasso_plot")

ggsave(filename = paste0(export_folder,iteration,export_name,".pdf"), plot = grid.arrange(
  plot_grid(plot1, plot2, ncol = 1),
  plot3,
  ncol = 2,
  widths = c(0.5, 0.5),
  layout_matrix = rbind(c(1, 2), c(1, 2))
), width = 18, height = 9, units = "in")


################################################
### Make barplot ###############################

# Initialize an empty vector to store the selection frequency of each trait
trait_frequency <- vector()

# Loop over each trait in the counts list
for(trait in names(variable_counts)){
  # Calculate the selection frequency of the trait
  selection_frequency <- variable_counts[[trait]] / length(lasso_coefs_list)
  
  # Add the selection frequency to the trait_frequency vector
  trait_frequency <- c(trait_frequency, selection_frequency)
}

# Create a dataframe with the trait names and their selection frequencies
trait_frequency_df <- data.frame(Variable = names(variable_counts), SelectionFrequency = trait_frequency)

# Sort the dataframe by selection frequency
trait_frequency_df <- trait_frequency_df[order(trait_frequency_df$SelectionFrequency, decreasing = TRUE), ]

# Create a new variable 'fill' for coloring bars
trait_frequency_df$fill <- rep(c("deepskyblue4", "dodgerblue2"), length.out = nrow(trait_frequency_df))

# Create the plot
ggplot(trait_frequency_df, aes(x = reorder(Variable, SelectionFrequency), y = SelectionFrequency, fill = fill)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() + # Flipping the x and y axes
  theme_minimal() +
  theme(legend.position = "none") + # Removing legend
  scale_fill_manual(values = c("deepskyblue4", "dodgerblue2")) + # Apply colors manually
  labs(x = "Variable", y = "SelectionFrequency", 
       title = "Selection Frequency of Variables") # Labels

# Save frequencies
write.table(trait_frequency_df, 
            file.path(export_folder, paste0("00_",par_sex_strat,"_",prop_tresh*100,"_iter_",par_nb_rep,"_freq_lasso.txt")),
            sep = "\t", row.names = FALSE, quote = FALSE)

################################################
###### Linear model with selected variables ####

# Read file if already run
# Define the file name based on the same parameters
file_name <- file.path(export_folder, paste0("00_", par_sex_strat, "_", prop_tresh * 100, "_iter_", par_nb_rep, "_freq_lasso.txt"))

# Read the file into the dataframe
trait_frequency_df <- read.table(file_name, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Note: The `header = TRUE` argument indicates that the first row of the file contains column names.


# Construct the formula
formula_str <- paste("Zadj_TS ~", 
                     paste(trait_frequency_df$Variable[which(trait_frequency_df$SelectionFrequency > 0.95)], collapse = " + "))

# Fit the model
result <- lm(formula = as.formula(formula_str), data = df)
cat(paste("Formula:",formula_str,"\n\nAdj.r-squared", summary(result)$adj.r.squared))

# Create dataframe from your coefficients
forestplot_data <- as.data.frame(summary(result)$coefficients)

# Add variable names to the dataframe
forestplot_data$Variable <- rownames(forestplot_data)

# Rename the columns
names(forestplot_data) <- c("Estimate", "Std.Error", "t.value", "Pval", "Variable")
colnames(metadata)[which(colnames(metadata) == "pheno")] <- "Variable"

forestplot_data <- dplyr::filter(forestplot_data, Variable != "(Intercept)")
forestplot_data$Variable <- sub("df\\$", "", forestplot_data$Variable)
forestplot_data <- merge(forestplot_data, metadata[,c("Variable","description")], by = "Variable", all.x = TRUE)

# Rename labels
forestplot_data$description[which(forestplot_data$description == "Age completed full time education")] <- "EA"
forestplot_data$description[which(forestplot_data$description == "Sex")] <- "Female"
forestplot_data$description[which(forestplot_data$description == "Age when attended assessment centre")] <- "Age"
forestplot_data$description[which(forestplot_data$description == "Mean corpuscular haemoglobin")] <- "MCH"
forestplot_data$description[which(forestplot_data$description == "Waist-to-hip ratio (WHR)")] <- "WHR"

# Reorder
forestplot_data <- dplyr::arrange(forestplot_data, desc(Estimate))
 
# Use ggforestplot to create the plot
lasso_jlm_plot <- ggforestplot::forestplot(
  df = forestplot_data,
  name = description,
  estimate = Estimate,
  pvalue = Pval,
  xlab = expression(""*beta*" (95% CI)"),
  title = "",
  se = Std.Error,
  xlim = c(-0.2, 0.1)) +
  # annotate("text", x = -0.05, y = 10, label = paste("R-squared:", round(0.054, 4)), 
  #          hjust = 1, vjust = 1, size = 4) +
  theme(
    axis.text.y = element_text(size = 16, hjust = 1),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 17)) 
lasso_jlm_plot

require(svglite)
ggsave(file="lasso_jlm.svg", plot=lasso_jlm_plot, width=7, height=6)
