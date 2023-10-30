################################################################################
### Plotting Joint model with parental age                                   ###
### Author: Samuel Moix                                                      ###
### Date: 17.03.2023                                                         ###
################################################################################

### Load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(conflicted)
library(readr)


### Load data
data_folder <- "/SET/PATH/TO/DIRECTORY/"
export_folder <- "/SET/PATH/TO/DIRECTORY/"

# Sample data
df <- as.data.frame(fread(file.path(data_folder,  "filtered_DF.txt"), header = T)) 
metadata <- read.csv(file = "TL_metadata.csv", 
                     sep = ',', header = TRUE)

### Remove blood cancer patients 
if(TRUE){
  # Remove patients with blood cancer
  df <- df[which(df$blood_cancer == 0 | is.na(df$blood_cancer)),]
  df <- select(df, -c("blood_cancer"))
}

df <- df %>% select(c("Zadj_TS","age","age2","sex","father_age_at_birth","mother_age_at_birth",
                      "age_end_education","array"))


### Add corrected LTL
df$TLc <-  residuals(lm(df$Zadj_TS ~ df$age + df$sex + df$age:df$sex + df$age2 + df$age2:df$sex + df$array + df$array:df$sex))
df$TLc <- as.numeric(scale(df$TLc))

### Remove outliers --> None in that case
if(FALSE){
  for (i in c(2,3,5,6,7)) {
    upper <- mean(df[,i], na.rm = T) + sd(df[,i], na.rm = T)*5
    lower <- mean(df[,i], na.rm = T) - sd(df[,i], na.rm = T)*5
    df[,i][which(df[,i] > upper | df[,i] < lower)] <- NA
    rm(upper)
    rm(lower)
  }
}

### Remove LTL outliers
upper <- mean(df[,1], na.rm = T) + sd(df[,1], na.rm = T)*5
lower <- mean(df[,1], na.rm = T) - sd(df[,1], na.rm = T)*5
df <- df[which(df[,1] <= upper & df[,1] >= lower),]

### Scale data
if (TRUE){
  #### To put all traits as numeric
  df[,1:length(df)] <- lapply(df[,1:length(df)], as.numeric) ###
  ### Scale whole data
  df <- as.data.frame(scale(df)) ###
}

################################################################################
### Run analysis
model <- summary(lm(df$Zadj_TS ~ df$age + df$sex + df$age:df$sex + df$father_age_at_birth +
                      df$mother_age_at_birth + df$age_end_education))

# Prepare plot df
plot_df <- as.data.frame(model$coefficients)
plot_df$pheno <- gsub("df\\$", "",rownames(plot_df))

colnames(plot_df) <- c("beta","se","t_val","pval","pheno")
plot_df <- plot_df[-1,]


# Add metadata
plot_df <- merge(plot_df, metadata[,c("pheno","description")], 
                 by = "pheno", all.x = TRUE)

################################################################################
### Add univariable estimate
plot_df <- select(plot_df, -t_val)
plot_df$model <- "Multivariable"


# Create an empty dataframe to store the results
results_df <- data.frame(
  pheno = character(),
  beta = numeric(),
  se = numeric(),
  pval = numeric(),
  description = character(),
  stringsAsFactors = FALSE
)

# Iterate over each pheno in plot_df
for (pheno in plot_df$pheno) {
  if (pheno %in% setdiff(names(df), c('Zadj_TS', 'array', 'age2'))) {
    # Fit the linear regression model
    model <- lm(paste('Zadj_TS ~', pheno), data = df)
    
    # Extract the coefficient, standard error, and p-value for the pheno
    beta <- coef(model)[2] # The [2] is to get the coefficient of the predictor
    se <- summary(model)$coefficients[2, 2] # Standard error for the predictor
    pval <- summary(model)$coefficients[2, 4] # p-value for the predictor
    
    # Extract the description for the pheno from plot_df
    description <- plot_df$description[plot_df$pheno == pheno]
    
    # Append the results to results_df
    results_df <- rbind(results_df, data.frame(pheno, beta, se, pval, description))
  }
}

results_df$model <- "Univariable"
rownames(results_df) <- NULL

# Create an empty dataframe to store the results
results_df_2 <- data.frame(
  pheno = character(),
  beta = numeric(),
  se = numeric(),
  pval = numeric(),
  description = character(),
  stringsAsFactors = FALSE
)

# Iterate over each pheno in plot_df
for (pheno in plot_df$pheno) {
  if (pheno %in% setdiff(names(df), c('Zadj_TS', 'array', 'TLc','age2','age','sex'))) {
    # Fit the linear regression model
    model <- lm(paste('TLc ~', pheno), data = df)
    
    # Extract the coefficient, standard error, and p-value for the pheno
    beta <- coef(model)[2] # The [2] is to get the coefficient of the predictor
    se <- summary(model)$coefficients[2, 2] # Standard error for the predictor
    pval <- summary(model)$coefficients[2, 4] # p-value for the predictor
    
    # Extract the description for the pheno from plot_df
    description <- plot_df$description[plot_df$pheno == pheno]
    
    # Append the results to results_df
    results_df_2 <- rbind(results_df_2, data.frame(pheno, beta, se, pval, description))
  }
}

results_df_2$model <- "Univariable (LTL corrected)"
rownames(results_df_2) <- NULL


plot_df <- rbind(plot_df,results_df, results_df_2)

### Specific modifications
plot_df$description[which(plot_df$pheno == "age:sex")] <- "Age-sex interaction"
plot_df$description[which(plot_df$pheno == "age")] <- "Age"
plot_df$description[which(plot_df$pheno == "age_end_education")] <- "Educational attainment"
plot_df$description[which(plot_df$pheno == "sex")] <- "Female"


### Reorder dataframe by lm estimate
order_df <- plot_df %>%
  dplyr::filter(model == "Multivariable") %>%
  arrange(beta)

plot_df <- plot_df %>%
  group_by(pheno) %>%  
  arrange(factor(pheno, levels = order_df$pheno)) %>%
  ungroup()


################################################################################
### Plotting

plot_exp <- ggforestplot::forestplot(
  df = plot_df,
  name = description,
  estimate = beta,
  pvalue = pval,
  psignif = 0.05/141,
  xlab = expression(""*beta*" (95% CI)"),
  title = "",
  colour = model,
  se = se) +   
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values=c("#985277", "#779852", "black")) +
  theme(axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title=element_blank(),
        legend.position = "top")
plot_exp
### Export plot
ggsave(file=file.path(export_folder, "SupFigure2.png"),
       plot = plot_exp, width = 7, height = 3.5, dpi = 300)
ggsave(file=file.path(export_folder, "SupFigure2.svg"),
       plot = plot_exp, width = 7, height = 3.5)

### Computing pdiff.
# Values from linear regression 
mother_beta_LTL <- plot_df$beta[which(plot_df$pheno == "mother_age_at_birth" & plot_df$model == "Univariable (LTL corrected)")]
mother_se_LTL <- plot_df$se[which(plot_df$pheno == "mother_age_at_birth" & plot_df$model == "Univariable (LTL corrected)")]

father_beta_LTL <- plot_df$beta[which(plot_df$pheno == "father_age_at_birth" & plot_df$model == "Univariable (LTL corrected)")]
father_se_LTL <- plot_df$se[which(plot_df$pheno == "father_age_at_birth" & plot_df$model == "Univariable (LTL corrected)")]

# Extracting beta estimates
mother_beta_model <- plot_df$beta[which(plot_df$pheno == "mother_age_at_birth" & plot_df$model == "Multivariable")]
father_beta_model <- plot_df$beta[which(plot_df$pheno == "father_age_at_birth" & plot_df$model == "Multivariable")]

# Extracting standard errors
mother_se_model <- plot_df$se[which(plot_df$pheno == "mother_age_at_birth" & plot_df$model == "Multivariable")]
father_se_model <- plot_df$se[which(plot_df$pheno == "father_age_at_birth" & plot_df$model == "Multivariable")]

print("Mother")
2*pnorm(-abs((mother_beta_LTL-mother_beta_model)/sqrt(mother_se_LTL**2+mother_se_model**2)),mean = 0, sd = 1)

print("Father")
2*pnorm(-abs((father_beta_LTL-father_beta_model)/sqrt(father_se_LTL**2+father_se_model**2)),mean = 0, sd = 1)
