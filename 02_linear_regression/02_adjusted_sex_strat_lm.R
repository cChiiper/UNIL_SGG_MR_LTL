#################################################################################################################
### Table for LTL ~ age sex stratified adjusted linear models                                                 ###
### Author: Samuel Moix                                                                                       ###
### Date: 24.03.2022                                                                                          ###
#################################################################################################################


################################################
### Libraries ##################################
library(data.table) # Read data 
library(ggplot2) # Plotting basics
library(ggExtra) # Extra plotting XP
library(ggpubr) # Plotting regression stuff
library(dplyr)
library(conflicted)

################################################
### Working directories ########################
data_folder = "/SET/PATH/TO/DIRECTORY/"

################################################
### Parameters  ################################

### Whether to remove patients with blood cancer
par_rm_blood_cancer <- TRUE

### V to correct for lifestyles, H to correct for hormones, I to correct for triglycerides
par_vege <- "V" #V,H,I


################################################
### Data loading and filtering #################
df <- as.data.frame(fread(file.path(data_folder,  "filtered_DF.txt"), header = T)) 


### Remove blood cancer patients 
if(par_rm_blood_cancer){
  # Remove patients with blood cancer
  df <- df[which(df$blood_cancer == 0 | is.na(df$blood_cancer)),]
  df <- select(df, -c("blood_cancer"))
}

df$Zadj_TS <- scale(df$Zadj_TS)
df$sex <- as.factor(df$sex) # Set sex as factor

################################################
### Analysis of corrected sex-stratified lm ####
conflicts_prefer(dplyr::filter)

if(par_vege == "V"){
  df <- df %>%
    filter(!is.na(alcool_frq)) %>%
    filter(!is.na(alcool_weekly)) %>%
    filter(!is.na(smoking_status)) %>%
    filter(!is.na(fruit)) %>%
    filter(!is.na(beef)) %>%
    filter(!is.na(vegetable))
}

if(par_vege == "H"){
  df <- df %>%
    filter(!is.na(testosterone)) %>%
    filter(!is.na(SHBG))
}

if(par_vege == "I"){
  df <- df %>%
    filter(!is.na(TG)) 
}

df$Zadj_TS <- scale(df$Zadj_TS)
#df$age <- scale(df$age)

if(par_vege == "V"){
  ### Add corrected trait
  df$TLc <- residuals(lm(df$Zadj_TS ~ df$alcool_frq + df$alcool_weekly +
                                df$smoking_status + df$fruit + df$beef +
                                df$vegetable))
}

if(par_vege == "H"){
  df$TLc <- residuals(lm(df$Zadj_TS ~ df$testosterone + df$SHBG))
}

if(par_vege == "I"){
  df$TLc <- residuals(lm(df$Zadj_TS ~ df$TG))
}


df$TLc <- scale(df$TLc)


df_M_ls <- filter(df, sex == 1) %>% select(-sex)
df_F_ls <- filter(df, sex == 2) %>% select(-sex)

### Regression
# Non-adjusted regression
male_reg <- lm(df_M_ls$Zadj_TS ~ df_M_ls$age)
male_b <- male_reg$coefficients[[2]]
male_se <- coef(summary(male_reg))[2, 2]

female_reg <- lm(df_F_ls$Zadj_TS ~ df_F_ls$age)
female_b <- female_reg$coefficients[[2]]
female_se <- coef(summary(female_reg))[2, 2]

pdiff_1 <- 2*pnorm(-abs((male_b-female_b)/sqrt(male_se**2+female_se**2)),mean = 0, sd = 1)
print(paste("Male b:", round(male_b,3), "male se:", round(male_se,4),
            "Female b:", round(female_b,3), "female se:", round(female_se,4),
            "P_diff:", pdiff_1), sep=" ")


# Adjusted regression
male_reg_c <- lm(df_M_ls$TLc ~ df_M_ls$age)
male_b_c <- male_reg_c$coefficients[[2]]
male_se_c <- coef(summary(male_reg_c))[2, 2]

female_reg_c <- lm(df_F_ls$TLc ~ df_F_ls$age)
female_b_c <- female_reg_c$coefficients[[2]]
female_se_c <- coef(summary(female_reg_c))[2, 2]

pdiff_2 <- 2*pnorm(-abs((male_b_c-female_b_c)/sqrt(male_se_c**2+female_se_c**2)),mean = 0, sd = 1)
print(paste("Male b:", round(male_b_c,3), "male se:", round(male_se_c,4),
            "Female b:", round(female_b_c,3), "female se:", round(female_se_c,4),
            "P_diff:", pdiff_2), sep=" ")

table_1 <- data.frame(male_b  = c(round(male_b,3),round(male_b_c,3)),
                      male_se = c(round(male_se,4),round(male_se_c,4)),
                      female_b = c(round(female_b,3),round(female_b_c,3)),
                      female_se = c(round(female_se,4),round(female_se_c,4)),
                      pdiff = c(pdiff_1,pdiff_2))
rownames(table_1) <- c("TL","TLc")


### Difference between adjusted and non-adjusted
pdiff_m <- 2*pnorm(-abs((male_b_c-male_b)/sqrt(male_se_c**2+male_se**2)),mean = 0, sd = 1)
pdiff_f <- 2*pnorm(-abs((female_b_c-female_b)/sqrt(female_se_c**2+female_se**2)),mean = 0, sd = 1)

table_1 ; print(paste("pdiff males:", pdiff_m)) ; print(paste("pdiff females:", pdiff_f)) ; print(paste("N:", nrow((df))))

