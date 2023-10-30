#################################################################################################################
### Data preparation                                                                                          ###
### Author: Samuel Moix                                                                                       ###
### Date: 04.04.2022                                                                                          ###
#################################################################################################################


################################################
### Libraries ##################################
library(data.table) # Read data 
library(dplyr)
library(tibble)
library(tidyr)
library(conflicted)

################################################
### Working directories ########################
data_folder = "/SET/PATH/TO/DIRECTORY"
export_folder = "/SET/PATH/TO/DIRECTORY"


################################################
### Data loading and filtering #################
df <- as.data.frame(fread(file.path(data_folder,  "complete_DF.txt"), header = T)) 

### Remove PCs and arrays and set factors as factors
df <- df[,c(1:which(colnames(df) == "array"), 
            which(colnames(df) == "phasing_Y"):length(df))]
df <- df %>% 
  select(-"phasing_Y")

### Remove patients below 40 and above 70 years old
df <- dplyr::filter(df, age >= 40 & age <= 70 )

################################################
### LATER ADDED TRAITS #########################

### Relocate cholesterol medication with factors

df <- df %>% mutate(cholest_medication = coalesce(cholest_medication_M, cholest_medication_F)) %>%
  select(-c(cholest_medication_M,cholest_medication_F))


df <- relocate(df, cholest_medication, .after = alcool_frq)
### Set -7 to 0 as -7 corresponds to "None of the above"
df$cholest_medication[which(df$cholest_medication == -7)] <- 0
### Change into binary trait (written like that so that we see what there is)
df$cholest_medication[which(df$cholest_medication > 0 & 
                              df$cholest_medication != 2 & # 2: Blood pressure medication
                              df$cholest_medication != 3 & # 3: Insulin
                              df$cholest_medication != 4 & # 4: Hormone replacement therapy
                              df$cholest_medication != 5)] <- 1 # 5: Oral contraceptive pill or minipill
### Set 0 to 2 and 3 as they do not affect cholesterol levels
df$cholest_medication[which(df$cholest_medication > 1)] <- 0

### Weekly alcool use (attention average was taken if multiple values given,

## Merge alcool data into one trait
# df$alcool_weekly <- df %>% 
#   select(contains("alcool_0")) %>%
#   mutate_all(round, 1) %>% # As average has been calculated round numbers
#   rowSums(na.rm = T)  # Sum each consumption units

# Keep non-full answered
temp_alc <- df %>% 
  select(contains("alcool_0")) %>%
  mutate_all(round, 1)

# Get rows full of NA's
temp_full_na <- which(rowSums(is.na(temp_alc)) == ncol(temp_alc))
# Replace NA's by average
temp_alc <- temp_alc %>% 
  replace_na(replace = as.list(round(colMeans(temp_alc, na.rm = T),1)))
# Put full row NA's back 
temp_alc[temp_full_na,] <- NA
# If one wants to not count the NA's
#df$alcool_weekly <- rowSums(temp_alc, na.rm=TRUE) * ifelse(rowSums(is.na(temp_alc)) == ncol(temp_alc), NA, 1)

# Sum the alcool quantities
df$alcool_weekly <- rowSums(temp_alc)
rm(temp_alc)
rm(temp_full_na)

# Relocate trait and remove others  
df <- relocate(df, alcool_weekly, .after = alcool_frq)
#df <- select(df, !contains("alcool_0"))
df <- relocate(df, contains("alcool_0"), .after = alcool_frq)
df <- df %>% 
  mutate(across(c("alcool_01_bc", "alcool_02_cw", "alcool_03_fw", "alcool_04_o",
                  "alcool_05_rw", "alcool_06_s"), ~round(.x, 1)))

# Relocate handgrip
df <- relocate(df, hand_grip_r, .before = alcool_frq) 
df <- relocate(df, hand_grip_l, .before = alcool_frq) 

# Modify age_end_education
df$age_end_education[which(df$age_end_education < 0)] <- NA
df$age_end_education[which(df$age_end_education < 14)] <- 14
df$age_end_education[which(df$age_end_education > 19)] <- 19
df$age_end_education[which(df$qualifications == 1)] <- 19 # 1 == university
df <- select(df, -qualifications) # Remove qualification trait

### LATER ADDED TRAITS #########################
################################################

# Set NA's for Do not know and prefer not to answer
no_negatives_index <- c(1:length(df))
no_negatives_index <- no_negatives_index[!no_negatives_index %in% which(colnames(df) %in% c("Zadj_TS","TDI"))]
# Attention sets NA to all negative values to traits except TL and TDI
df[,no_negatives_index][df[,no_negatives_index] < 0] <- NA

# Set NAs not answered (not needed if already done during extraction)
df$water_intake[which(df$water_intake < 0)] <- NA # loosing -10 which = less than one glass/day


# Set variables as factors
fact_start_pos <- which(colnames(df) == "alcool_frq")
df[,fact_start_pos:length(df)] <- lapply(df[,fact_start_pos:length(df)] , factor)
df$sex <- as.factor(df$sex) # Set sex as factor

### Put NAs in Sex specific traits
# In diseases c("BC", "OC", "PC", "endometriosis", "menstruation")
# Menopause, menarche already ok

df$BC[which(df$sex == 1)] <- NA
df$OC[which(df$sex == 1)] <- NA
df$endometriosis[which(df$sex == 1)] <- NA
df$menstruation[which(df$sex == 1)] <- NA
df$PC[which(df$sex == 2)] <- NA

### Set NA's for had menopause and relocate
# 2	Not sure - had a hysterectomy
# 3	Not sure - other reason
df <- relocate(df, had_menopause, .before = menopause)
df[df[, "had_menopause"] %in% c(2, 3), "had_menopause"] <- NA
df$had_menopause <- as.integer(df$had_menopause)

### Composite traits  

### Add hand grip average and waist to hip ratio
df <- df %>%
  add_column(hand_grip_avg = (df$hand_grip_l+df$hand_grip_r)/2, .after = "blood_cancer") %>%
  select(-hand_grip_l, -hand_grip_r) %>%
  add_column(WHR = (df$waist_cir/df$hip_cir), .after = "hip_cir")

### Add age squared
df <- df %>% add_column(age2 = (df$age)**2 , .after = "age") 

### Add reproductive life span
df <- df %>% add_column(reproductive_lp = df$menopause-df$menarche, .after = "menarche") 

# Stratified WHRadjBMI 
df_F <- subset(df, sex == 2 & !is.na(df$WHR) & !is.na(df$BMI))
res_F <- residuals(lm(df_F$WHR ~ df_F$BMI, na.action=na.exclude))
df[which(df$sex == 2 & !is.na(df$WHR) & !is.na(df$BMI)), "WHRadjBMI"] <- res_F
rm(df_F);rm(res_F)
df_M <- subset(df, sex == 1 & !is.na(df$WHR) & !is.na(df$BMI))
res_M <- residuals(lm(df_M$WHR ~ df_M$BMI, na.action=na.exclude))
df[which(df$sex == 1 & !is.na(df$WHR) & !is.na(df$BMI)), "WHRadjBMI"] <- res_M
df$WHRadjBMI[is.na(df$WHR) | is.na(df$BMI)] <- NA # Not necessary but to be sure
rm(df_M);rm(res_M)
df <- df %>%
  relocate(WHRadjBMI, .after = "WHR")


# Define and relocate father and mother age at birth and remove parental age
df$father_age <- as.integer(as.character(df$father_age))
df$father_age_at_birth <- df$father_age - df$age
df <- relocate(df, father_age_at_birth, .before = father_age_death) 
df <- select(df, -father_age) # Remove father_age as not needed 

df$mother_age <- as.integer(as.character(df$mother_age))
df$mother_age_at_birth <- df$mother_age - df$age
df <- relocate(df, mother_age_at_birth, .before = mother_age_death) 
df <- select(df, -mother_age) # Remove father_age as not needed 

# UNCOMMENT Not stratified
# df <- df %>% add_column(WHRadjBMI_G = residuals(lm(df$WHR ~ df$BMI, na.action=na.exclude)), .after = "WHR") 

# ### Add proxied lifespan (version based on sex)
# mean_to_male <- mean(df$father_age_death, na.rm = T)-mean(c(df$father_age_death, df$mother_age_death), na.rm = T)
# female_to_male <- mean(df$father_age_death, na.rm = T)-mean(df$mother_age_death, na.rm = T)
# mean_to_female <- mean(df$mother_age_death, na.rm = T)-mean(c(df$father_age_death, df$mother_age_death), na.rm = T)
# 
# 
# 
# df <- df %>% add_column(prox_life_span = ifelse(!is.na(df$age_death), df$age_death,
#                                                ifelse(df$sex == 1, 
#                                                       ifelse((!is.na(df$father_age_death) & !is.na(df$mother_age_death)),
#                                                              (df$father_age_death+df$mother_age_death)/2+mean_to_male,
#                                                              ifelse(!is.na(df$father_age_death),
#                                                                     df$father_age_death,
#                                                                     ifelse(!is.na(df$mother_age_death), df$mother_age_death+female_to_male, NA))), 
#                                                       ifelse((!is.na(df$father_age_death) & !is.na(df$mother_age_death)),
#                                                              (df$father_age_death+df$mother_age_death)/2+mean_to_female,
#                                                              ifelse(!is.na(df$mother_age_death),
#                                                                     df$mother_age_death,
#                                                                     ifelse(!is.na(df$father_age_death), df$father_age_death-female_to_male, NA)))))
#                         , .after = "age_death") 

### Proxied lifespan scaled by sex

# Scale between females
df$age_death_female_scale <- df$age_death
df$age_death_female_scale[which(df$sex == 1)] <- NA
df$age_death_female_scale <- as.vector(scale(df$age_death_female_scale))
# Scale between males
df$age_death_male_scale <- df$age_death
df$age_death_male_scale[which(df$sex == 2)] <- NA
df$age_death_male_scale <- as.vector(scale(df$age_death_male_scale))
# Merge
df <- mutate(df, age_death_scale = coalesce(age_death_male_scale, age_death_female_scale)) 
# Scale parental death
df$father_age_death_scale <- as.vector(scale(df$father_age_death))
df$mother_age_death_scale <- as.vector(scale(df$mother_age_death))

# Take age death if present, mean if two parental deaths, or either that is not NA
df <- df %>% add_column(prox_life_span = ifelse(!is.na(df$age_death_scale), df$age_death_scale,
                                                ifelse(!is.na(df$father_age_death_scale) & !is.na(df$mother_age_death_scale),
                                                       (df$father_age_death_scale + df$mother_age_death_scale)/2,
                                                       ifelse(!is.na(df$father_age_death_scale), df$father_age_death_scale,
                                                              ifelse(!is.na(df$mother_age_death_scale),df$mother_age_death_scale,NA))
                                                      )
                                                )
                                              
                        , .after = "age_death")

# Remove computing columns
df <- df %>% select(-c("age_death_scale","father_age_death_scale","mother_age_death_scale",
                 "age_death_female_scale","age_death_male_scale"))



### Add age of first birth of primiparous women
df$age_birth_primi <- as.integer(as.character(df$age_birth_primi))
df <- mutate(df, age_first_birth = coalesce(age_first_birth, age_birth_primi)) 
df <- mutate(df, age_last_birth = coalesce(age_last_birth, age_birth_primi)) 
df <- select(df, -age_birth_primi)

### Add smoking cessation based on smoking status
# As defined "Binary phenotype with current smokers coded as "2" and former smokers 
# coded as "1", and never smokers are coded as missing"
df <- df %>%
  mutate(smoking_cessation = smoking_status) %>%
  mutate(smoking_cessation = as.numeric(as.character(smoking_cessation))) %>%
  mutate(smoking_cessation = na_if(smoking_cessation, 0)) %>%
  mutate(smoking_cessation = as.factor(smoking_cessation)) %>%
  relocate(smoking_cessation, .after = smoking_status)

### Adapt non binary factorial traits

# Household income
# 1	Less than 18,000  # 2	18,000 to 30,999 # 3	31,000 to 51,999 
# 4	52,000 to 100,000 # 5	Greater than 100,000
df$household_income <- as.integer(df$household_income)

# Time owning a mobile phone
# 0	Never used mobile phone at least once per week
# 1	One year or less # 2	Two to four years # 3	Five to eight years # 4	More than eight years
df$mobile_phone <- as.integer(df$mobile_phone)

# Beef and cheese intake
# 0	Never # 1	Less than once a week # 2	Once a week # 3	2-4 times a week
# 4	5-6 times a week # 5	Once or more daily
df$beef <- as.integer(df$beef)
df$cheese <- as.integer(df$cheese)

# Alcool frequency 
# 5	Daily or almost daily # 4	Three or four times a week # 3	Once or twice a week
# 2	One to three times a month # 1	Special occasions only
val <- c(5:1)
df$alcool_frq <- as.integer(df$alcool_frq)
df$alcool_frq <- val[df$alcool_frq]

# Playing computer games
# 0	Never/rarely # 1	Sometimes # 2	Often
df$computer_games <- as.integer(df$computer_games)

# Use of sun/uv protection
# 1	Never/rarely # 2	Sometimes # 3	Most of the time # 4	Always # 5	Do not go out in sunshine
df$sun_protection <- as.integer(df$sun_protection)

#Relative age of first facial hair
# 1	Younger than average # 2	About average age # 3	Older than average
df$facial_hair <- as.integer(df$facial_hair)

# Balding 
# 1	Pattern 1 # 2	Pattern 2 # 3	Pattern 3 # 4	Pattern 4 (1 not, 4 is really bald)
df$balding <- as.integer(df$balding)

# Diseases burden
df$DiseaseBurden <- as.integer(df$DiseaseBurden)

# Set 0 and 1 for remaining factors
df$blood_cancer <- df$blood_cancer*1
df$blood_cancer <- as.factor(df$blood_cancer)

df$array <- recode_factor(df$array, UKBL = "0", UKBB = "1")


# Set NAs to blood cancer (patients with another cancer) as this wasn't done initially
df$blood_cancer[which(df$blood_cancer == 0 & (is.na(df$KC) | df$KC == 2))] <- NA

### Write data to table for further analysis:

fwrite(df, file.path(export_folder,  "filtered_DF.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
