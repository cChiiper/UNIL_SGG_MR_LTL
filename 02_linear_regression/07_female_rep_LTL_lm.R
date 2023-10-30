################################################################################
### Analysis to test whether pregnancy has an effect on TL                   ###
### Author: Samuel Moix                                                      ###
### Date: 11.12.2022                                                         ###
################################################################################

################################################
### Libraries ##################################
library(data.table) # Read data 
library(ggplot2) # Plotting basics
library(dplyr)
library(conflicted)

################################################
### Working directories ########################
data_folder = "/SET/PATH/TO/DIRECTORY/"
export_folder = "/SET/PATH/TO/DIRECTORY/"

################################################
### Parameters  ################################

###  Number of independent traits
indep_trait <- 141
### Whether to remove patients with blood cancer
par_rm_blood_cancer <- TRUE
### Remove outliers
par_outlier_filter <- TRUE

################################################
### Load Data ##################################

df <- as.data.frame(fread(file.path(data_folder,  "filtered_DF.txt"), header = T)) 
metadata <- read.csv(file = "TL_metadata.csv", sep = ',', header = TRUE)


### Remove blood cancer patients 
if(par_rm_blood_cancer){
  # Remove patients with blood cancer
  df <- df[which(df$blood_cancer == 0 | is.na(df$blood_cancer)),]
  df <- select(df, -c("blood_cancer"))
}

### Linear regression model TL corrected 
df$TLc <- residuals(lm(df$Zadj_TS ~ df$age + df$sex + df$age:df$sex + df$age2 + df$age2:df$sex + df$array + df$array:df$sex))

df <- df %>%
  select(c("Zadj_TS","TLc","sex","age","nb_birth", "age_last_birth", "age_first_birth","had_menopause","menopause"))

# Outliers defined as ?5SDs
if(par_outlier_filter){
  print("Outliers filtered")
  for (i in c("Zadj_TS","TLc","nb_birth", "age_last_birth", "age_first_birth","had_menopause","menopause")) {
    upper <- mean(df[,i], na.rm = T) + sd(df[,i], na.rm = T)*5
    lower <- mean(df[,i], na.rm = T) - sd(df[,i], na.rm = T)*5
    df[,i][which(df[,i] > upper | df[,i] < lower)] <- NA
    rm(upper)
    rm(lower)
  }
}

### Remove TL outliers
df <- df[which(!is.na(df$Zadj_TS)),]

################################################
### Scale data and add column ##################

df$Zadj_TS <- as.vector(scale(df$Zadj_TS))
df$TLc <- as.vector(scale(df$TLc))
df$age_last_birth_scaled <- as.vector(scale(df$age_last_birth))


### Add ever had a pregnancy
df$ever_pregnant <- ifelse(df$nb_birth > 0, 1, 0)

################################################
### Calculating regression coef ################

### P-value function
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

### Regress against ever had kids or not
df$ever_pregnant_scaled <- as.vector(scale(df$ever_pregnant))
model <- lm(df$TLc ~ df$ever_pregnant_scaled)
print(lmp(model))
print(model)
print(confint(model)[2,1])
print(confint(model)[2,2])

# Welch two-sample t-test
wtst_test <- t.test(TLc ~ ever_pregnant, data = df)
print(wtst_test)
wtst_mean_0 <- wtst_test$estimate[[1]]
wtst_mean_1 <- wtst_test$estimate[[2]]
(wtst_mean_0-wtst_mean_1)/wtst_mean_1*100


################################################
### TL ~ age_last_birth by age #################

df_young_female <- df[which(df$sex == 2 & df$age <= 50),]
df_older_female <- df[which(df$sex == 2 & df$age > 50),]

### Check all together
lm(df$TLc ~ df$age_last_birth_scaled)

### Compare models
young_reg <- lm(df_young_female$TLc ~ df_young_female$age_last_birth_scaled)
young_b <- young_reg$coefficients[[2]]
young_se <- coef(summary(young_reg))[2, 2]

older_reg <- lm(df_older_female$TLc ~ df_older_female$age_last_birth_scaled)
older_b <- older_reg$coefficients[[2]]
older_se <- coef(summary(older_reg))[2, 2]

### Show results
print("TLc ~ age last birth by age group")
pdiff <- 2*pnorm(-abs((young_b-older_b)/sqrt(young_se**2+older_se**2)),mean = 0, sd = 1)
print(paste("Young b:", round(young_b,3), "Young se:", round(young_se,4),
            "Older b:", round(older_b,3), "Older se:", round(older_se,4),
            "P_diff:", pdiff), sep=" ")

rm(df_young_female) 
rm(df_older_female)
rm(young_reg,young_b,young_se,older_reg,older_b,older_se)

################################################
### TL ~ age by sex (corrected by age last birth

df_M <- df[which(df$sex == 1),]
df_F <- df[which(df$sex == 2),]

### Adjust female length by age_last_birth

df_F <- df_F[which(!is.na(df_F$age_last_birth_scaled)),]
df_F$Zadj_TS_corrected <- residuals(lm(df_F$Zadj_TS ~ df_F$age_last_birth_scaled))

### Scale TL
df_M$Zadj_TS <- as.vector(df_M$Zadj_TS)
df_F$Zadj_TS <- as.vector(df_F$Zadj_TS)
df_F$Zadj_TS_corrected <- as.vector(df_F$Zadj_TS_corrected)

### Compare models
male_reg <- lm(df_M$Zadj_TS ~ df_M$age)
male_b <- male_reg$coefficients[[2]]
male_se <- coef(summary(male_reg))[2, 2]

female_reg <- lm(df_F$Zadj_TS_corrected ~ df_F$age)
#female_reg <- lm(df_F$Zadj_TS ~ df_F$age)
female_b <- female_reg$coefficients[[2]]
female_se <- coef(summary(female_reg))[2, 2]

### Show results
print("Female adjusted for age_last_birth")
pdiff <- 2*pnorm(-abs((male_b-female_b)/sqrt(male_se**2+female_se**2)),mean = 0, sd = 1)
print(paste("Male b:", round(male_b,3), "male se:", round(male_se,4),
            "Female b:", round(female_b,3), "female se:", round(female_se,4),
            "P_diff:", pdiff), sep=" ")

rm(male_reg,male_b,male_se,female_reg,female_b,female_se)

################################################
### TL ~ age by sex (without pregnant women) ###

df_M <- df[which(df$sex == 1),]
df_F <- df[which(df$sex == 2 & df$ever_pregnant == 0),]

### Scale TL
df_M$Zadj_TS <- as.vector(df_M$Zadj_TS)
df_F$Zadj_TS <- as.vector(df_F$Zadj_TS)

### Compare models
male_reg <- lm(df_M$Zadj_TS ~ df_M$age)
male_b <- male_reg$coefficients[[2]]
male_se <- coef(summary(male_reg))[2, 2]

female_reg <- lm(df_F$Zadj_TS ~ df_F$age)
female_b <- female_reg$coefficients[[2]]
female_se <- coef(summary(female_reg))[2, 2]

### Show results
print("No pregnant women")
pdiff <- 2*pnorm(-abs((male_b-female_b)/sqrt(male_se**2+female_se**2)),mean = 0, sd = 1)
print(paste("Male b:", round(male_b,3), "male se:", round(male_se,4),
            "Female b:", round(female_b,3), "female se:", round(female_se,4),
            "P_diff:", pdiff), sep=" ")

rm(male_reg,male_b,male_se,female_reg,female_b,female_se)


################################################
### Add time since last birth ##################

df$time_since_last_birth <- df$age - df$age_last_birth
df$time_since_last_birth_scaled <- as.vector(scale(df$time_since_last_birth))

model <- lm(df$TLc ~ df$time_since_last_birth_scaled)
print(lmp(model))
print(model)


################################################
### TL ~ age had kids vs never had Kids ########

df_F_kid <- df[which(df$ever_pregnant == 1),]
#df_F_kid <- df_F_kid[sample(1:141000,30000),]
df_F_nokid <- df[which(df$ever_pregnant == 0),]

### Scale TL
df_F_nokid$Zadj_TS <- as.vector((df_F_nokid$Zadj_TS))
df_F_kid$Zadj_TS <- as.vector((df_F_kid$Zadj_TS))

### Compare models
nokid_reg <- lm(df_F_nokid$Zadj_TS ~ df_F_nokid$age)
nokid_b <- nokid_reg$coefficients[[2]]
nokid_se <- coef(summary(nokid_reg))[2, 2]

kid_reg <- lm(df_F_kid$Zadj_TS ~ df_F_kid$age)
kid_b <- kid_reg$coefficients[[2]]
kid_se <- coef(summary(kid_reg))[2, 2]

### Show results
print("No kids vs kids")
pdiff <- 2*pnorm(-abs((nokid_b-kid_b)/sqrt(nokid_se**2+kid_se**2)),mean = 0, sd = 1)
print(paste("No kid b:", round(nokid_b,4), "No kid se:", round(nokid_se,4),
            "Kid b:", round(kid_b,4), "kid se:", round(kid_se,4),
            "P_diff:", pdiff), sep=" ")

# ### Groups same as ever pregnant
# # Boxplot
# ggplot(df, aes(x = factor(ever_pregnant), y = TLc)) +
#   geom_boxplot() +
#   labs(x = "Had Children", y = "Telomere Length")
# 
# # Then, you could perform a t-test to see if the difference in means is statistically significant
# df$TLc <- scale(df$TLc)
# t.test(TLc ~ ever_pregnant, data = df)

rm(nokid_reg,nokid_b,nokid_se,kid_reg,kid_b,kid_se)

################################################
### Select women around menopause ##############

# Plot distribution of menopausal age
hist(df$menopause, breaks = 43)
mean(df_F$menopause, na.rm = T)

# Filter women between 48-52 with menopausal information
df_age_50 <- df %>%
  dplyr::filter(age <= 52 & age >= 48) %>%
  dplyr::filter(sex == 2 & !is.na(had_menopause))

# df_age_50 <- dplyr::filter(df, sex == 2 & !is.na(had_menopause))

# Boxplot
ggplot(df_age_50, aes(x = factor(had_menopause), y = TLc)) +
  geom_boxplot() +
  labs(x = "Had menopause", y = "Telomere Length")

# T-test
df_age_50$TLc <- scale(df_age_50$TLc)
t.test(TLc ~ had_menopause, data = df_age_50)

### Separate pre and post menopausal women
df_post_menop <- dplyr::filter(df_age_50, had_menopause == 2)
df_pre_menop <- dplyr::filter(df_age_50, had_menopause == 1)

### Scale TL
df_post_menop$Zadj_TS <- as.vector((df_post_menop$Zadj_TS))
df_pre_menop$Zadj_TS <- as.vector((df_pre_menop$Zadj_TS))

### Compare models
pre_reg <- lm(df_pre_menop$Zadj_TS ~ df_pre_menop$age)
pre_b <- pre_reg$coefficients[[2]]
pre_se <- coef(summary(pre_reg))[2, 2]

post_reg <- lm(df_post_menop$Zadj_TS ~ df_post_menop$age)
post_b <- post_reg$coefficients[[2]]
post_se <- coef(summary(post_reg))[2, 2]

### Show results
print("Premonpausal vs postmenopausal")
pdiff <- 2*pnorm(-abs((pre_b-post_b)/sqrt(pre_se**2+post_se**2)),mean = 0, sd = 1)
print(paste("Premenopausal b:", round(pre_b,4), "Premenopausal se:", round(pre_se,4),
            "Postmenopausal b:", round(post_b,4), "Postmenopausal se:", round(post_se,4),
            "P_diff:", pdiff), sep=" ")

################################################
### Life periods in women ######################

# Re-select only women

df_F <- dplyr::filter(df, sex == 2)

# Generate the groups
df_F_grouped <- df_F %>%
  dplyr::filter(!is.na(had_menopause) | (had_menopause == 2 & is.na(nb_birth))) %>%
  mutate(group = case_when(
    nb_birth == 0 & had_menopause == 1 ~ "No kids and premenopausal",
    nb_birth > 0 & had_menopause == 1 ~ "Kids and premenopausal",
    had_menopause == 2  ~ "Postmenopausal" #, & nb_birth == 0
    #had_menopause == 2 & nb_birth > 0 ~ "Postmenopausal kids"
  ))

# Number of women per group
table(df_F_grouped$group)
df_F_grouped$TLc <- as.vector(scale(df_F_grouped$TLc))

# Perform the ANOVA
anova_result <- aov(TLc ~ group, data = df_F_grouped)
summary(anova_result)

# Perform the post-hoc test
posthoc_result <- TukeyHSD(anova_result)
print(posthoc_result)

# Form groups
df_group_1 <- df_F_grouped[which(df_F_grouped$group == "No kids and premenopausal"),]
df_group_2 <- df_F_grouped[which(df_F_grouped$group == "Kids and premenopausal"),]
df_group_3 <- df_F_grouped[which(df_F_grouped$group == "Postmenopausal"),]

### Scale TL
df_group_1$Zadj_TS <- scale(as.vector((df_group_1$Zadj_TS)))
df_group_2$Zadj_TS <- scale(as.vector((df_group_2$Zadj_TS)))
df_group_3$Zadj_TS <- scale(as.vector((df_group_3$Zadj_TS)))

### Compare models
group_1_reg <- lm(df_group_1$Zadj_TS ~ df_group_1$age)
group_1_b <- group_1_reg$coefficients[[2]]
group_1_se <- coef(summary(group_1_reg))[2, 2]

group_2_reg <- lm(df_group_2$Zadj_TS ~ df_group_2$age)
group_2_b <- group_2_reg$coefficients[[2]]
group_2_se <- coef(summary(group_2_reg))[2, 2]

group_3_reg <- lm(df_group_3$Zadj_TS ~ df_group_3$age)
group_3_b <- group_3_reg$coefficients[[2]]
group_3_se <- coef(summary(group_3_reg))[2, 2]

### Show results
print("Group_1 vs Group_2")
pdiff <- 2*pnorm(-abs((group_1_b-group_2_b)/sqrt(group_1_se**2+group_2_se**2)),mean = 0, sd = 1)
print(paste("Group_1 b:", round(group_1_b,4), "Group_1 se:", round(group_1_se,4),
            "Group_2 b:", round(group_2_b,4), "Group_2 se:", round(group_2_se,4),
            "P_diff:", pdiff), sep=" ")

print("Group_1 vs Group_3")
pdiff <- 2*pnorm(-abs((group_1_b-group_3_b)/sqrt(group_1_se**2+group_3_se**2)),mean = 0, sd = 1)
print(paste("Group_1 b:", round(group_1_b,4), "Group_1 se:", round(group_1_se,4),
            "Group_3 b:", round(group_3_b,4), "Group_3 se:", round(group_3_se,4),
            "P_diff:", pdiff), sep=" ")

print("Group_2 vs Group_3")
pdiff <- 2*pnorm(-abs((group_2_b-group_3_b)/sqrt(group_2_se**2+group_3_se**2)),mean = 0, sd = 1)
print(paste("Group_2 b:", round(group_2_b,4), "Group_2 se:", round(group_2_se,4),
            "Group_3 b:", round(group_3_b,4), "Group_3 se:", round(group_3_se,4),
            "P_diff:", pdiff), sep=" ")


################################################
### Age since menopause ########################

df_F_grouped$time_since_menopause <- df_F_grouped$age - df_F_grouped$menopause
df_F_grouped$time_since_menopause[which(df_F_grouped$had_menopause == 1)] <- 0 

hist(df_F_grouped$time_since_menopause)


### We define three periods
df_per <- df_F_grouped %>% 
  dplyr::filter(is.na(age_first_birth) | is.na(menopause) | (age_first_birth <= menopause)) %>% # Remove women with postmenopausal kids
  dplyr::filter(!is.na(had_menopause)) %>% #                              Remove women with NA in whether they had menopause or not
  dplyr::filter(!(had_menopause == 2 & is.na(menopause))) %>% #           Remove women with menopause but with NA menopausal age
  dplyr::filter(!is.na(ever_pregnant)) %>% #                              Remove women with NA in whether they were ever pregnant
  dplyr::filter(!(ever_pregnant == 1 & is.na(age_first_birth))) %>% #     Remove women that were pregnant without age of first birth
  mutate(A = ifelse(ever_pregnant == 0, #                                 Period A = pre menopausal time until first kid
                    ifelse(had_menopause == 1, age, menopause), age_first_birth)) %>% 
  mutate(B = ifelse(ever_pregnant == 0, 0, #                              Period B = pre menopausal time after first kid
                    ifelse(had_menopause == 2,menopause - age_first_birth, age - age_first_birth))) %>%
  mutate(C = ifelse(had_menopause == 1, 0, age - menopause)) %>% #        Period C = post menopausal time
  select(c("Zadj_TS","TLc","age","nb_birth", "ever_pregnant","age_first_birth","had_menopause","menopause","group","A","B","C")) %>%
  mutate(agecheck = A + B + C)

model <- (lm(Zadj_TS ~ A + B + C + had_menopause + ever_pregnant, df_per))
summary(model)

# Distribution plot
ggplot(df_per, aes(x=A, fill="A")) + 
  geom_histogram(alpha=0.5, binwidth=1, color="black") +
  geom_histogram(data=df_per, aes(x=B, fill="B"), alpha=0.5, binwidth=1, color="black") +
  geom_histogram(data=df_per, aes(x=C, fill="C"), alpha=0.5, binwidth=1, color="black") +
  scale_fill_manual(values=c("red", "green", "blue")) + 
  labs(
    fill="Variable",
    title="Overlapped Histogram of A, B, and C",
    x="Years spent in period",
    y=""
  ) +
  theme_minimal()

# Extract coefficients and standard errors
coef_vals <- coef(model)
std_errs <- coef(summary(model))[, 2]

# Define a function to compute the p-value for a pair of variables
compute_p_value <- function(b1, b2, se1, se2) {
  z_value <- abs(b1 - b2) / sqrt(se1^2 + se2^2)
  p_value <- 2 * pnorm(-abs(z_value), mean = 0, sd = 1)
  return(p_value)
}

# Compute p-values for the pairs
p_value_A_vs_B <- compute_p_value(coef_vals["A"], coef_vals["B"], std_errs["A"], std_errs["B"])
p_value_A_vs_C <- compute_p_value(coef_vals["A"], coef_vals["C"], std_errs["A"], std_errs["C"])
p_value_B_vs_C <- compute_p_value(coef_vals["B"], coef_vals["C"], std_errs["B"], std_errs["C"])

# Print the p-values
cat("p-value for A vs B:", p_value_A_vs_B, "\n")
cat("p-value for A vs C:", p_value_A_vs_C, "\n")
cat("p-value for B vs C:", p_value_B_vs_C, "\n")

### Control for non linear effects
df_per$age2 <- df_per$age**2
summary(lm(df_per$Zadj_TS ~ df_per$age + df_per$age2))

df_group_1$age2 <- df_group_1$age**2
summary(lm(df_group_1$Zadj_TS ~ df_group_1$age + df_group_1$age2))

df_group_2$age2 <- df_group_2$age**2
summary(lm(df_group_2$Zadj_TS ~ df_group_2$age + df_group_2$age2))

df_group_3$age2 <- df_group_3$age**2
summary(lm(df_group_3$Zadj_TS ~ df_group_3$age + df_group_3$age2))


################################################
### Age since menopause ########################

# Linear model for females
model1 <- lm(Zadj_TS ~ A + B + C, data = df_per)
coeff1 <- summary(model1)$coefficients

# model1 <- lm(Zadj_TS ~ A + I(A^2) + B + I(B^2)+ C + I(C^2), data = df_per)
# coeff1 <- summary(model1)$coefficients

# Adjusted Linear model for age vs. Zadj_TS in males to include age squared
model2 <- lm(Zadj_TS ~ age + I(age^2), data = df_M)
coeff2 <- summary(model2)$coefficients

# Linear model for female
model3 <- lm(Zadj_TS ~ age + I(age^2), data = df_per)
coeff3 <- summary(model3)$coefficients

# Extract coefficients and standard errors from model1
Intercept <- coeff1["(Intercept)", "Estimate"]
A <- coeff1["A", "Estimate"]
B <- coeff1["B", "Estimate"]
C <- coeff1["C", "Estimate"]

# Create a new data frame for prediction
data <- data.frame(age = 0:70) 
# data <- data.frame(age = c(0:26, 26:50, 50:70))

# Get predictions with confidence intervals
preds <- predict(model2, newdata = data, interval = "confidence")

# Bind the predictions to the new data
data <- cbind(data, preds)

# Get predictions with confidence intervals for model3
preds3 <- predict(model3, newdata = data, interval = "confidence")

# Bind the predictions for model3 to the new data
data$fit3 <- preds3[,1]
data$lwr3 <- preds3[,2]
data$upr3 <- preds3[,3]

# Calculate predictions and confidence intervals for each segment
years1 <- 0:26
years2 <- 27:50
years3 <- 51:70

# Calculate predictions and confidence intervals for each segment
preds1 <- predict(model1, newdata = data.frame(A = years1, B = 0, C = 0, 
                                               had_menopause = 0, ever_pregnant = 0), interval = "confidence")
preds2 <- predict(model1, newdata = data.frame(A = 26, B = years2 - 26, C = 0, 
                                               had_menopause = 0, ever_pregnant = 1), interval = "confidence")
preds3 <- predict(model1, newdata = data.frame(A = 26, B = 24, C = years3 - 50, 
                                               had_menopause = 1, ever_pregnant = 1), interval = "confidence")

# Combine the predictions and confidence intervals
LTL_fit <- c(preds1[,1], preds2[,1], preds3[,1])
LTL_lwr <- c(preds1[,2], preds2[,2], preds3[,2])
LTL_upr <- c(preds1[,3], preds2[,3], preds3[,3])

data$LTL_fit <- LTL_fit
data$LTL_lwr <- LTL_lwr
data$LTL_upr <- LTL_upr

p <- ggplot(data, aes(x = age)) +
  # Adding background highlight
  geom_rect(aes(xmin = 40, xmax = 70, ymin = -Inf, ymax = Inf), fill = "#fefae0", inherit.aes = FALSE) +
  # Adding vertical lines
  geom_vline(aes(xintercept = 26), linetype="dotted", size = 1, color = "#6c757d") +
  geom_vline(aes(xintercept = 50), linetype="dotted", size = 1, color = "#6c757d") +
  # Plotting for model2
  geom_line(aes(y = fit, color = "Male"), size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = "Male"), alpha = 0.2) +
  
  geom_line(aes(y = LTL_fit, color = "Female"), size = 1) +
  geom_ribbon(aes(ymin = LTL_lwr, ymax = LTL_upr, fill = "Female"), alpha = 0.2) +
  
  # Adding labels
  annotate("text", x = 13, y = 0.55, label = sprintf("bold('\u03B2 = %0.3f')", A), parse = TRUE, color = "#C0392B") +
  annotate("text", x = 38, y = 0.5, label = sprintf("bold('\u03B2 = %s')", sprintf("%.3f", B)), parse = TRUE, vjust = -0.5, color = "#C0392B") +
  annotate("text", x = 60, y = 0.12, label = sprintf("bold('\u03B2 = %s')", round(C, 3)), parse = TRUE, vjust = -0.5, color = "#C0392B") +
  
  annotate("text", x = 15, y = 1.15, 
           label = sprintf("bold(beta[age] == %s)", 
                           round(coeff2["age", "Estimate"], 3)), 
           parse = TRUE, vjust = -0.5, color = "#2471A3") +
  
  annotate("text", x = 15, y = 1.05, 
           label = sprintf("bold(beta[age^2] == %s)", 
                           round(coeff2["I(age^2)", "Estimate"], 5)), 
           parse = TRUE, vjust = -0.5, color = "#2471A3") +
  
  geom_label(aes(x = 26, y = min(lwr)+0.1), label = "Age of first birth", vjust = 1.5, hjust = 0.5, fill = "white", label.size = NA) +
  geom_label(aes(x = 50, y = min(lwr)+0.1), label = "Age at menopause", vjust = 1.5, hjust = 0.5, fill = "#fefae0", label.size = NA) +
  # Axis
  xlab("Age") +
  ylab("Predicted LTL") +
  scale_x_continuous(breaks = seq(0, 70, 10)) +
  # Legend
  scale_color_manual(values = c("Male" = "#2471A3", "Female" = "#C0392B"), name = "") +
  guides(fill = "none") + # Takes legend of geom_point away
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=12),
        #axis.text.y = element_blank(),
        axis.title = element_text(size=14),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = c(0.15, 0.3)
  )

print(p)


ggsave(file="reproductive_periods.svg", 
       plot = p, width=7, height=4)
ggsave(file="reproductive_periods.png", 
       plot = p, width = 7, height = 4, dpi = 600)


################################################
### Plot mock for Figure 1 #####################

p_mock <- ggplot(data, aes(x = age)) +
  geom_line(aes(y = fit, color = "Male"), size = 2) + # Adjusted size for thicker line
  geom_line(aes(y = LTL_fit, color = "Female"), size = 2) + # Adjusted size for thicker line
  xlab("Age") +
  ylab("LTL") +
  scale_x_continuous(breaks = seq(0, max(data$age, na.rm = TRUE), 20)) + # Reverted to 10-year intervals
  scale_color_manual(values = c("Male" = "#2471A3", "Female" = "#C0392B"), name = "") +
  geom_vline(aes(xintercept = 26), linetype="dotted", size = 1, color = "#6c757d") +
  geom_vline(aes(xintercept = 50), linetype="dotted", size = 1, color = "#6c757d") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=12, face="bold"), # Made axis text bold
    axis.title = element_text(size=14, face="bold"), # Made axis title bold
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "none" 
  )

print(p_mock)
ggsave(file="mock_reproductive_periods.svg", 
       plot = p_mock, width=2.5, height=2)