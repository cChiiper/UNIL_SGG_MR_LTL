#################################################################################################################
### Plotting of Telomere length with age, sex and cancer                                                      ###
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
data_folder = "/SET/PATH/TO/DIRECTORY"

################################################
### Parameters  ################################

### Whether to remove patients with blood cancer
par_rm_blood_cancer <- TRUE

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
### Add corrected phenotypes ###################

df$TLc_age = residuals(lm(df$Zadj_TS ~ df$age))
df$TLc_age_sex <- residuals(lm(df$Zadj_TS ~ df$age * df$sex))


################################################
### Plotting ###################################

### Boxplot TL ~ age regression line ########### ########### ########### ########### ########### ###########

lm_plot_age <- ggplot(df, aes(x = age, y = Zadj_TS)) +
  geom_point(alpha = 0) + # Has to be set to fit regression line alpha = 0  to not show points
  geom_boxplot(data = df, aes(x = age, y = Zadj_TS, group = factor(age)), 
               alpha = 0.7, lwd = 0.7, outlier.alpha = 0.05, fill = "#E3E4B4") + 
  labs(title = "Telomere length vs age",  # Editing labs
       x = "Age [years]", 
       y = "Z-adjusted T/S log [-]",
       fill = "Sex") +
  scale_x_continuous(breaks = seq(40, 70, by = 2))+ # Editing x scale
  geom_smooth(method = "lm", aes(x = age, y = Zadj_TS), formula = y ~ x, color = "#CB0000") + # Draw lm
  stat_regline_equation( # Add lm equation
    aes(x = age, y = Zadj_TS, label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    formula = y ~ x, show.legend = FALSE, color = "#CB0000") +
  theme(legend.position="left",  # Change layout
        plot.title = element_text(size=20),
        axis.text = element_text(size = 12),
        axis.title=element_text(size=14),
        legend.key.height= unit(1, 'cm'), 
        legend.key.width= unit(1, 'cm'), 
        legend.title = element_text(face = "bold") 
  ) +
  ylim(-5,5)

ggMarginal(lm_plot_age, margins = "both", type = "density", size = 10, fill = "#E3E4B4")


# plot(df$age, df$TLc_age, ylab="Residuals", xlab="Age", main="TL ~ age") 
# abline(0, 0)  
  
### Boxplot TL ~ age by sex regression line ########### ########### ########### ########### ########### ####

lm_plot_age <- ggplot(df, aes(x = age, y = Zadj_TS, group = factor(sex), color = factor(sex))) +
  geom_point(alpha = 0) + # Has to be set to fit regression line alpha = 0  to not show points
  geom_boxplot(data = df, aes(x = age, y = Zadj_TS, group = interaction(age, factor(sex)), fill = factor(sex)), 
               alpha = 0.7, lwd = 0.7, outlier.alpha = 0.05) +
  labs(title = "",  # Editing labs #Telomere length vs age by sex
       x = "Age [years]", 
       y = "LTL", # Z-adjusted T/S log [-]
       fill = "Sex") +
  scale_x_continuous(breaks = seq(40, 70, by = 5))+ # Editing x scale
  geom_smooth(method = "lm", aes(x = age, y = Zadj_TS, group=sex, color=sex), formula = y ~ x) + # Draw lm
  stat_regline_equation( # Add lm equation
    aes(x = age, y = Zadj_TS, color = sex, label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    formula = y ~ x, show.legend = FALSE) +
  scale_fill_manual(values=c("#94CDE1", "#EA6262"), labels = c("males", "females")) + #Set colors
  scale_color_manual(values=c("#2471A3","#C0392B")) +
  guides(colour = "none") + # Takes legend of geom_point away
  theme(legend.position="bottom",  # Change layout
        legend.direction = "horizontal",
        plot.title = element_text(size=18),
        axis.text = element_text(size = 14),
        axis.title=element_text(size=16),
        legend.key.height= unit(1.2, 'cm'), 
        legend.key.width= unit(1, 'cm'), 
        legend.text = element_text(size=14),
        legend.title = element_blank() 
  ) +
  ylim(-5,5)

# Add marginal plots (distribution)
sex_plot <- ggMarginal(lm_plot_age, margins = "both", type = "density", groupColour = T, groupFill = T, size = 10)

require(svglite)
ggsave(file="tl_age_sex_reg.svg", plot=sex_plot, width=8, height=7)

################################################
### Statistics #################################

### Function to calculate p-value
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

### TL ~ age all regression ####################
TL_age_model = lm(df$Zadj_TS ~ df$age)
summary(TL_age_model)
lmp(TL_age_model)

### TL ~ scale(sex)
df$sex_scale <- scale(as.numeric(df$sex))
TL_ss_model = lm(df$Zadj_TS ~ df$sex_scale)
summary(TL_ss_model)
lmp(TL_ss_model)

# residu <- residuals(TL_age_model)
# par(mfrow=c(1,3))
# plot(df$age, residu) ; abline(lm(residu~df$age))
# plot(df$age, abs(residu)) ; abline(lm(abs(residu)~df$age))
# qqnorm(residu) ; qqline(residu)
# cor.test(abs(residu), df$age, method="spearman") # p-value = 0.03148 but already log transformed

par(mfrow=c(2,2)) 
plot(TL_age_model)

### TL ~ age Male Female regression ############

# df$age = scale(df$age)
# df_M <- df[which(df$sex == 1),]
# df_F <- df[which(df$sex == 2),]
# df_M$Zadj_TS <- scale(df_M$Zadj_TS)
# df_F$Zadj_TS <- scale(df_F$Zadj_TS)
# 
# lm_male = lm(df_M$Zadj_TS ~ df_M$age)
# lm_female = lm(df_F$Zadj_TS ~ df_F$age)


lm_male = lm(df$Zadj_TS[which(df$sex == 1)] ~ df$age[which(df$sex == 1)])
lm_female = lm(df$Zadj_TS[which(df$sex == 2)] ~ df$age[which(df$sex == 2)])

print(summary(lm_male))
print(summary(lm_female))

print("Pdiff males vs females:")
2*pnorm(-abs((lm_male$coefficients[[2]]-lm_female$coefficients[[2]])/sqrt(summary(lm_male)$coefficients[4][1]**2+summary(lm_female)$coefficients[4][1]**2)),mean = 0, sd = 1)

### Pearson's cor test between age and TL #####
cor.test(df$age, df$Zadj_TS)

### Average TL between male and female
ggplot(df, aes(x=as.factor(sex), y=Zadj_TS)) + 
  geom_boxplot(fill="slateblue", alpha=0.2,outlier.shape = NA) + 
  coord_cartesian(ylim = c(-3,3)) +
  xlab("Sex") +
  ylab("TL") +
  scale_x_discrete(labels = c('M','F'))


t.test(Zadj_TS ~ sex, data = df) 