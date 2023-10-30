################################################################################
### Generate reproductive lifespan trait from Menarche and Menopause         ###
### Author: Samuel Moix                                                      ###
### Date: 02.08.2022                                                         ###
################################################################################

################################################
### Libraries ##################################
require(data.table)
library(dplyr)

################################################
### Working directories ########################
data_folder = "/SET/PATH/TO/DIRECTORY/"
export_folder = "/SET/PATH/TO/DIRECTORY/"


################################################
### Load data ##################################
df_MNC <- as.data.frame(data.table::fread(file = file.path(data_folder,  
                                             "MENARCHE_gwas_summary_uk10kck.ma"), 
                            sep = '\t', header = TRUE))

df_MNP <- as.data.frame(data.table::fread(file = file.path(data_folder,  
                                                           "MENOPAUSE_gwas_summary_uk10kck.ma"), 
                                          sep = '\t', header = TRUE))


################################################
### Calculate new trait data ###################

### Merge needed values
df <- merge(df_MNP[,which(colnames(df_MNP) != "p")],df_MNC[,c("SNP","b","se", "N")],by="SNP", all.x = FALSE)
rm(df_MNP)
rm(df_MNC)

### .x for menopause and .y for menarche
# SD = 1 because of effect size standardization
sdmena <- 1.628354
sdmeno <- 5.433486
df <- df %>% 
  mutate(b = sdmeno*b.x - sdmena*b.y) %>%
  mutate(se = sqrt(sdmeno**2*se.x**2 + sdmena**2*se.y**2)) %>%
  mutate(p = 2*pnorm(-abs((b/se)), mean = 0, sd = 1)) %>%
  mutate(N = min(N.x,N.y)) %>%
  select(c("SNP","A1","A2","Freq","b","se","p","N"))

### Export file
fwrite(df, file.path(export_folder,  "REPLP_gwas_summary_uk10kck.ma"), col.names = T, row.names = F, quote = F, sep = "\t")

