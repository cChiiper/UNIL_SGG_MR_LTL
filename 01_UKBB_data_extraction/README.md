## 01_UKBB_data_extraction

### 00_extract_UKBB_metadata.ipynb
Python script to extract metadata from the UK Biobank portal (https://biobank.ndph.ox.ac.uk/ukb/index.cgi)
Takes a metadata dataframe containing a column named "FieldID" corresponding to UK Biobank data-fields (e.g., 2724 for Had menopause). Primarily used to get trait descriptions.

### 01_sample_filtering.R
R script to filter individuals of interest i.e participants with known sex, age, and LTL after the exclusion of individuals of non-white and non-British ancestry (self-reported + genetically defined), relatives, and gender mismatches (see UKBB Resource 531), as well as those who retracted their participation.

### 02_trait_extraction.R
R script to extract UK Biobank traits based on data-fields
Takes a metadata dataframe containing "FieldID" for UK Biobank data-fields, "pheno" for abbreviation, "type" for the type of trait (factor, integer, continuous). 

### 03_disease_extraction.R
R script to extract disease traits based on an inclusion and exclusion list. See https://doi.org/10.1101/2023.07.31.23293408.

### 04_merge_dataframes.R
R script to merge non-disease and disease traits together.

### 05_data_formatting.R
R script to adjust the format of certain traits, add composite traits, and adjust some variables to have a clean dataframe for observational analyses.