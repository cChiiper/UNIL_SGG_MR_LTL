## 03_MR

### 00_uk10k_snps_filter.sh
Bash script to extract SNPs of specific region from the uk10k reference panel

### 00_GWAS_prep
Contains scripts to format GWAS summary statistics into the same format (https://yanglab.westlake.edu.cn/software/gcta/#COJO).
* GWAS_prep_composed: Scripts to generate the reproductive lifespan GWAS.
* GWAS_prep_Neale: Scripts to format GWAS from the Neale lab (http://www.nealelab.is/uk-biobank).
* GWAS_prep_PAN_UKB: Scripts to format GWAS from the Pan-UK Biobank ([link](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288)).
* GWAS_prep_special: Scripts to format particular GWAS summary statistics.

More information can be found at https://github.com/masadler/smrivw/wiki/Univariable-MR

### 01_univariable_MR
Contains scripts to run univariable two-sample MR.
* 01_MR_TSMR: Scripts to run univariable MR [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/) version 0.5.7 package.
* 02_MR_APSS: Scripts to run [MR-APSS] (https://github.com/YangLabHKUST/MR-APSS)
* 03_MR_PRESSO: Scripts that run [MR-PRESSO](https://github.com/rondolab/MR-PRESSO) with the TwoSampleMR wrapper.
* 04_sensitivity_IV: Scripts to select instrumental variables (IVs) based on a stringent Steiger filter.
* 05_sensitivtiy_MR: Scripts to run univariable MR with the IVs from 04_sensitivity_IV.

### 02_MVMR
Contains scripts to run multivariable MR (MVMR).
* 01_EXP_MED_OUT: Scripts to run mediation analysis with set exposure, mediator, and outcome.
* 02_EXPS_OUT: Scripts to run MVMR with multiple exposures on one outcome.