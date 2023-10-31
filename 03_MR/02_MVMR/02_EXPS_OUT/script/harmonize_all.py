### Import modules
print("Harmonize all script start")
import pandas as pd 

### Load outcome dataframe
snps_df = pd.read_csv(snakemake.input[0], sep = '\t', usecols = ["SNP"])
print("First dataframe loaded", snps_df.head(3))
print("Number of first trait SNPs:", len(snps_df.index))

### Loop over files to get SNPs present in all files
for file in snakemake.input[1:len(snakemake.input)]:
    snps_med = pd.read_csv(file, sep = '\t', usecols = ["SNP"])
    snps_df = snps_med[snps_med.SNP.isin(snps_df.SNP)]

print("Number of harmonized SNPs:", len(snps_df.index))
### Save list of SNPs
snps_df.to_csv(snakemake.output[0], index = False, header = False)