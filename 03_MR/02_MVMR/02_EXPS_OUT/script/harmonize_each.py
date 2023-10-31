### Import modules
print("Harmonize each script start")
import pandas as pd 

### Load outcome dataframe
snps_all = pd.read_csv(snakemake.input[0], sep = '\t', names = ["SNP"])
snps_each = pd.read_csv(snakemake.input[1], sep = '\t', names = ["SNP"])

snps_df = snps_each[snps_each.SNP.isin(snps_all.SNP)]
print("Number of harmonized SNPs:", len(snps_df.index))

### Save list of SNPs
snps_df.to_csv(snakemake.output[0], index = False, header = False)