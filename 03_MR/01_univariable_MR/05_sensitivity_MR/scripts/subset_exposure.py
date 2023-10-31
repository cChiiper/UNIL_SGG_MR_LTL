import pandas as pd

# Subset exposure data
expo_df = pd.read_table(snakemake.input["expo_df"], header=0)
print("Exposure df is in memory.")

snps = pd.read_table(snakemake.input["snps"], header=None)
snps.columns = ['SNP']

expo_df = expo_df[expo_df['SNP'].isin(snps['SNP'])]
print("Exposure df is subsetted to clumped SNPs.")

############################################
### Save file ##############################
expo_df[['SNP', 'A1', 'A2', 'Freq', 'b', 'se', 'p', 'N']].to_csv(snakemake.output[0], index = False, header = True,  sep='\t')