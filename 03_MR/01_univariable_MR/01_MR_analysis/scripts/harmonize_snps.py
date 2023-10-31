import pandas as pd 

df = pd.read_csv(snakemake.input[1], sep = '\t')
snps_df = pd.read_csv(snakemake.input[0], sep = '\t', names = ['SNP'])

snps_df = snps_df[snps_df.SNP.isin(df.SNP)]

snps_df.to_csv(snakemake.output[0], index = False, header = False)

