import pandas as pd 

outcome_df = pd.read_csv(snakemake.input["outcome_df"], sep = '\t')
print("Outcome DF")
print(outcome_df.head())

med_df = pd.read_csv(snakemake.input["med_df"], sep = '\t')
print("Mediator DF")
print(med_df.head())

snps_df = pd.read_csv(snakemake.input["snps"], sep = '\t', names = ['SNP'])
print("Exposure SNPs")
print(snps_df.head())

snps_df = snps_df[snps_df.SNP.isin(outcome_df.SNP)]
print("subset 1")
print(snps_df.head())

snps_df = snps_df[snps_df.SNP.isin(med_df.SNP)]
print("subset 2")
print(snps_df.head())


snps_df.to_csv(snakemake.output[0], index = False, header = False)

