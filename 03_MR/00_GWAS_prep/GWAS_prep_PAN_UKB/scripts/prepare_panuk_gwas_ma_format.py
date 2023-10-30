import numpy as np
import pandas as pd

pheno = snakemake.wildcards["trait"]
ncases = snakemake.params["ncases"]
ncontrols = snakemake.params["ncontrols"]

variant_df = pd.read_csv(snakemake.input["variant_file"], compression='gzip', sep = '\t', usecols = ["af_EUR", "rsid", "chrom", "pos", "ref", "alt"])
gwas_df = pd.read_csv(snakemake.params["path"], compression='gzip', sep = '\t', usecols = ["chr", "pos", "ref", "alt", "beta_EUR", "se_EUR", "pval_EUR"])

gwas_df["SNP"] = gwas_df["chr"].astype(str) + gwas_df["pos"].astype(str)
variant_df["SNP"] = variant_df["chrom"].astype(str) + variant_df["pos"].astype(str)

gwas_df = gwas_df.merge(variant_df, on = 'SNP', suffixes=("_g", "_var"))
del variant_df
print(gwas_df.head())

gwas_df["freq"] = gwas_df["af_EUR"]
mask = ((gwas_df["ref_g"] == gwas_df["alt_var"]) & (gwas_df["alt_g"] == gwas_df["ref_var"]))
gwas_df["freq"][mask] = 1-gwas_df["freq"]

gwas_df = gwas_df[~gwas_df.beta_EUR.isna()]
gwas_df = gwas_df[~np.isinf(gwas_df.beta_EUR)]

gwas_df['N'] = (4*ncases*ncontrols)/(ncases+ncontrols)
gwas_df['b'] = gwas_df['beta_EUR']/gwas_df['se_EUR']/np.sqrt(gwas_df['N'])
gwas_df['se'] = 1/np.sqrt(gwas_df['N'])
gwas_df['p'] = gwas_df['pval_EUR'].astype(np.float)
gwas_df['p'] = np.exp(gwas_df.pval_EUR)

# Ref:A2 and Alt:A1 because the effect size refers to alt and the frequency as well
gwas_df = gwas_df.drop(["SNP"], axis = 1)
gwas_df = gwas_df.rename(columns = {'ref_g': 'A2', 'alt_g': 'A1','rsid': 'SNP'})
gwas_filtered_df = gwas_df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].drop_duplicates()
gwas_filtered_df = gwas_filtered_df.dropna()

print('Number of entries before and after duplicate removing:', gwas_df.shape[0], gwas_filtered_df.shape[0])

gwas_filtered_df.to_csv(snakemake.output[0], index = False, sep = '\t')

