import numpy as np
import pandas as pd

pheno = snakemake.wildcards["trait"]

variant_df = pd.read_csv('variants.tsv.gz', compression='gzip', sep = '\t')
gwas_df = pd.read_csv(snakemake.input[0], compression='gzip', sep = '\t')

gwas_df = gwas_df.merge(variant_df, on = 'variant')
gwas_df = gwas_df[~gwas_df.beta.isna()]

gwas_df['b'] = gwas_df['beta']/gwas_df['se']/np.sqrt(gwas_df['n_complete_samples'])
gwas_df['se'] = 1/np.sqrt(gwas_df['n_complete_samples'])

# Ref:A2 and Alt:A1 because the effect size refers to alt and the frequency as well
gwas_df = gwas_df.rename(columns = {'n_complete_samples':'N', 'pval':'p', 'ref': 'A2', 'alt': 'A1', 'AF':'freq','rsid': 'SNP'})
gwas_filtered_df = gwas_df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].drop_duplicates()
gwas_filtered_df = gwas_filtered_df.dropna(subset = ['b'])

print('Number of entries before and after duplicate removing:', gwas_df.shape[0], gwas_filtered_df.shape[0])

gwas_filtered_df.to_csv(snakemake.output[0], index = False, sep = '\t')

