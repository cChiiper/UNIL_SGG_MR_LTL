import numpy as np
import pandas as pd

pheno = snakemake.wildcards["trait"]

df = pd.read_csv(snakemake.input[0], compression='gzip', sep = '\t')

print(df.head(2))

df = df.rename(columns = {'variant_id':'SNP', 'effec_allele': 'A1', 'other_allele': 'A2', 'effect_allele_frequency':'freq', 'p_value':'p'})

df['N'] = 472174

df['b'] = df['beta']/df['standard_error']/np.sqrt(df['N'])
df['se'] = 1/np.sqrt(df['N'])

df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(snakemake.output[0], index = False, sep = '\t')