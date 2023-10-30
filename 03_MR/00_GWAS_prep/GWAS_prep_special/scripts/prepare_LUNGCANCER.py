### Load packages
import numpy as np
import pandas as pd

### Load file GCST90041821
# http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90041001-GCST90042000/GCST90041821/ 
pheno = snakemake.wildcards["trait"]
df = pd.read_csv(snakemake.input[0], compression='gzip', sep = '\t')
print(df.head(2))

### Rename columns accordingly
df = df.rename(columns = {'variant_id':'SNP','effect_allele':'A1','other_allele':'A2','effect_allele_frequency':'freq','p_value':'p'})

### Effect size standardization
df['b'] = df['T']/df['SE_T']/np.sqrt(df['N'])
df['se'] = 1/np.sqrt(df['N'])

### Export to .ma
df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(snakemake.output[0], index = False, sep = '\t')