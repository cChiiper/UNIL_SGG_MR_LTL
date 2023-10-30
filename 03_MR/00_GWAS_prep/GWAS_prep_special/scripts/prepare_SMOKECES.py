### Load packages
import numpy as np
import pandas as pd

### Load file
# https://conservancy.umn.edu/handle/11299/201564 
pheno = snakemake.wildcards["trait"]
df = pd.read_csv(snakemake.input[0], compression='gzip', sep = '\t')
print(df.head(2))

### Rename columns accordingly
df = df.rename(columns = {'RSID':'SNP', 'ALT': 'A1', 'REF': 'A2', 'AF':'freq','PVALUE':'p'})

### Effect size standardization
df['b'] = df['BETA']/df['SE']/np.sqrt(df['N'])
df['se'] = 1/np.sqrt(df['N'])

### Export to .ma
df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(snakemake.output[0], index = False, sep = '\t')