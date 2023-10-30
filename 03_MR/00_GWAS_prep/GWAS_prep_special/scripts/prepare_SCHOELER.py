### Load packages
import numpy as np
import pandas as pd

### Load file
pheno = snakemake.wildcards["trait"]
df = pd.read_csv(snakemake.input[0], sep = '\t')
print(df.head(2))

### Rename columns accordingly
df = df.rename(columns = {'MAF':'freq', 'P':'p'})

### Effect size standardization
df['b'] = df['BETA']/df['SE']/np.sqrt(df['N'])
df['se'] = 1/np.sqrt(df['N'])

### Export to .ma
df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(snakemake.output[0], index = False, sep = '\t')