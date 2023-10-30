### Load packages
import numpy as np
import pandas as pd

### Load file
pheno = snakemake.wildcards["trait"]
df = pd.read_csv(snakemake.input[0], compression='zip', sep = '\t')
print(df.head(2))

### Rename columns accordingly
df = df.rename(columns = {'ALLELE1':'A1', 'ALLELE0':'A2', 'A1FREQ':'freq', 'P_BOLT_LMM':'p'}) # P_inf is shrinking things more

### Add N based on fraction of individuals with missing genotype (F_MISS) 
df['N'] = (205327 - 205327*df['F_MISS']).apply(np.floor).astype(int)

### Changing beta of X chromosome (23) as it is haploid for males
df.loc[df["CHR"] == 23, "BETA"] = df["BETA"]/2
### Effect size standardization
df['b'] = df['BETA']/df['SE']/np.sqrt(df['N'])
df['se'] = 1/np.sqrt(df['N'])

### Export to .ma
df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(snakemake.output[0], index = False, sep = '\t')