### Load packages
import numpy as np
import pandas as pd

### Load file GCST90091033
# http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90091001-GCST90092000/GCST90091033/

df = pd.read_csv(snakemake.input[0], compression='gzip', sep = '\t')
print(df.head(2))

### Rename columns accordingly
df = df.rename(columns = {'variant_id':'SNP','effect_allele':'A1','other_allele':'A2','effect_allele_frequency':'freq','p_value':'p'})

# Set N from (4*ncases*ncontrols)/(ncases+ncontrols)
df['N'] = (4*8434*770180)/(8434+770180)

### Effect size standardization
df['b'] = df['beta']/df['standard_error']/np.sqrt(df['N'])
df['se'] = 1/np.sqrt(df['N'])

### Export to .ma
df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(snakemake.output[0], index = False, sep = '\t')

