### Load packages
import numpy as np
import pandas as pd
import pandas as pd
import fastparquet 

### Load file
df = pd.read_parquet(snakemake.input[0])
# snp  pos_name chromosome  position a1 a2  beta  se   z pval af1  n_iids
print(df.head(2))

### Rename columns accordingly
df = df.rename(columns = {'snp':'SNP','a1':'A1', 'a2':'A2', 'af1':'freq', 'pval':'p', 'n_iids':'N'}) 

### Effect size standardization
df['b'] = df['beta']/df['se']/np.sqrt(df['N'])
df['se'] = 1/np.sqrt(df['N'])

### Export to .ma
df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(snakemake.output[0], index = False, sep = '\t')