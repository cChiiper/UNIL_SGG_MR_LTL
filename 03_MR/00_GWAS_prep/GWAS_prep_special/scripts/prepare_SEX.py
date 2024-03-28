### Load packages
import numpy as np
import pandas as pd

### Load file GCST90013474
# http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013474/
# GWAS of sex in UK biobank. Coding: 1= Female, 0 = Male.

df = pd.read_csv(snakemake.input[0], compression='gzip', sep = '\t')
print(df.head(2))

### Rename columns accordingly
df = df.rename(columns = {'variant_id':'SNP','effect_allele':'A1','other_allele':'A2','effect_allele_frequency':'freq','p_value':'p'})

# Set N as not present
df['N'] = 452302

### Effect size standardization
df['b'] = df['beta']/df['standard_error']/np.sqrt(df['N'])
df['se'] = 1/np.sqrt(df['N'])

### Export to .ma
df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(snakemake.output[0], index = False, sep = '\t')


##########################################################################################
### Later 3 SNPs were removed from the file with the following command ###################
# awk -F'\t' '{empty=0; for(i=1; i<=NF; i++) if($i=="") empty=1; if(empty==0) print $0}' filename > modified_filename
# variant_id	chromosome	base_pair_location	effect_allele	other_allele	effect_allele_frequency	INFO	beta	standard_error	chisq	p_value
# rs80343720	6	109603416	A	G	0.999999	0.960556	-2819860000	NA	-1e+09	1
# rs78213867	3	85876178	T	C	0.999999	0.988221	2397200000	NA	-1e+09	1
# rs116540629	6	109603450	T	G	0.999999	0.960649	-2809550000	NA	-1e+09	1
##########################################################################################