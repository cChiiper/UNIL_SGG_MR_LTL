import numpy as np
import pandas as pd

# Data from  DOI: 10.1016/j.xhgg.2023.100201 https://www.sciencedirect.com/science/article/pii/S2666247723000337?via%3Dihub
# Thanks for the early access to the data to Patrick Allaire et al. (2023)
# File headers are:
# "SNP"	"BETA"	"P"	"Direction"	"SE"	"TEST"	"OTHER"	"FREQ_TEST"	"RS"	"CHR"	"POS"	"Gene"	"Func"	"ExonicFunc"
# "10:87086125:A:G"	-0.01	0.1169	"--"	0.0064	"A"	"G"	0.7019	"rs12779881"	10	87086125	"LINC01519;LOC101929646"	"intergenic"	""

df = pd.read_csv(snakemake.input[0], compression='gzip', sep = '\t')

pd.set_option('display.max_columns', None)
print(df.head(2))

# Select specific columns to avoid confusion with SNP
columns_to_keep = ['RS', 'TEST', 'OTHER', 'FREQ_TEST', 'P', 'BETA', 'SE']
df = df[columns_to_keep]

df = df.rename(columns = {'RS':'SNP', 'TEST': 'A1', 'OTHER': 'A2', 'FREQ_TEST':'freq', 'P':'p'})

df['N'] = 62271

df['b'] = df['BETA']/df['SE']/np.sqrt(df['N'])
df['se'] = 1/np.sqrt(df['N'])

print("Formated df head:")
print(df.head(2))
df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(snakemake.output[0], index = False, sep = '\t')