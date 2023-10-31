#############################################
### Import libraries ########################
import pandas as pd

#############################################
### Load data ###############################
exp_snps = pd.read_csv(snakemake.input['snps'], sep = '\t', names = ['SNP'])
print("Exposure SNPs loaded")
print(exp_snps.head())

### Read by chuncksizes
chunksize = 10 ** 6
dtype = {
    "SNP": str,
    "A1": str,
    "A2": str,
    "freq": float,
    "b": float,
    "se": float,
    "p": float,
    "N": float,
}

reader = pd.read_csv(snakemake.input['gwas'], chunksize=chunksize, sep = '\t', dtype=dtype, float_precision='high')

# Initializing an empty dataframe with the specified columns
filtered_data = pd.DataFrame(columns=dtype.keys())

# Appending filtered data for each chunk
for chunk in reader:
    filtered_chunk = chunk[chunk["SNP"].isin(exp_snps['SNP'])]
    filtered_data = pd.concat([filtered_data, filtered_chunk], ignore_index=True)

#############################################
### Save file ###############################
filtered_data[['SNP', 'A1', 'A2', 'Freq', 'b', 'se', 'p', 'N']].to_csv(snakemake.output[0], index = False, header = True,  sep='\t')
