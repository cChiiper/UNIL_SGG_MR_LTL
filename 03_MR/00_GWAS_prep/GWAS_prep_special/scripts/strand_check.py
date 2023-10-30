import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='strand-files')
parser.add_argument('-g','--gwas', help='GWAS summary data, COJO format', required=True)
parser.add_argument('-r','--reference', help='reference BIM file data', required=True)
parser.add_argument('-o','--output', help='Output file', required=True)

args = parser.parse_args()

print("Reading GWAS data...")
df = pd.read_csv(args.gwas, sep = '\t')
df['A1'] = df['A1'].str.upper()
df['A2'] = df['A2'].str.upper()

if 'freq' in list(df.columns):
    df = df.rename(columns = {'freq':'Freq'})

print("Number of SNPs:", df.shape[0])

print(df.head(2))

print("Removing A/T and C/G SNPs...")

df = df[~(((df.A1 == 'A') & (df.A2 == 'T')) | ((df.A1 == 'T') & (df.A2 == 'A')))]
df = df[~(((df.A1 == 'C') & (df.A2 == 'G')) | ((df.A1 == 'G') & (df.A2 == 'C')))]

print('Number of SNPs after filtering:', df.shape[0])

print("Reading reference data...")

ref_df = pd.read_csv(args.reference, sep = '\t', names = ['Chr', 'SNP', 'dum', 'bp_r', 'A1_r', 'A2_r'])

print(ref_df.head(2))

df = df.merge(ref_df, on = 'SNP')

print("Number of SNPs after merging with reference data", df.shape[0])

print("Detecting strand flipped SNPs...")
dic = {'C':'G', 'A':'T', 'G': 'C', 'T': 'A'}

def strand_flip(x):
    strand_check = True
    strand_flip = False
    if (x['A1'] == x['A1_r']) & (x['A2'] == x['A2_r']):
        s = '+'
        strand_check = False
    elif (x['A1'] == x['A2_r']) & (x['A2'] == x['A1_r']):
        s = 'i'
        strand_check = False
    if strand_check:
        for k in dic.keys():
            if (x['A1'] == k) & (x['A1_r'] == dic[k]):
                for k in dic.keys():
                    if (x['A2'] == k) & (x['A2_r'] == dic[k]):
                        s = '-'
                        strand_flip = True
    if strand_check & ~strand_flip:
        s = 'n'
        
    return s

df['strand'] = df.apply(lambda x: strand_flip(x), axis = 1)

df_count = df[['SNP', 'strand']].groupby('strand').count()

print("Strand detection results; +: same strand and order, -: inverse strand, i: same strand and inverse order, n: not found")

print(df_count)

print("Eliminate n, reverse i and correct -")

df = df[df.strand!='n']

df['b'] = df.apply(lambda x: x['b'] if (x['strand'] == '+') | (x['strand'] == '-') else -x['b'], axis = 1)
df['Freq'] = df.apply(lambda x: x['Freq'] if (x['strand'] == '+') | (x['strand'] == '-') else 1-x['Freq'], axis = 1)

df['A1'] = df['A1_r']
df['A2'] = df['A2_r']

print("Write Reference data checked GWAS summary data to file...")
df = df.drop_duplicates(subset = ['SNP'])

df[['SNP', 'A1', 'A2', 'Freq', 'b', 'se', 'p', 'N']].to_csv(args.output, sep = '\t', index = False)










