"""
Py Version: 3.9.16
Date: 13/11/2023
Author: D. Kaptijn

Input QTL results and output a file per cell type
"""

import pandas as pd
import numpy as np
import sys

# variables defined in snakemake
infile = sys.argv[1]
outfile = sys.argv[2]
chrom = sys.argv[3]

print(f"Input file: {infile}")
print(f"Output path: {outfile}")
print(f"Chromosome: {chrom}")

# find number of values in header, subsequent lines can have more values than expected
header = open(infile, 'r').readline().split('\t')
n_cols = len(header)
print(f"There are {n_cols} values in header of the file")

# read in file and separate into table per cell type
df = pd.read_csv(infile,sep='\t',usecols=[i for i in range(n_cols)])

print(df)

dict = {}

for ct in np.unique(df["feature_id"]):
    better_ct = ct.replace(" ", "_")
    dict[better_ct] = df[df["feature_id"]==ct]


print(f"Number of features in results: {len([i for i in dict])}")
# write result
for key in dict:
    dict[key].to_csv(f"{outfile}/chr{chrom}_{key}.txt.gz", sep='\t', index=None, compression="gzip")

# also save an empty file which snakemake expects as output
open(f"{outfile}/ct-split-done.txt","w")

print("done.")

