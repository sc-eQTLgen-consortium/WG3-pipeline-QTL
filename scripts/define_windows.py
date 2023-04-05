"""
Author:     D.Kaptijn
Date:       31/03/2023
PyVersion:  3.9.12
About:      This code defines windows to calculate LD in
"""

"""
MODULES TO IMPORT
"""

import sys
import numpy as np
import pandas as pd
import gzip
import pickle

CHROM = str(sys.argv[3])

window_size = 5000000 # size of sliding window
gene_window = 1000000 # window size around a gene

filepath = str(sys.argv[1])
filepath_gene = str(sys.argv[2])
filepath_out = str(sys.argv[4])
input_filepath_feature = f"{filepath_gene}"

## OUTFILE DEFINED AT END OF SCRIPT


"""
CODE
"""
print(f"\nRunning LD window script on chromosome {CHROM}\n")
print('reading gene features file..')

with gzip.open(input_filepath_feature, 'r') as f:
    genes=[(line.decode()).strip() for line in f.readlines()]

genes_dict = {}
for i in genes[0].split('\t'):
    genes_dict[i] = []

col_names = [i for i in genes_dict]

genes = genes[1:] # header is no longer required
for i in range(len(genes)):
    temp = genes[i].split('\t')
    for j in range(len(col_names)):
        if col_names[j] == 'chromosome' and str(temp[j]) == CHROM:
            for k in range(len(col_names)):
                genes_dict[col_names[k]].append(temp[k])
# genes_dict is a dictionary with a key for each heading in the genes features (e.g. feature_ID) file and a list of values that correspond to each key

print("Creating LD windows..")

### defining the windows for LD
windows = {}
n = 0 # variable which defines key in windows dictionary
start = False # variable required for first loop
for i in range(len(genes_dict['start'])):
    if start == False:
      win_start = int(genes_dict['start'][i]) - gene_window
      block = win_start + window_size
      start = True
      end = int(genes_dict['end'][i])
    mid_point_gene = int(genes_dict['start'][i]) + (int(genes_dict['end'][i]) - int(genes_dict['start'][i]))/2
    if mid_point_gene < block and end < int(genes_dict['end'][i]):
      end = int(genes_dict['end'][i])
    elif mid_point_gene > block: # start of next window
      win_end = end + gene_window
      windows[str(n)] = {}
      windows[str(n)]['range'] = [win_start, win_end]
      windows[str(n)]['SNPs'] = []
      n+=1
      win_start = int(genes_dict['start'][i]) - gene_window
      block = win_start + window_size
      end = int(genes_dict['end'][i])


##### WRITE EACH WINDOW TO FILE #####

for i in windows:
    with open(f"{filepath_out}/chr_{CHROM}_window_{i}_{windows[i]['range'][0]}_{windows[i]['range'][1]}.pkl", 'wb') as fp:
        pickle.dump(windows[i], fp)

