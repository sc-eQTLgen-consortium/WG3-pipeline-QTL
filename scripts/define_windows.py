"""
Author:     D.Kaptijn & M.J.Bonder
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

CHROM = str(sys.argv[2])

window_size = 5000000 # size of sliding window
gene_window = 1000000 # window size around a gene

filepath_gene = str(sys.argv[1])
filepath_out = str(sys.argv[3])
input_filepath_feature = f"{filepath_gene}"

## INDIVIDUAL OUTFILES DEFINED AT END OF SCRIPT


"""
CODE
"""
print(f"\nRunning LD window script on chromosome {CHROM}\n")
print('reading and preparing gene features file..')

gene_info = pd.read_table(input_filepath_feature)
gene_info["chromosome"] = gene_info["chromosome"].astype(str)
gene_info = gene_info.loc[gene_info["chromosome"]==CHROM]

##Add midpoint
gene_info["mid"] = gene_info["start"] +((gene_info["end"]-gene_info["start"])/2)
##Sort on start.
gene_info = gene_info.sort_values(["start","mid"])

print("Creating LD windows..")

### defining the windows for LD
windows = {}
n = 0 # variable which defines key in windows dictionary
start = False # variable required for first loop.
for i in range(gene_info.shape[0]):
    if start == False:
      win_start = gene_info['start'].iloc[i] - gene_window
      block = win_start + window_size
      start = True
      end = int(gene_info['end'].iloc[i])
    
    if gene_info['mid'].iloc[i] < block and end < int(gene_info['end'].iloc[i]):
      end = int(gene_info['end'].iloc[i])
    elif gene_info['mid'].iloc[i] > block: # start of next window
      win_end = end + gene_window
      windows[str(n)] = {}
      windows[str(n)]['range'] = [win_start, win_end]
      windows[str(n)]['SNPs'] = []
      n+=1
      win_start = gene_info['start'].iloc[i] - gene_window
      block = win_start + window_size
      end = int(gene_info['end'].iloc[i])
    
    if i == (gene_info.shape[0]-1):
        ##Adding the last block
        win_end = end + gene_window
        windows[str(n)] = {}
        windows[str(n)]['range'] = [win_start, win_end]
        windows[str(n)]['SNPs'] = []

##### WRITE EACH WINDOW TO FILE #####

for i in windows:
    with open(f"{filepath_out}/chr_{CHROM}_window_{i}_{windows[i]['range'][0]}_{windows[i]['range'][1]}.pkl", 'wb') as fp:
        pickle.dump(windows[i], fp)
