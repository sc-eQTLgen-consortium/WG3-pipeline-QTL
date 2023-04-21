"""
Author:     D.Kaptijn
Date:       31/03/2023
PyVersion:  3.9.12
About:      This code loads windows to calculate LD within
"""

"""
MODULES TO IMPORT
"""

import sys
import numpy as np
import pandas as pd
import gzip
import pickle
from ldstore.bdose import bdose
from sklearn.datasets import fetch_openml
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score

"""
GLOBAL VARIABLES
"""
CHROM = str(sys.argv[4])

window_size = 5000000 # size of sliding window
gene_window = 1000000 # window size around a gene

filepath = str(sys.argv[1])
filepath_out = str(sys.argv[5])
input_filepath_window = f"{str(sys.argv[3])}"
input_filepath_bdose  = f"{str(sys.argv[2])}"

outfile = filepath_out[:-4] #removes file extenstion (.pkl)

outfile_low_dim = f'{outfile}_low_dim.pkl'
outfile_components = f'{outfile}_components.pkl'
outfile_means = f'{outfile}_means.pkl'


"""
CODE
"""

### Loading calculated windows
print("Loading windows..")
with open(input_filepath_window, 'rb') as fp:
    window = pickle.load(fp)


### Sorting bdose file
print("Reading bdose file..")
myBdose = bdose(input_filepath_bdose)
CHROM=str(myBdose.getMeta()['chromosome'][1])
SNP_position_by_chromosome = {}
SNP_position_by_chromosome[CHROM] = []
for i in range(len(myBdose.getMeta()['chromosome'])):
    SNP_position_by_chromosome[CHROM].append(myBdose.getMeta()['position'].iloc[i])

# SNP_position_by_chromosome is a dictionary in which the keys relate to the chromosome and the values relate to the position of the SNP


### loop to collect SNPs per window
print("Collecting SNPs per window..")
myBdoseChrom = myBdose.getMeta()[myBdose.getMeta()['chromosome']==CHROM]
total_SNPs = len(SNP_position_by_chromosome[CHROM])
n=0
for SNP in SNP_position_by_chromosome[CHROM]:
    progress = f"Progress:" + f" {round((n/total_SNPs)*100)}%"
    n+=1
    print(progress, end='\r')
    if int(SNP) >= window['range'][0] and int(SNP) <= window['range'][1]:
        index_SNP = myBdoseChrom[myBdoseChrom['position'] == SNP].index[0]
        window['SNPs'].append(int(index_SNP))


### calculating LD and saving output file
print("Calculating LD and saving output file..")
f = open(outfile_low_dim, 'w') # this will clear the contents of the file
f.close()
f = open(outfile_components, 'w') # this will clear the contents of the file
f.close()
f = open(outfile_means, 'w') # this will clear the contents of the file
f.close()

if len(window["SNPs"])>0 :
    n_samples = int(myBdose.getNumOfSamples()) # required to select number of PCs
    no_components = min([len(window['SNPs']), n_samples-1])
    pca = PCA(no_components)
    LD_matrix = myBdose.computeCorr(window['SNPs'])
    LD_matrix = LD_matrix.to_numpy()
    lower_dimensional_data = pca.fit_transform(LD_matrix)
    components = pca.components_
    means = pca.mean_
    print(f"Shape of lower_dim_data: {lower_dimensional_data.shape}")
    with open(outfile_low_dim, 'ab') as outp:
        pickle.dump(lower_dimensional_data, outp)
    outp.close()
    print(f"Shape of components: {components.shape}")
    with open(outfile_components, 'ab') as outp:
        pickle.dump(components, outp)
    outp.close()
    print(f"Shape of means: {means.shape}\n")
    with open(outfile_means, 'ab') as outp:
         pickle.dump(means, outp)
    outp.close()

print("\nScript finished.")

