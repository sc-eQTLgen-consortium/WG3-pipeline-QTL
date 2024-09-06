"""
Author:     D.Kaptijn and a bit of M.Vochteloo
Date:       31/03/2023
PyVersion:  3.9.12
About:      This code loads windows to calculate LD within
"""

"""
MODULES TO IMPORT
"""

import argparse
import re
import gzip
import pickle
from ldstore.bdose import bdose
from sklearn.decomposition import PCA

"""
GLOBAL VARIABLES
"""

"""input"""
parser = argparse.ArgumentParser(description="")
parser.add_argument("--bdose", required=True, type=str, help="")
parser.add_argument("--outfile", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

outfile_low_dim = f'{args.out}{args.outfile}_low_dim.pkl.gz'
outfile_components = f'{args.out}{args.outfile}_components.pkl.gz'
outfile_means = f'{args.out}{args.outfile}_means.pkl.gz'
outfile_window = f'{args.out}{args.outfile}.pkl.gz'

match = re.match("chr_([0-9]{1,2}|X|Y|MT)_window_([0-9]+)_([0-9]+)_([0-9]+)", args.outfile)
chr = str(match.group(1))
num = int(match.group(2))
start = int(match.group(3))
end = int(match.group(4))

### Sorting bdose file
print("Reading bdose file..")
myBdose = bdose(args.bdose)

print("Collecting SNPs per window..")
snp_indices = []
snp_positions = []
for index, row in myBdose.getMeta().iterrows():
    bdose_chr = str(row["chromosome"])
    bdose_position = int(row["position"])
    if bdose_chr != chr or bdose_position < start or bdose_position > end:
        continue

    snp_indices.append(index)
    snp_positions.append(bdose_position)

### calculating LD and saving output file
print("Calculating LD and saving output file..")
f = open(outfile_window, 'w')  # this will clear the contents of the file
f.close()
f = open(outfile_low_dim, 'w')  # this will clear the contents of the file
f.close()
f = open(outfile_components, 'w')  # this will clear the contents of the file
f.close()
f = open(outfile_means, 'w')  # this will clear the contents of the file
f.close()

print(f"  Number of variants: {len(snp_positions)}")
with gzip.open(outfile_window, 'ab') as outp:
    pickle.dump({"chr": chr, "num": num, "start": start, "end": end, "SNPs": snp_positions}, outp)
outp.close()


if len(snp_indices) > 0:
    n_samples = int(myBdose.getNumOfSamples())  # required to select number of PCs
    no_components = min([len(snp_indices), n_samples-1])

    print(f"  Number of samples: {n_samples}")
    print(f"  Number of components: {no_components}")

    pca = PCA(no_components)

    LD_matrix = myBdose.computeCorr(snp_indices)
    LD_matrix = LD_matrix.to_numpy()
    lower_dimensional_data = pca.fit_transform(LD_matrix)

    print(f"Shape of lower_dim_data: {lower_dimensional_data.shape}")
    with gzip.open(outfile_low_dim, 'ab') as outp:
        pickle.dump(lower_dimensional_data, outp)
    outp.close()

    print(f"Shape of components: {pca.components_.shape}")
    with gzip.open(outfile_components, 'ab') as outp:
        pickle.dump(pca.components_, outp)
    outp.close()

    print(f"Shape of means: {pca.mean_.shape}\n")
    with gzip.open(outfile_means, 'ab') as outp:
        pickle.dump(pca.mean_, outp)
    outp.close()

print("\nScript finished.")
