"""
Author:     D.Kaptijn & M.J.Bonder and a bit of M.Vochteloo
Date:       31/03/2023
PyVersion:  3.9.12
About:      This code defines windows to calculate LD in
"""

"""
MODULES TO IMPORT
"""

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="")
parser.add_argument("--feature_file", required=True, type=str, help="")
parser.add_argument("--window_size", required=True, type=int, default=5000000, help="size of sliding window")
parser.add_argument("--gene_window", required=True, type=int, default=1000000, help="window size around a gene")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

print(f"\nRunning LD window script\n")
print('reading and preparing gene features file..')

gene_info = pd.read_table(args.feature_file)
gene_info["chromosome"] = gene_info["chromosome"].astype(str)

##Add midpoint
gene_info["mid"] = gene_info["start"] + ((gene_info["end"] - gene_info["start"]) / 2)

lines = []
for chr in range(1, 23):
    chr = str(chr)
    print(f"  Creating LD windows for chromosome {chr}..")

    chr_gene_info = gene_info.loc[gene_info["chromosome"] == chr].copy()
    chr_gene_info = chr_gene_info.sort_values(["start", "mid"])

    ### defining the windows for LD
    n = 0  # variable which defines key in windows dictionary
    start = False  # variable required for first loop.
    for i in range(chr_gene_info.shape[0]):
        if start == False:
            win_start = chr_gene_info['start'].iloc[i] - args.gene_window
            block = win_start + args.window_size
            start = True
            end = int(chr_gene_info['end'].iloc[i])

        if chr_gene_info['mid'].iloc[i] < block and end < int(chr_gene_info['end'].iloc[i]):
            end = int(chr_gene_info['end'].iloc[i])
        elif chr_gene_info['mid'].iloc[i] > block:  # start of next window
            win_end = end + args.gene_window
            lines.append(f"{chr}:{n}:{max(0, win_start)}-{win_end}")
            n += 1
            win_start = chr_gene_info['start'].iloc[i] - args.gene_window
            block = win_start + args.window_size
            end = int(chr_gene_info['end'].iloc[i])

        if i == (chr_gene_info.shape[0] - 1):
            ##Adding the last block
            win_end = end + args.gene_window
            lines.append(f"{chr}:{n}:{max(0, win_start)}-{win_end}")

##### WRITE WINDOWS TO FILE #####
print(f"Created {len(lines)} chunks")

with open(f"{args.out}LDChunkingFile.txt", 'w') as f:
    for line in lines:
        f.write(line + "\n")
f.close()
