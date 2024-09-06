#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import gzip
import glob
import json
import numpy as np
import pandas as pd
# from natsort import natsort_keygen
import time
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True, type=str, help="")
parser.add_argument("--output", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


print("Reading run_qtl_mapping logfiles")
data = []
for fpath in glob.glob(args.input + "*"):
    found = False
    with gzopen(fpath, mode='r') as f:
        for line in f:
            if line.startswith("{") and prev_line == "saving log!\n":
                runtime = json.loads(line.rstrip("\n").replace("'", "\""))
                for key in runtime:
                    if not runtime[key]:
                        runtime[key] = [np.nan, np.nan]

                df = pd.DataFrame(runtime).T
                df.reset_index(drop=False, inplace=True)
                df.columns = ["Gene", "Time", "Mean"]
                found = True
                break

            prev_line = line
    f.close()
    print("\t{}\t{}".format(os.path.basename(fpath), found))

    if found:
        # run_qtl_mapping.{ancestry}.{cell_level}.{cell_type}.chr_{qtl_chunk}.log
        rule, ancestry, cell_level, cell_type, chunk, _ = os.path.basename(fpath).split(".")
        chromosome, start, end = chunk.lstrip("chr_").split("_")

        df["Ancestry"] = ancestry
        df["Cell level"] = cell_level
        df["Cell type"] = cell_type
        df["Chromosome"] = int(chromosome)
        df["Start"] = int(start)
        df["End"] = int(end)
        df["Genomic range"] = "{}:{}-{}".format(chromosome, start, end)
        data.append(df)

df = pd.concat(data, axis=0)
if len(data) == 0 or df.shape[0] == 0:
    print("\nNo data to combine.")
    exit()

print("\nN-genes:\t{:,}".format(df.shape[0]))
print("Min time:\t{:.2f} seconds".format(df["Time"].min()))
print("Max time:\t{:.2f} seconds".format(df["Time"].max()))
print("Average time:\t{:.2f} seconds".format(df["Time"].mean()))
print("Total time:\t{}".format(time.strftime('%H hours %M minutes %S seconds', time.gmtime(df["Time"].sum()))))

print("\nTime per gene:")
# df.sort_values(by=["Chromosome", "Start", "Gene"], key=natsort_keygen(), ascending=True, inplace=True)
df.sort_values(by=["Chromosome", "Start", "Gene"], ascending=True, inplace=True)
print(df)

print("\nTime per chunk:")
print(df.groupby([col for col in df.columns if col not in ["Gene", "Time", "Mean"]]).sum(["Time", "Mean"]).reset_index(drop=False))

print("\nSaving output file")
df.to_csv(args.output, sep="\t", header=0, index=False, compression="gzip")

print("Done")
