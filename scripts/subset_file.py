#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--data", required=True, type=str, help="")
parser.add_argument("--n_columns", required=True, type=int, help="")
parser.add_argument("--no_index", dest="has_index", action="store_false", default=True, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

n_columns = args.n_columns
if args.has_index:
    n_columns += 1

sep = "\t"


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

# Doing this check first to prevent the output file from being half created and snakemake
# thinking it is done without issues.
fhin = gzopen(args.data, mode="r")
total_columns = len(fhin.readline().strip("\n").split(sep))
fhin.close()

print("Input file has {:,} columns.".format(total_columns))
if n_columns > total_columns:
    print("Error, selecting the first {} columns which is higher than the total number of columns.".format(n_columns))
    exit()

fhin = gzopen(args.data, mode="r")
fhout = gzopen(args.out, mode="w")
for line in fhin:
    values = line.strip("\n").split(sep)
    fhout.write(sep.join(values[:n_columns]) + "\n")
fhin.close()
fhout.close()

print("Done")
