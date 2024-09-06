#!/usr/bin/env python
# Author: M. Vochteloo

import argparse
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--wg1_metadata", required=True, type=str, help="")
parser.add_argument("--individual_aggregate", required=False, default="Assignment", type=str, help="")
parser.add_argument("--ancestry", required=True, type=str, help="")
parser.add_argument("--output_dir", required=True, type=str, help="")
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


print("Extracting {} values from WG1 metadata:".format(args.individual_aggregate))

individual_col = None
individuals = set()
with gzopen(args.wg1_metadata, mode='r') as f:
    for line in f:
        values = line.strip("\n").split("\t")
        if individual_col is None:
            individual_col = values.index(args.individual_aggregate)
            continue

        if values[individual_col] not in individuals:
            individuals.add(values[individual_col])
f.close()

print("\tFound {:,} individuals".format(len(individuals)))
if len(individuals) == 0:
    exit()

with gzopen(os.path.join(args.output_dir, args.ancestry + "_individuals.txt"), mode='w') as f:
    for individual in individuals:
        f.write(individual + "\n")
f.close()

print("Done")
