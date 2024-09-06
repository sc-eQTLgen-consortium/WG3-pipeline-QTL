#!/usr/bin/env python
# Author: M. Vochteloo (Adapted from 'https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/blob/master/scripts/process_output.R')

import argparse
import os
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True, type=str, help="")
parser.add_argument("--maf", required=True, type=float, default=0.01, help="")
parser.add_argument("--r2", required=True, type=float, default=0.6, help="")
parser.add_argument("--call", required=True, type=float, default=0.99, help="")
parser.add_argument("--out_dir", required=True, type=str, help="")
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


print("Reading and filtering variants ...")
fhin = gzopen(args.input, mode='r')
fhout = gzopen(os.path.join(args.out_dir, os.path.basename(args.input).replace("_stats", "_stats_filtered")), mode='w')

chr_variants = {}
coi = {"ID": None, "CHR": None, "MAF": None, "MACH_R2": None, "CALL": None}
maf_filter = 0
r2_filter = 0
call_filter = 0
n_all_pass = 0
i = 0
for i, line in enumerate(fhin):
    if (i == 0) or ((i % 1000000) == 0):
        print("\tParsed {:,} lines\tMAF removed: {:,}\tMACH_R2 removed: {:,}\tCALL removed: {:,}\tPassed: {:,}".format(i, maf_filter, r2_filter, call_filter, n_all_pass), end='\r',  flush=True)
    values = line.strip("\n").split("\t")
    if i == 0:
        coi["ID"] = values.index("ID")
        coi["CHR"] = values.index("CHR")
        coi["MAF"] = values.index("MAF")
        coi["MACH_R2"] = values.index("MACH_R2")
        coi["CALL"] = values.index("CALL")
        fhout.write(line)
        continue

    all_pass = True
    if float(values[coi["MAF"]]) < args.maf:
        maf_filter += 1
        all_pass = False

    if float(values[coi["MACH_R2"]]) < args.r2:
        r2_filter += 1
        all_pass = False

    if float(values[coi["CALL"]]) < args.call:
        call_filter += 1
        all_pass = False

    if all_pass:
        n_all_pass += 1
        fhout.write(line)

        if values[coi["CHR"]] in chr_variants:
            if values[coi["ID"]] in chr_variants[values[coi["CHR"]]]:
                print("Error, duplicate ID '{}'".format(values[coi["ID"]]))
                exit()
            chr_variants[values[coi["CHR"]]].add(values[coi["ID"]])
        else:
            chr_variants[values[coi["CHR"]]] = {values[coi["ID"]]}
print("\tParsed {:,} lines\tMAF removed: {:,}\tMACH_R2 removed: {:,}\tCALL removed: {:,}\tPassed: {:,}".format(i, maf_filter, r2_filter, call_filter, n_all_pass), flush=True)
fhin.close()
fhout.close()

if n_all_pass == 0:
    print("Error, no variants passed filter.")
    exit()

print("\nWriting variants per chromosome")
for chr, variants in chr_variants.items():
    if len(variants) == 0:
        print("\tChromosome: {}\tNo variants.".format(chr))
        continue
    print("\tChromosome: {}\t{:,} variants.".format(chr, len(variants)))
    with gzopen(os.path.join(args.out_dir, os.path.basename(args.input).replace("_stats.vars.gz", "_inclusion_{}.vars".format(chr))), mode='w') as f:
        for variant in variants:
            f.write(variant + "\n")
    f.close()

print("Done")
