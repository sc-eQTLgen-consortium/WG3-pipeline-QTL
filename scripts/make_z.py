"""
Author:     D.Kaptijn
Date:       14/10/2022
PyVersion:  3.8.8
About:      This code is designed to generate a z file as described in LDstore2
"""

""" MODULE IMPORTS """

import argparse
from bgen_reader import read_bgen
import gzip


""" GLOBAL VARIABLES """

"""input"""
parser = argparse.ArgumentParser(description="")
parser.add_argument("--in_filepath", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

"""output"""
filepath_output_z      = f'{args.out}.z'
filepath_output_sample = f'{args.out}.sample'
filepath_output_master = f'{args.out}_master.txt'

""" FILE OPEN FUNCTION """


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


""" FILE LOCATIONS FOR MASTER FILE """

z_file            = f"{filepath_output_z}"
bgen_file         = f"{args.in_filepath}"
bgi_file          = f"{args.in_filepath}.bgi"
sample_file       = f"{filepath_output_sample}"
bdose_file_output = f"{args.out}.bdose"


""" READ INPUT VCF """

bgen = read_bgen(args.in_filepath, verbose=False)
variant = bgen["variants"].compute()

data_z = {}


""" CREATE Z FILE """

data_z['rsid'] = [variant['rsid'].iloc[i] for i in range(len(variant['rsid']))]
data_z['chromosome'] = [variant['chrom'].iloc[i] for i in range(len(variant['rsid']))]
data_z['position'] = [variant['pos'].iloc[i] for i in range(len(variant['rsid']))]
data_z['allele1'] = [variant['allele_ids'].iloc[i].split(',')[0] for i in range(len(variant['rsid']))]
data_z['allele2'] = [variant['allele_ids'].iloc[i].split(',')[1] for i in range(len(variant['rsid']))]

print(f"size of data {len(data_z['rsid'])}")
print("writing output z file..")

with gzopen(filepath_output_z, mode='w') as f:
    f.write(f'rsid chromosome position allele1 allele2 \n')
    for i in range(len(data_z['rsid'])):
        f.write(f'{data_z["rsid"][i]} {data_z["chromosome"][i]} {data_z["position"][i]} {data_z["allele1"][i]} {data_z["allele2"][i]} \n')
f.close()



""" CODE TO CREATE SAMPLE FILE """

data_sample = {}

data_sample['ID_1'] = [0]
data_sample['ID_2'] = [0]
for i in range(len(bgen['samples'])):
    data_sample['ID_1'].append(bgen['samples'].iloc[i])
    data_sample['ID_2'].append(bgen['samples'].iloc[i])
data_sample['missing'] = [0 for i in range(len(bgen['samples'])+1)]


print("writing output sample file..")

with gzopen(filepath_output_sample, mode='w') as f:
    f.write(f'ID_1 ID_2 missing \n')
    for i in range(len(data_sample['ID_1'])):
        f.write(f"{data_sample['ID_1'][i]} {data_sample['ID_2'][i]} {data_sample['missing'][i]} \n")
f.close()


""" CODE TO WRITE MASTER FILE """

data_master = {}

n_samples = len(bgen['samples'])

data_master['z'] = z_file
data_master['bgen'] = bgen_file
data_master['bgi'] = bgi_file
data_master['sample'] = sample_file
data_master['bdose'] = bdose_file_output
data_master['n_samples'] = n_samples

print("writing output master file..")

with gzopen(filepath_output_master, mode='w') as f:
    for key in data_master:
        f.write(f'{key};')
    f.write('\n')
    for key in data_master:
        f.write(f"{data_master[key]};")
f.close()
