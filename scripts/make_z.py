"""
Author:     D.Kaptijn
Date:       14/10/2022
PyVersion:  3.8.8
About:      This code is designed to generate a z file as described in LDstore2
"""

""" MODULE IMPORTS """

import sys
import gzip
from bgen_reader import read_bgen


""" GLOBAL VARIABLES """

"""input"""
in_filepath = str(sys.argv[1])
out_filepath = str(sys.argv[2])
infilename = in_filepath.split("/")[len(in_filepath.split("/"))-1]

filepath_vcf = in_filepath
filepath_bgen = in_filepath

"""output"""
filepath_output_z      = f'{out_filepath}outputs/{infilename}.z'
filepath_output_sample = f'{out_filepath}outputs/{infilename}.sample'
filepath_output_master = f'{out_filepath}outputs/{infilename}_master.txt'


""" FILE LOCATIONS FOR MASTER FILE """

z_file            = f"{out_filepath}outputs/{infilename}.z"
bgen_file         = f"{in_filepath}"
bgi_file          = f"{in_filepath}.bgi"
sample_file       = f"{out_filepath}outputs/{infilename}.sample"
bdose_file_output = f"{out_filepath}outputs/test_dataset.bdose"


""" READ INPUT VCF """

try:
    with gzip.open(filepath_vcf, 'r') as f:
        vcf=[(line.decode()).strip() for line in f.readlines()]
    verify = "vcf"
except:
    bgen = read_bgen(filepath_bgen, verbose=False)
    verify = 'bgen'
    variant = bgen["variants"].compute()

print(f"file type: {verify}")

data_z = {}

if verify == "vcf":
    n=0
    for line in vcf:
        if '#' not in line:
            start_line = n
            break
        n+=1

    vcf=vcf[start_line-1:]
    vcf_dict = {}

    for i in vcf[0].split('\t'):
        vcf_dict[i] = []

    keys = [i for i in vcf_dict]
    vcf = vcf[1:]

    for i in range(len(vcf)):
        temp = vcf[i].split('\t')
        for j in range(len(keys)):
            vcf_dict[keys[j]].append(temp[j])


""" CREATE Z FILE """

if verify == "vcf":
    data_z = {}
    data_z['rsid'] = [vcf_dict['ID'][i] for i in range(len(vcf_dict['ID']))]
    data_z['chromosome'] = [vcf_dict['#CHROM'][i] for i in range(len(vcf_dict['ID']))]
    data_z['position'] = [vcf_dict['POS'][i] for i in range(len(vcf_dict['ID']))]
    data_z['allele1'] = [vcf_dict['REF'][i] for i in range(len(vcf_dict['ID']))]
    data_z['allele2'] = [vcf_dict['ALT'][i]for i in range(len(vcf_dict['ID']))]

if verify == "bgen":
    data_z['rsid'] = [variant['rsid'].iloc[i] for i in range(len(variant['rsid']))]
    data_z['chromosome'] = [variant['chrom'].iloc[i] for i in range(len(variant['rsid']))]
    data_z['position'] = [variant['pos'].iloc[i] for i in range(len(variant['rsid']))]
    data_z['allele1'] = [variant['allele_ids'].iloc[i].split(',')[0] for i in range(len(variant['rsid']))]
    data_z['allele2'] = [variant['allele_ids'].iloc[i].split(',')[1] for i in range(len(variant['rsid']))]

print(f"size of data {len(data_z['rsid'])}")
print("writing output z file..")

with open(filepath_output_z, 'w') as f:
    f.write(f'rsid chromosome position allele1 allele2 \n')
    for i in range(len(data_z['rsid'])):
        f.write(f'{data_z["rsid"][i]} {data_z["chromosome"][i]} {data_z["position"][i]} {data_z["allele1"][i]} {data_z["allele2"][i]} \n')
    f.close()



""" CODE TO CREATE SAMPLE FILE """

data_sample = {}

if verify == "vcf":
    list_of_samples = [i for i in vcf_dict if i not in ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']]
    data_sample['ID_1'] = [0] 
    data_sample['ID_2'] = [0]
    for i in list_of_samples:
        data_sample['ID_1'].append(i)
        data_sample['ID_2'].append(i)
    data_sample['missing'] = [0 for i in range(len(list_of_samples)+1)]

if verify == 'bgen':
    data_sample['ID_1'] = [0]
    data_sample['ID_2'] = [0]
    for i in range(len(bgen['samples'])):
        data_sample['ID_1'].append(bgen['samples'].iloc[i])
        data_sample['ID_2'].append(bgen['samples'].iloc[i])
    data_sample['missing'] = [0 for i in range(len(bgen['samples'])+1)]


print("writing output sample file..")

with open(filepath_output_sample, 'w') as f:
    f.write(f'ID_1 ID_2 missing \n')
    for i in range(len(data_sample['ID_1'])):
        f.write(f"{data_sample['ID_1'][i]} {data_sample['ID_2'][i]} {data_sample['missing'][i]} \n")
    f.close()


""" CODE TO WRITE MASTER FILE """

data_master = {}

if verify == "vcf":
    n_samples = len(list_of_samples)

if verify == 'bgen':
    n_samples = len(bgen['samples'])

data_master['z'] = z_file
data_master['bgen'] = bgen_file
data_master['bgi'] = bgi_file
data_master['sample'] = sample_file
data_master['bdose'] = bdose_file_output
data_master['n_samples'] = n_samples

print("writing output master file..")

with open(filepath_output_master, 'w') as f:
    for key in data_master:
        f.write(f'{key};')
    f.write('\n')
    for key in data_master:
        f.write(f"{data_master[key]};")
    f.close()


