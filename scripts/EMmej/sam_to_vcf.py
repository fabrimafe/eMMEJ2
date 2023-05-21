"""
This script converts a .sam file into .vcf-like format.
vcf-like format:
CHR     POS     ALT     REF
"""
import re
import argparse

import pandas as pd
import numpy as np
import csv

# Construct an argument parser
all_args = argparse.ArgumentParser()
# Add arguments to the parser
all_args.add_argument("-sam", "--inputsamfile", required=True,
   help="path to input sam file (string)")
all_args.add_argument("-O", "--outputfile", required=True,
      help="path to output file (string)")
all_args.add_argument("-r", "--ref", required=True,
   help="path to reference genome in fasta format (string)")
args = vars(all_args.parse_args())


with open(args['ref']) as f:
    lines = f.readlines()
for lin in lines:
    if lin.startswith('>'): continue
    else: ref = lin.split('\n')[0]
print(ref)
# ref = 'GGAAAAAATTCGTACTTTGGAGTACGAAATGCGTCGTTTAGAGCAGCAGCCGAATTCGGTACATTACCCTGTTATCCCTAGCGGCCGCATAGGCCACTAGTGGATCTGGATCCTCTAGAGTCGACCTCGAACGTTAACGTTAACGTAACGTTAACTCG'
# ref = args['ref']

path_to_df = args['inputsamfile']
# reading sam
df = [i for i in csv.reader(open(path_to_df, "r"), delimiter = '\t')]

# convert sam to pd.DataFrame
data = pd.DataFrame(columns=[
    'QNAME', 'FLAG', 'RNAME', 'align_POS',
    'MAPQ', 'CIGAR', 'SEQ', 'QUAL'])

# print(data)

for idx,i in enumerate(df[2:]): # skipping the 1st two rows of the sam
    data.loc[idx, 'QNAME']= i[0]
    data.loc[idx, 'FLAG']= i[1] 
    data.loc[idx, 'RNAME']= i[2]
    data.loc[idx, 'align_POS']= i[3]
    data.loc[idx, 'MAPQ']= i[4]
    data.loc[idx, 'CIGAR']= i[5]
    data.loc[idx, 'SEQ']= i[9]
    data.loc[idx, 'QUAL']= i[10]
# print(data)

def parse_cigar(ref: str, read: str, cigar: str):
    """
    A function that reconstruct the REF and ALT
    sequenses based on the CIGAR srting
    Args:
        res (str): reference sequence
        read (str): read sequence
        cigar (str): CIGAR string
    returns:
        POS, REF, ALT as in vcf format
    """
    matches = re.findall(r'(\d+)([A-Z]{1})', cigar)
    cigar_split = [(m[1],int(m[0])) for m in matches]
    index = 0
    # conditioning cigar parsing to operate only on simple cigars
    cigar_chr_count = [cigar.count(c) for c in ['I', 'D']]
    # either both I and D, or more then 1 I or D
    complex_cigar = (('I' in cigar) & ('D' in cigar) | (max(cigar_chr_count)>1))
    # only M
    total_match = all([c == 0 for c in cigar_chr_count])
    if complex_cigar == total_match == False:
        for i in cigar_split:
            if i[0] == 'M': 
                index += i[1]
                continue
            elif i[0] == 'D':
                REF = ref[index:(index+i[1]+1)]
                ALT = read[index:(index+1)]
                POS = index
                index += i[1]
            elif i[0] == 'I':
                REF = ref[(index-1):(index)]
                ALT = read[(index-1):(index+i[1])]
                POS = index
                index += i[1]

        return [POS, REF, ALT]
    else: return [np.nan, np.nan , np.nan]

print(data.info())
print(data.loc[data['CIGAR'].isin(data['CIGAR'].unique()), :].info())

data[['POS', 'REF', 'ALT']] = data.apply(
    lambda row: parse_cigar(ref=ref, read=row['SEQ'], 
            cigar=row['CIGAR']), axis=1, result_type='expand')

vcf = data.loc[
        (data['REF'].isna() == False), 
        ['RNAME', 'POS', 'REF', 'ALT']].reset_index(drop=True)

vcf['POS'] = vcf['POS'].astype('int')
# print(vcf)
path_to_output = args['outputfile']
# exporting vcf
vcf.to_csv(f'{path_to_output}.vcf', sep='\t', header=False, index=False)


# Running command example: python3 sam_to_vcf.py -sam /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/raw_data_to_SAM/1806_20220811_R0_WT/bwa_output/bwa_output.sam -r GGAAAAAATTCGTACTTTGGAGTACGAAATGCGTCGTTTAGAGCAGCAGCCGAATTCGGTACATTACCCTGTTATCCCTAGCGGCCGCATAGGCCACTAGTGGATCTGGATCCTCTAGAGTCGACCTCGAACGTTAACGTTAACGTAACGTTAACTCG -O output_path