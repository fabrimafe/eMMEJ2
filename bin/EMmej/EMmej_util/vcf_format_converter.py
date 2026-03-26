"""
This script converts the formats:
REF     ALT
A       {.,-,''}

into:
REF     ALT
AC      A
"""

import argparse
import regex as re

import pandas as pd
import numpy as np

from pysam import FastaFile

# setting pandas display options
pd.options.display.max_colwidth = 3100
pd.set_option('display.max_columns', 500)

# Construct an argument parser
all_args = argparse.ArgumentParser()
# Add arguments to the parser
all_args.add_argument("-v", "--vcf", required=True,
   help="path to input vcf (string)")
all_args.add_argument("-o", "--outputfile", required=True,
      help="path to output file (string)")
all_args.add_argument("-r", "--ref", required=True,  
   help="path to reference genome in fasta format (string)")
all_args.add_argument("-af", "--allele_freq", required=True, 
      help="path to allele frequency data if available, if not, input=False (str)")


args = vars(all_args.parse_args())

# Set a FastaFile object
filename = args['ref']
refFA=FastaFile(args['ref'])

vcf = pd.read_csv(args['vcf'], sep='\t' ,comment='#', header=None)
vcf = vcf.iloc[:, 0:5]
columns = {0:'CHR', 1:'POS', 2:'ID', 3:'REF', 4:'ALT'}
vcf.rename(columns=columns, inplace=True)
vcf.drop(columns=['ID'], inplace=True)

# Only for the test data:
# vcf['CHR'] = 'NT_033779.5'

def get_ref_allele(chr:str, pos:int,
                    ref:str, alt:str):
    """
    A function that gets the base at position
    pos-1 and reconstruct the ALT and REF alleles
    in the right format
    Args:
        chr (str): chromosome number as in reference genome
        pos (int): position in chromosome
        ref (str): reference allele in vcf
        alt (str): alternative allele in vcf
    Returns:
        recostructed ref and alt
    """
    ref_base =  refFA.fetch(chr, (pos-1),(pos))
    if ref in ['', '.', '-']:
        ref,alt = ref_base, (ref_base+alt)
    else: ref,alt = ref,alt
    if alt in ['', '.', '-']:
        ref,alt = (ref_base+ref), ref_base
    else: ref,alt = ref,alt
    return (ref,alt)

corrected_ref_alt = vcf.apply(lambda row: get_ref_allele(chr=row['CHR'], pos=row['POS'],
                    ref=row['REF'], alt=row['ALT']), axis=1, result_type='expand')
vcf.loc[:, 'REF'] = corrected_ref_alt.iloc[:,0]
vcf.loc[:, 'ALT'] = corrected_ref_alt.iloc[:,1]


# Generating an ancestral state indication column
if args['allele_freq'] != 'False':
    columns = ['CHR', 'POS', 'N_ALLELES', 'N_CHR', 'FREQ_REF', 'FREQ_ALT']
    af_df = pd.read_csv(args['allele_freq'], sep='\t', names=columns, skiprows=1)
    af_df['ancestral_state'] = np.nan
    af_df.loc[(af_df['FREQ_REF'] > af_df['FREQ_ALT']), 'ancestral_state'] = 0
    af_df.loc[(af_df['FREQ_ALT'] > af_df['FREQ_REF']), 'ancestral_state'] = 1
    af_df.loc[(af_df['FREQ_ALT'] == af_df['FREQ_REF']), 'ancestral_state'] = 2
    vcf['ancestral_state'] = af_df['ancestral_state'].astype(int)


path_to_output = args['outputfile']
vcf.to_csv(f'{path_to_output}/EMmej_formated_data.vcf', sep='\t', header=False, index=False)

