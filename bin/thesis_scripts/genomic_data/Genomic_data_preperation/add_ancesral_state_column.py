"""
A script that adds an "ancestral_state" indication
column and outputs a vcf-like file in a EMmej format
"""
import argparse
# import regex as re

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
all_args.add_argument("-af", "--allele_freq", required=True, 
      help="path to allele frequency data if available, if not, input=False (str)")


args = vars(all_args.parse_args())

vcf = pd.read_csv(args['vcf'], sep='\t' ,comment='#', header=None)
vcf = vcf.iloc[:, 0:5]
columns = {0:'#CHR', 1:'POS', 2:'ID', 3:'REF', 4:'ALT'}
vcf.rename(columns=columns, inplace=True)
vcf.drop(columns=['ID'], inplace=True)
# filter out instences were ref and alt has equal length
vcf = vcf.loc[(vcf['REF'].str.len() != vcf['ALT'].str.len()), :]
# Generating an ancestral state indication column
if args['allele_freq'] != 'False':
    columns = ['#CHR', 'POS', 'N_ALLELES', 'N_CHR', 'FREQ_REF', 'FREQ_ALT']
    af_df = pd.read_csv(args['allele_freq'], sep='\t', names=columns, skiprows=1)
    af_df['FREQ_REF'] = af_df['FREQ_REF'].str.split(pat=':', expand=True).loc[:,1].astype(float)
    af_df['FREQ_ALT'] = af_df['FREQ_ALT'].str.split(pat=':', expand=True).loc[:,1].astype(float)
    af_df['ancestral_state'] = np.nan
    af_df.loc[(af_df['FREQ_REF'] > af_df['FREQ_ALT']), 'ancestral_state'] = 0
    af_df.loc[(af_df['FREQ_ALT'] > af_df['FREQ_REF']), 'ancestral_state'] = 1
    af_df.loc[(af_df['FREQ_ALT'] == af_df['FREQ_REF']), 'ancestral_state'] = 2
    vcf['ancestral_state'] = af_df['ancestral_state'].astype(int)

vcf.sort_values(['#CHR', 'POS'], inplace=True)
vcf['#CHR'] = vcf.apply(lambda row: f"Chr{row['#CHR']}", axis=1)

vcf = vcf.loc[(~vcf['#CHR'].isin(['ChrX', 'ChrY', 'ChrMT', 'Chr4'])), :]

vcf.loc[(vcf['#CHR'] == 'Chr2L'), '#CHR'] = 'NT_033779.5'
vcf.loc[(vcf['#CHR'] == 'Chr2R'), '#CHR'] = 'NT_033778.4'
vcf.loc[(vcf['#CHR'] == 'Chr3L'), '#CHR'] = 'NT_037436.4'
vcf.loc[(vcf['#CHR'] == 'Chr3R'), '#CHR'] = 'NT_033777.3'

# filtering out wierd chromosomes
vcf.sort_values(['#CHR','POS'],inplace=True)
# filtering out deletions of 1bp
vcf['indel_len'] = vcf['ALT'].str.len() - vcf['REF'].str.len()
vcf = vcf.loc[vcf['indel_len'] != (-1), :]
vcf.drop(columns=['indel_len'], inplace=True)

path_to_output = args['outputfile']
vcf.to_csv(f'{path_to_output}/DGRP_EMmej_formated_data.vcf', sep='\t', index=False)

# python3 add_ancesral_state_column.py -v /home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/raw_data_to_EMmej_format/dgrp2_indels.vcf.gz -af /home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/raw_data_to_EMmej_format/allele_freq.frq -o /home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/raw_data_to_EMmej_format