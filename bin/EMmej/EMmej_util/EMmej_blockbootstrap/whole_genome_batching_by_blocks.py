"""
This script takes a full multi-chromosome dataset
and batch it into smaller batches with a given size.
Since this is step that one would perform befor
block-bootstraping, we refer here batches as blocks,
but practicaly they are the same.
This script was adapted from a previous script that
performs block-bootstraping.
"""
import os
import argparse
import datetime

import pandas as pd
import numpy as np

pd.options.display.max_colwidth = 3100
pd.set_option('display.max_columns', 500)

# Construct an argument parser
all_args = argparse.ArgumentParser()
# Add arguments to the parser
all_args.add_argument("-v", "--vcf", required=True,
   help="path to input vcf (string)")
all_args.add_argument("-o", "--outputfile", required=True,
      help="path to output file (string)")
all_args.add_argument("-bs", "--Blocksize", required=True,
   help="number of blocks to devide per chromosome (int)")
all_args.add_argument("-anc", "--ancestral", 
    required=False, action='store_true',default=False, 
   help="indicate whether ancestral state column is available (bool)")

args = vars(all_args.parse_args())

# set arguments
path_to_output = args['outputfile']
path_to_df = args['vcf']
blocksize = int(args['Blocksize'])
try: bootN = int(args['bootN'])
except: bootN = 1
# ancestral = True if args['ancestral'] == "True" else False

block_size = 2*10**blocksize
date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
dir = os.path.join(f"{path_to_output}/{date}_{block_size}bp_blocks")
if not os.path.exists(dir):
    os.mkdir(dir)

path_to_output = f"{path_to_output}/{date}_{block_size}bp_blocks"

# load data
if args['ancestral']: df = pd.read_csv(
   path_to_df, names=['CHR', 'POS', 'REF', 'ALT', 'ANC'], sep='\t', comment='#') # 'END',
else: df = pd.read_csv(
   path_to_df, names=['CHR', 'POS', 'REF', 'ALT'], sep='\t', comment='#') # 'END',
# df.drop(columns=['END'], inplace=True)
df.rename(columns={'CHR': '#CHR'}, inplace=True)
print(df)
MM_window_size = 1100
df['block_N'] = df.apply(lambda row: 
               int(round((row['POS']/block_size), ndigits=0)), axis=1)

# blocks_metadata = pd.DataFrame(columns=['#CHR','blockN', 'bootN', 'vcf_size'])
# # main loop
# for chrom in df['#CHR'].unique():
#    chr_df = df.loc[(df['#CHR'] == chrom), :]
#    for block in chr_df['block_N'].unique():
#       for boot in range(bootN):
#          samp_block_N = chr_df.loc[chr_df['block_N'] == block, :].sample(
#             chr_df.loc[chr_df['block_N'] == block, :].shape[0], replace=True)
         
#          fname = path_to_df.split('/')[-1].split('.')[0]
#          samp_block_N.drop(columns=['block_N'], inplace=True)
#          samp_block_N.to_csv(
#             f"{path_to_output}/{date}_block_size_{block_size}bp_chr_{chrom}_blockN{block}_sampls_bootN{boot}_{fname}.vcf",
#                   sep='\t', index=False)
#          samp_metadata = pd.DataFrame(columns=['#CHR','blockN', 'bootN', 'vcf_size'])
#          samp_metadata.loc[0,:] = chrom, block, boot, samp_block_N.shape[0]
#          blocks_metadata = pd.concat([blocks_metadata, samp_metadata])
         
# blocks_metadata.to_csv(f"{path_to_output}/{date}_blocks_metadata.tsv", sep='\t')


# # Creating README
# README = f"""This folder contains the output of: 
# /home/labs/alevy/guyta/guy_master_project/scripts/EMmej/EMmej_util/whole_genome_block_bootrapping.py
# The inputs are:
# input data: {path_to_df}
# Bootsrap parameters;
#    Each cromosome was devided to blocks with blocksize = {block_size}bp,
#    each block was then sampled {bootN} times with sample size == to the 
#    number of data points in the block, with replacment.
#    Each sampled data is saved seperatly into a VCF file.
# """
# with open(f"{path_to_output}/{date}_whole_gemome_blockbootstrap_README.txt", 'w') as f:
#    f.write(README)


blocks_metadata = pd.DataFrame(columns=['#CHR','blockN', 'vcf_size'])
# main loop
for chrom in df['#CHR'].unique():
   chr_df = df.loc[(df['#CHR'] == chrom), :]
   for block in chr_df['block_N'].unique():
      block_N_df = chr_df.loc[chr_df['block_N'] == block, :].copy()
      fname = path_to_df.split('/')[-1].split('.')[0]
      block_N_df.drop(columns=['block_N'], inplace=True)
      block_N_df.to_csv(
         f"{path_to_output}/{date}_block_size_{block_size}bp_chr_{chrom}_blockN{block}_{fname}.vcf",
               sep='\t', index=False)
      block_metadata = pd.DataFrame(columns=['#CHR','blockN', 'vcf_size'])
      block_metadata.loc[0,:] = chrom, block, block_N_df.shape[0]
      blocks_metadata = pd.concat([blocks_metadata, block_metadata])
         
blocks_metadata.to_csv(f"{path_to_output}/{date}_blocks_metadata.tsv", sep='\t')

# Creating README
README = f"""This folder contains the output of: 
/home/labs/alevy/guyta/guy_master_project/scripts/EMmej/EMmej_util/whole_genome_batching_by_blocks.py
The inputs are:
input data: {path_to_df}
Batching parameters;
   Each cromosome was devided to blocks with blocksize = {block_size}bp,
   then rach block data is saved seperatly into a VCF file.
"""
with open(f"{path_to_output}/{date}_whole_genome_batching_by_blocks_README.txt", 'w') as f:
   f.write(README)

# terminal comand example:
# python3 whole_genome_batching_by_blocks.py -v /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/stratification_by_annotation/GDL_indel_EMmej_format_non_repeated_regions_exon_no_MT.vcf.gz -o /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/stratification_by_annotation/block_sampaled_vcfs/GDL_indel_EMmej_format_non_repeated_regions_exon_no_MT -bn 5 -bs 6 -anc True