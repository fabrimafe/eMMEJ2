# 29.06.2022

# importing libreries
import sys
sys.path.append('/home/labs/alevy/guyta/guy_master_project/src')
from io import StringIO
import datetime
from time import gmtime, strftime, localtime
import argparse
import logging

import pandas as pd
import numpy as np
from pysam import FastaFile

from pysam_getfasta import *
from data_exploration_util import *

# setting pandas display options
pd.options.display.max_colwidth = 3500
pd.set_option("display.max_columns", None)
buf, buf2, buf3 = StringIO(),StringIO(), StringIO()
date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
local_h = strftime("%H%M", localtime())
mem_reporter()

# Construct an argument parser
all_args = argparse.ArgumentParser()
# Add arguments to the parser
all_args.add_argument("-v", "--vcf", required=True,
   help="input vcf")
all_args.add_argument("-o", "--outputfile", required=True,
      help="output file in vcf-like format")
all_args.add_argument("-w", "--windowsize", required=True,
      help="size (in bp) of the context window")
all_args.add_argument("-r", "--ref", required=True,
   help="reference genome in fasta format")
all_args.add_argument("-d", "--date", required=True,
      help="date and time of job submision")
all_args.add_argument("-sim", "--simulations", required=True,
      help="a flag to indicate wether the vcf file is originated using simulations or not")
args = vars(all_args.parse_args())

fastafile=args['ref']
refFA=FastaFile(fastafile)
fname = args["vcf"].split('.')[0]

print('############################# preprocessing starts')
# setting logging config
logging.basicConfig(filename=f'{args["outputfile"]}/docs/{args["date"]}_pre_detection_processing_log.txt',
                         level=logging.DEBUG, 
                        format='%(asctime)s:%(levelname)s%(message)s')

logging.info(f'############ FILE NAME: {args["vcf"]} ############')
logging.info('-Merge VCF with reference context STARTS')
print('Merge VCF with reference context STARTS')
logging.info(mem_reporter())
# Loading data

path_to_data = f'{args["vcf"]}.vcf'
logging.info(f'-path to data: {path_to_data}')

# Loading data
Dtypes={'CHR':str, 'POS':int, 'original_pos': str,'REF':str, 'ALT':str}
col_names = ['CHR', 'POS', 'REF', 'ALT',  'original_pos']
df = pd.read_csv(path_to_data, sep="\t", header=None, dtype=Dtypes, 
                    skiprows=1, names=col_names)
df.info(buf=buf)
logging.info(f'-data loaded info: \n{buf.getvalue()}')

# making sure that there are no SNPs in the data
df = df.loc[(df.loc[:,'ALT'].str.len() != df.loc[:,'REF'].str.len()),:]
# getting rid of the 'N' in ref and NaN in ALT
df = df.loc[(df['REF'] != 'N'), :]
df = df.loc[(df['ALT'].isna() == False), :]

df.loc[:, 'ALT'] = df.loc[:, 'ALT'].str.upper()
df.loc[:, 'REF'] = df.loc[:, 'REF'].str.upper()

# defining indel type
df.loc[:, 'indel_type'] = np.nan
df.loc[(df['REF'].str.len() > df['ALT'].str.len()), 'indel_type'] = 'DEL'
df.loc[(df['REF'].str.len() < df['ALT'].str.len()), 'indel_type'] = 'INS'

# defining indel_len
df.loc[:, 'indel_len'] = np.nan
df.loc[:, 'ref_len'] = df.loc[:,'REF'].str.len()
df.loc[:, 'alt_len'] = df.loc[:,'ALT'].str.len()
df.loc[(df['indel_type'] == 'INS'), 'indel_len'] = df.loc[:, 'alt_len'] - df.loc[:, 'ref_len']
df.loc[(df['indel_type'] == 'DEL'), 'indel_len'] = df.loc[:, 'ref_len'] - df.loc[:, 'alt_len']

# isolating the indel itself (i.e. the actuall sequence that get inserted or deleted)
df.loc[:,'REF_corrected'] = np.nan
df.loc[:,'ALT_corrected'] = np.nan
df.loc[(df['indel_type'] == 'DEL'), 'REF_corrected'] = df.loc[(df['indel_type'] == 'DEL'), :].apply(lambda x: x['REF'][x['alt_len']:], axis=1)
df.loc[(df['indel_type'] == 'INS'), 'ALT_corrected'] = df.loc[(df['indel_type'] == 'INS'), :].apply(lambda x: x['ALT'][x['ref_len']:], axis=1)
df['REF'].str.slice(start=df['alt_len'])

df.drop(columns=['alt_len', 'ref_len'], inplace=True)

# Getting the reference genome around the indel
window = int(args["windowsize"])
df.loc[:, f'ref_context_seq_{window}bp'] = df.apply(lambda x: getfasta_context(refFA=refFA, chr=x['CHR'], indel_pos=x['POS'],
                     context_window_size=window), axis=1)

df.loc[:, f'ref_context_seq_{window}bp'] = df.loc[:, f'ref_context_seq_{window}bp'].str.upper()

# Reconstructing the accession context
df.loc[:, 'accession_context'] = df.apply(
    lambda x: accession_context_generator(reference_contex=x[f'ref_context_seq_{window}bp'],
                            alt=x['ALT'], ref=x['REF'], window_size=window), axis=1)
df.loc[:, 'accession_context'] = df.loc[:, 'accession_context'].str.upper()

df.drop(columns=['REF_corrected', 'ALT_corrected'], inplace=True)

# Output hendeling
path_to_output = args["outputfile"]
date = args["date"]

# saving relevant columns
col_to_save = ['CHR', 'POS', 'original_pos', 'REF', 
            'ALT', 'indel_type', 'indel_len', 
            f'ref_context_seq_{window}bp', 'accession_context']


df.loc[:, col_to_save].info(buf=buf3)
logging.info(f'-Output data info: \n{buf3.getvalue()}')
df.loc[:, col_to_save].to_csv(f'{path_to_output}/{date}_MMEJ_detection_formated_data.csv',
                                index=False, sep='\t')

logging.info('-Merge VCF with reference context DONE')
print('Merge VCF with reference context DONE')
logging.info(mem_reporter())

print('############################# preprocessing ENDS')

# python3 pre_detection_processing.py -v /home/labs/alevy/guyta/guy_master_project/data/indels_simulations/20220620_sims/mixed_files/20220629_00_MMEJ_10_NHEJ_mix.vcf -r /home/labs/alevy/guyta/guy_master_project/data/indels_simulations/from_fabrizio/indel_sumlations/simulations/GWHAAEV00000000.1.genome.fasta -w 1500 -d 20220629 -sim "True" -o "some_out"