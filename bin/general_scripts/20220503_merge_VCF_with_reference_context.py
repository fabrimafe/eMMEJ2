"""
This script takes a VCF file and merge it with the corresponding
fasta context, given the rigth context coordinates
"""
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

from data_exploration_util import *

# setting pandas display options
pd.options.display.max_colwidth = 2200
pd.set_option("display.max_columns", None)
buf, buf2, buf3 = StringIO(),StringIO(), StringIO()
mem_reporter()
date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
local_h = strftime("%H%M", localtime())

def add_fasta_context_to_dataframe(context_file: str, context_window: int):
    """
    This function opens the output file of the bedtools getfasta
    function and output it as a dataframe
    Args:
    returns:
    """
    with open(context_file, 'r') as fasta_context: # loading the file
        # create a dataframe out of the file
        context_df = pd.read_csv(fasta_context, header=None)
        # take only the rows that coresponse to 
        # context location in the fasta reference
        fasta_context_position = np.array(context_df.iloc[::2, 0])
        # take only the rows that coresponse to context sequance in the fasta reference
        fasta_context_seq = np.array(context_df.iloc[1::2, 0]) 
        # transforming into 1D vector
        fasta_context_position = pd.Series(fasta_context_position) 
        # transforming into 1D vector
        fasta_context_seq = pd.Series(fasta_context_seq)        
        return (pd.DataFrame({'ref_context_position': fasta_context_position,
                                   'ref_context' : fasta_context_seq}))


def accession_context_generator(reference_contex : str,
                                alt : str, ref : str, window_size: int):
    # matching window_size to be 0 indexed
    window_size = window_size-1
    # # Chacking that tha reference context match the REF
    reference_contex_REF = reference_contex[
        window_size:(len(ref)+window_size)]

    assert (reference_contex_REF == ref), f'REF dose not match the reference context, expected: {ref}, got: {reference_contex_REF}, Indel position (0 based): {window_size}, reference_contex: {reference_contex}'
    # Generating accession context
    per_indel_seq = reference_contex[:(window_size)]
    post_indel_seq = reference_contex[(window_size + len(ref)):]
    
    # Chacking that accsession context was is correct
    accession_context = f'{per_indel_seq}{alt}{post_indel_seq}'
    accession_context_check = accession_context[
        window_size:(len(alt)+window_size)]
    assert (accession_context_check == alt), f'accession_context_check dose not match ALT, expected: {alt}, got: {accession_context_check}'
    
    return accession_context

# Construct an argument parser
all_args = argparse.ArgumentParser()
# Add arguments to the parser
all_args.add_argument("-v", "--vcf", required=True,
   help="input vcf")
all_args.add_argument("-o", "--outputfile", required=True,
      help="output file in vcf-like format")
all_args.add_argument("-w", "--windowsize", required=True,
      help="size (in bp) of the context window")
all_args.add_argument("-f", "--fastacontext", required=True,
      help="a bed file with the fasta context")
all_args.add_argument("-d", "--date", required=True,
      help="date and time of job submision")
all_args.add_argument("-ra", "--realaigned", required=True,
      help="a flag to indicate wether the vcf file went through re-alaingment process")
all_args.add_argument("-sim", "--simulations", required=True,
      help="a flag to indicate wether the vcf file is originated using simulations or not")


args = vars(all_args.parse_args())

fname = args["vcf"].split('.')[0]
# setting logging config
logging.basicConfig(filename=f'{args["outputfile"]}/{args["date"]}_merge_VCF_with_reference_context_log.txt',
                         level=logging.DEBUG, 
                        format='%(asctime)s:%(levelname)s%(message)s')

logging.info(f'############ FILE NAME: {args["vcf"]} ############')
logging.info('Merge VCF with reference context STARTS')
print('Merge VCF with reference context STARTS')
logging.info(mem_reporter())
# Loading data
path_to_data = f'{args["vcf"]}.vcf'
logging.info(f'path to data: {path_to_data}')

if args["simulations"]:
    Dtypes={'CHR':str, 'POS':int, 'REF':str, 'ALT':str, 'ground_true': str}
    df = pd.read_csv(path_to_data, sep="\t", header=None, dtype=Dtypes, 
                skiprows=1, names=['CHR', 'POS', 'REF', 'ALT', 'ground_true', 'original_pos'])    
else:
    Dtypes={'CHR':str, 'POS':int, 'REF':str, 'ALT':str, }
    df = pd.read_csv(path_to_data, sep="\t", header=None, dtype=Dtypes, 
                    skiprows=1, names=['CHR', 'POS', 'REF', 'ALT', 'original_pos'])


df.info(buf=buf)
logging.info(f'data loaded info: \n{buf.getvalue()}')
# making sure that there are no SNPs in the data
df = df.loc[(df.loc[:,'ALT'].str.len() != df.loc[:,'REF'].str.len()),:]
# getting rid of the 'N' in ref and NaN in ALT
df = df.loc[(df['REF'] != 'N'), :]
df = df.loc[(df['ALT'].isna() == False), :]

df.loc[:, 'ALT'] = df.loc[:, 'ALT'].str.upper()
df.loc[:, 'REF'] = df.loc[:, 'REF'].str.upper()

window = int(args["windowsize"])
df.loc[:, f'ref_chr_format_{window}bp'] = df.apply(
    lambda x: (f">{x['CHR']}:{x['POS'] - window}-{x['POS'] + window}"), axis=1)

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

context_1500 = args["fastacontext"]
context_1500bp_df = add_fasta_context_to_dataframe(context_1500, window)
context_1500bp_df.rename(columns=
                        {'ref_context_position': f'ref_chr_format_{window}bp', 
                         'ref_context' : f'ref_context_seq_{window}bp'}, 
                       inplace=True)
context_1500bp_df.info(buf=buf2)
logging.info(f'context_1500bp_df info: \n{buf2.getvalue()}')

context_1500bp_df.drop_duplicates(inplace=True)


df = pd.merge(df,context_1500bp_df, on=[f'ref_chr_format_{window}bp'], how='inner')
df.loc[:, f'ref_context_seq_{window}bp'] = df.loc[:, f'ref_context_seq_{window}bp'].str.upper()

df.loc[:, 'accession_context'] = df.apply(
    lambda x: accession_context_generator(reference_contex=x[f'ref_context_seq_{window}bp'],
                            alt=x['ALT'], ref=x['REF'], window_size=window), axis=1)
df.loc[:, 'accession_context'] = df.loc[:, 'accession_context'].str.upper()

path_to_output = args["outputfile"]
date = args["date"]
# saving relevant columns

if args["simulations"]:
    col_to_save = ['CHR', 'POS', 'ground_true', 'REF', 
            'ALT', 'indel_type', 'indel_len', 
            f'ref_context_seq_{window}bp', 'accession_context']
else:
    col_to_save = ['CHR', 'POS', 'REF', 
            'ALT', 'indel_type', 'indel_len', 
            f'ref_context_seq_{window}bp', 'accession_context']

df.loc[:, col_to_save].info(buf=buf3)
logging.info(f'Output data info: \n{buf3.getvalue()}')
df.loc[:, col_to_save].to_csv(f'{path_to_output}/{date}_MMEJ_detection_formated_data.csv',
                                index=False, sep='\t')

logging.info('Merge VCF with reference context DONE')
print('Merge VCF with reference context DONE')
logging.info(mem_reporter())
