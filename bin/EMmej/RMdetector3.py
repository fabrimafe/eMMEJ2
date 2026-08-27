
# python3 RMdetector.py -v /home/labs/alevy/guyta/guy_master_project/data/indels_simulations/20220620_sims/mixed_files/20220629_07_MMEJ_03_NHEJ_mix_unlabled.vcf -o /home/labs/alevy/guyta/guy_master_project/scripts/general_scripts/RM_wraped_software -w 1500 -ra True -anc False -r /home/labs/alevy/guyta/guy_master_project/data/indels_simulations/from_fabrizio/indel_sumlations/simulations/GWHAAEV00000000.1.genome.fasta

from ast import arg
import sys
sys.path.append('src')
sys.path.append('/home/fabio/git/eMMEJ2/src')
from io import StringIO
import string
import datetime
from time import gmtime, strftime, localtime
import argparse
import logging
import regex as re

import pandas as pd
import numpy as np
from pysam import FastaFile
import Bio
from Bio import Align #from Bio import pairwise2. # pairwise2 is now deprecated in Biopython >1.80
from Bio import pairwise2

from MicroHomology_module_v5 import emMEJrealignment
from realignment_module import * 
from pysam_getfasta import *

# setting pandas display options
# pd.options.display.max_colwidth = 3500
pd.set_option("display.max_columns", None)
buf, buf2, buf3 = StringIO(),StringIO(), StringIO()
date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
local_h = strftime("%H%M", localtime())

# Construct an argument parser
all_args = argparse.ArgumentParser()
# Add arguments to the parser
all_args.add_argument("-v", "--vcf", required=True,
   help="path to input vcf (string)")
all_args.add_argument("-o", "--outputfile", required=False,
      help="path to output file (string)", default="RMdetector_output.tsv")
all_args.add_argument("-r", "--ref", required=True,
   help="path to reference genome in fasta format (string)")
all_args.add_argument("-w", "--windowsize", required=False,
     default=150, type=int,
      help="size (in bp) of the context window (int)")

# optional arguments, turned on by using the flag
all_args.add_argument("-mhl", "--MH_lengths", required=False ,default=2, 
      help="a flag that specify the specific MH length that one wants to analyze, default is all MH lengths")
all_args.add_argument("-cr", "--CRISPR", required=False, action='store_true',default=False, 
      help="a flag to indicate if the algorithm performes on CRISPR data")
all_args.add_argument("-anc", "--ancestral", required=False, action='store_true',default=False, 
      help="a flag to indicate whether an ancestral state is available or not")
all_args.add_argument("-ic", "--include_context", required=False,  default=False, 
      action='store_true',
      help="a flag to indicate whether to include context in output or not")
all_args.add_argument("-vr", "--verbose", required=False, 
        action='store_true',default=False, 
      help="a flag to turn on verbose mode")
all_args.add_argument("-ex", "--extension", required=False,
        action='store_true',default=False,
      help="a flag to turn on the extension of MH and SD")

args = vars(all_args.parse_args())

# ------------------ PREPARING DATA ------------------------------------
#
# In this section
# 1) cleaning of the data based on the xo
# 2) define ancestral state
# 3) if genomic data realign
#
# 
# Loading data

MH_lengths=args['MH_lengths'].split(",")
MH_lengths = [int(x) for x in MH_lengths]

path_to_data = args["vcf"]
Dtypes={'CHR':str, 'POS':int,'variant_id':str, 'REF':str, 'ALT':str, 'ANCESTRAL': str}
col_names = ['CHR', 'POS','variant_id', 'REF', 'ALT', 'ANCESTRAL']

if not args['ancestral']: col_names.remove('ANCESTRAL')
df = pd.read_csv(path_to_data, sep="\t", header=None, dtype=Dtypes, 
                    names=col_names, comment='#') 

# filtering out cases in which length of REF and ALT are equal

#df = df.loc[(df['REF'].str.len() != df['ALT'].str.len()), :]
if args['verbose']: print(f"Input data dimentions:{df.shape}")
if args['ancestral']: set_ancestral_state_indel(df=df)
else: df.rename(columns={'REF': 'ANC', 'ALT': 'DER'}, inplace=True)
#  reading fasta file using pysam
fastafile=args['ref']
refFA=FastaFile(fastafile)

# Filtering out indels with non-ACGT characters (e.g. N) in their context.
df = df.loc[(df['POS'] > 1000), :]

problematic_nuc = list(string.ascii_uppercase)
problematic_nuc.remove('A')
problematic_nuc.remove('C')
problematic_nuc.remove('G')
problematic_nuc.remove('T')


df['context_contains_N'] = df.apply(lambda row: any(c in str(refFA.fetch(row['CHR'],
            row['POS']-1000 ,row['POS']+1000)).upper() for c in problematic_nuc), axis=1)

if args['verbose']:
    print(
        f"""Filtering {df.loc[df['context_contains_N'], :].shape[0]} indels that contain N in their reference context""")

df = df.loc[(df['context_contains_N'] == False), :]
df.drop(columns=['context_contains_N'], inplace=True)
print(df[['DER', 'ANC']])
# df = df.loc[df['CHR'] == 'pol_slip_monodirectional_0',:]

#if genomic data realign indels
SUB = (
        (df['ANC'].str.len() >= 1) & 
        (df['DER'].str.len() >= 1) &
        (df['ANC'].str[0] != df['DER'].str[0])).any()

if not args['CRISPR'] and not SUB:
    df=df.apply(lambda row : vcf2realignedvcfs_pairwise2(refFA,row['CHR'],row['POS'], row['ANC'],row['DER'],150), axis = 1)
    df=flatten_2list(df.tolist())
    df=pd.DataFrame(df,columns = ['CHR','POS','ANC',"DER","original_pos"])
    df.loc[:, 'POS'] = df.loc[:, 'POS'] + 1
else:
    df['original_pos'] = df['POS']

    
"""
#now this new lines allow to work only with the variant on the original vfc (input)
if not args['CRISPR']:
    orig_ANC = df['ANC'].iloc[0]
    orig_DER = df['DER'].iloc[0]

    df = df.apply(lambda row : vcf2realignedvcfs_pairwise2(refFA,row['CHR'],row['POS'],row['ANC'],row['DER'],150),axis = 1)

df = flatten_2list(df.tolist())

df = pd.DataFrame(df,columns=['CHR', 'POS', 'ANC', 'DER', 'original_pos'])
df = df[(df['ANC'] == orig_ANC) & (df['DER'] == orig_DER)]
"""

# -------------------------- PREPARING THE DATA TO MOTIF FINDING STEP -------------------------
# making sure that there are no SNPs in the data
#df = df.loc[(df.loc[:,'DER'].str.len() != df.loc[:,'ANC'].str.len()),:]


print(df[['DER', 'ANC']])
# Removing duplicates --------------------------
df.drop_duplicates(subset=['CHR','POS','ANC',"DER"],inplace=True)
df.reset_index(drop=True,inplace=True)
print(df[['DER', 'ANC']])  
# getting rid of the 'N' in ref and NaN in DER
if args['verbose']:
    print(
        f"Filtering {df.loc[df['ANC'].str.contains(pat='N'), :].shape[0]} indels that contain N in their indel sequence")
df = df.loc[(df['ANC'].str.contains(pat='N') == False), :]
df = df.loc[(df['DER'].isna() == False), :]

# create unique variant name for each original indel
def create_variant_id(chrom, original_pos):
    """Define name (id) for a set of potential realignments of the same original indel
    """
    return(chrom+"_"+str(original_pos).split(".")[0])

if args['CRISPR']: df['original_pos'] = df['POS']
if df['original_pos'].isna().any():
    df['variant_id'] = df.apply(lambda x: create_variant_id(x['CHR'], (x['POS'])), axis=1)
else:
    df['variant_id'] = df.apply(lambda x: create_variant_id(x['CHR'], (x['original_pos'])), axis=1)

df.loc[:, 'DER'] = df.loc[:, 'DER'].str.upper()
df.loc[:, 'ANC'] = df.loc[:, 'ANC'].str.upper()


# defining indel type and length
df.loc[:, 'indel_type'] = np.nan
df.loc[
        (df['ANC'].str.len() > 1) & 
        (df['DER'].str.len() == 1) &
        (df['ANC'].str[0] == df['DER']),
         'indel_type'] = 'DEL'



df.loc[
        (df['ANC'].str.len() == 1) &
        (df['DER'].str.len() > 1) &
        (df['ANC'] == df['DER'].str[0]),
        'indel_type'] = 'INS'

#df.loc[(df['ANC'].str.len() < df['DER'].str.len()), 'indel_type'] = 'INS'

df.loc[
        (df['ANC'].str.len() >= 1) &
        (df['DER'].str.len() >= 1) &
        (df['ANC'].str[0] != df['DER'].str[0]),
        'indel_type'] = 'SUB'


df.loc[:, 'indel_len'] = np.nan
df.loc[:, 'ref_len'] = df.loc[:,'ANC'].str.len()
df.loc[:, 'alt_len'] = df.loc[:,'DER'].str.len()
df.loc[(df['indel_type'] == 'INS'), 'indel_len'] = df.loc[:, 'alt_len'] - df.loc[:, 'ref_len']
df.loc[(df['indel_type'] == 'DEL'), 'indel_len'] = df.loc[:, 'ref_len'] - df.loc[:, 'alt_len']

# --------------- Repair mechanism detection ---------------------------------------------

indel_position = args['windowsize']

#Define each row as a emMEJrealignment class, applyng class to every row (realignment) of dataset, which is define in MicroHomology_module_v3. In this class, methods for detection and manipulation.


#First step, run 5'->3' (flip=False)---------->>>>>>>>>>
microhomology_object = df.apply( lambda x: emMEJrealignment(
                            DER=x['DER'], # change name to derived seq
                            ANC=x['ANC'],
                            pos_on_chr=x['POS'],
                            indel_type=x['indel_type'],
                            flip=False,
                            include_context=args["include_context"],
                            extension=args["extension"],
                            MH_lengths=MH_lengths,
                            windowsize=args['windowsize'], 
                            refFA=refFA, chrom=x['CHR']), axis = 1)
print("detection done")
#create a dataframe (df2) by appending all the exported data (ex_data) from microhomology_object
df2 = None
for i in microhomology_object.index:
    if (i == 0):
        df2 = microhomology_object[i].ex_data
    else:
        df2 = pd.concat([df2, microhomology_object[i].ex_data])


print(df2.columns)       
df2.reset_index(drop=True, inplace=True)
df.reset_index(drop=True, inplace=True)
df = df.join(df2, how = 'right')
df['direction'] = 0

print("appending done")
#2nd step, run on flipped seq (3'->5'):--------->>>>>>>>
#if crispr, look for ALL patterns ALSO on the flipped seq:-->
if args['CRISPR']:
    df_fl = df.loc[:, ['CHR', 'POS', 
        'ANC', 'DER', 'original_pos',
        'variant_id', 'indel_type', 'indel_len' ]].copy() 

    microhomology_object = df_fl.apply( lambda x: emMEJrealignment(
                            DER=x['DER'],
                            ANC=x['ANC'],
                            pos_on_chr=x['POS'],
                            indel_type=x['indel_type'],
                            flip=True,
                            include_context=args["include_context"],
                            extension=args["extension"],
                            MH_lengths=MH_lengths,
                            windowsize=args['windowsize'], 
                            refFA=refFA, chrom=x['CHR']), axis = 1)

    tmp_df = None
    for i in microhomology_object.index:
        if (i == 0):
            tmp_df = microhomology_object[i].ex_data
        else:
            tmp_df = pd.concat([tmp_df, microhomology_object[i].ex_data])
            
    tmp_df.reset_index(drop=True, inplace=True)
    df_fl.reset_index(drop=True, inplace=True)
    df_fl = df_fl.join(tmp_df, how='right')
    df_fl['direction'] = 1
    df = pd.concat([df, df_fl], axis=0)
    # in the case of crispr, make variable ID uniq per
    # direction, this is necessary for proper EM functioning
    df['variant_id'] = df.apply(
        lambda row: f"{row['variant_id']}.{row['direction']}", axis=1)
    
# if not crispr, look for patterns (only for insertions, no deletion) on the fliped seq (3'->5'):--> 
"""
print("start insertions")
if not args['CRISPR']:
    df_ins = df.loc[(df['indel_type'] == 'INS'), ['CHR', 'POS', 
        'ANC', 'DER', 'original_pos',
        'variant_id', 'indel_type', 'indel_len']].copy()
    if df_ins.shape[0] > 0:
        microhomology_object = df_ins.apply( lambda x: emMEJrealignment(
                                DER=x['DER'],
                                ANC=x['ANC'],
                                pos_on_chr=x['POS'],
                                indel_type=x['indel_type'],
                                flip=True,
                                include_context=args["include_context"],
                                extension=args["extension"],
                                MH_lengths=MH_lengths,
                                windowsize=args['windowsize'], 
                                refFA=refFA, chrom=x['CHR']), axis = 1)
        
        tmp_df = None
        for i in microhomology_object.index:
            if (i == 0):
                tmp_df = microhomology_object[i].ex_data
            else:
                tmp_df = pd.concat([tmp_df, microhomology_object[i].ex_data])
                
        tmp_df.reset_index(drop=True, inplace=True)
        df_ins.reset_index(drop=True, inplace=True)
        df_ins = df_ins.join(tmp_df, how='right')
        df_ins['direction'] = 1
        df = pd.concat([df, df_ins], axis=0)

df.sort_values(by=['CHR', 'POS'], axis=0, ascending=True, inplace=True)
df.reset_index(drop=True, inplace=True)
"""
print("save results")
#Use logging module to have a nice log output
logging.debug(f'## MMEJ detection is done({args["vcf"]})')
df.info(buf=buf)
logging.info(buf.getvalue())

#Rename some output variables
#df['snap_mh2_len'] = df['snap_mh2'].str.len()
#df['loop_mh2_len'] = df['loop_mh2'].str.len()

#df.loc[(df['del_mmej'] == True), 'del_motif_d'] = df.loc[(df['del_mmej'] == True), 'indel_len'] - df.loc[(df['del_mmej'] == True), 'del_mmej_cand_len'].apply(lambda x: x[0]) #useful for markov model. but easier to do it later, after I reconvert string into lists

#Create output table of positions. Unfortunately, it works only for one pattern for mechanism.
"""
mechanism = ['del_mmej', 'SD_loop_out', 'SD_snap_back']
patterns = ['del_mmej_cand', 'loop_repeat_pat', 'snap_repeat_pat']
cols = [['del_mmej_motif_pos', 'del_mmej_freq_small_window', 'del_mmej_freq_large_window'],
        ['loop_motif_pos', 'loop_freq_small_window', 'loop_freq_large_window'],
        ['snap_motif_pos', 'snap_freq_small_window', 'snap_freq_large_window']]
for mech,pat,col in zip(mechanism, patterns, cols):
    df.loc[:, col] = np.nan, np.nan, np.nan
    tmp_df = df.loc[
            (df[mech] == True),:].apply(lambda row: get_motifs_freqs(
                ref=refFA, CHR=row['CHR'], POS=row['POS'], large_window=1000,
                small_window=args['windowsize'],
                    motifs=row[pat], indel_type=row['indel_type']),
                    axis=1, result_type='expand')
    print("1done")
    tmp_df.rename(columns={0:col[0],1:col[1],2:col[2]},inplace=True)
    df.loc[
        (df[mech] == True),col] = tmp_df
"""

# ----------- Make this more efficient -------------------------------------

df['snap_repeat_pat_len'] = np.nan
if 'SD_inverted_insertion' in df.columns:
    df.loc[df['SD_inverted_insertion']==True,'SD_II_repeat_pat_len'] = df.loc[df['SD_inverted_insertion']==True, 'SD_II_repeat_pat'].str.len()
#if 'SD_loop_out' in df.columns:
    #df['loop_repeat_pat_len'] = np.nan
#df.loc[df['SD_loop_out']==True,'loop_repeat_pat_len'] = df.loc[df['SD_loop_out']==True, 'loop_repeat_pat'].str.len()

if 'SD_direct_insertion' in df.columns:
    df['SD_DI_repeat_pat_len'] = np.nan
df.loc[df['SD_direct_insertion']==True,'SD_DI_repeat_pat_len'] = df.loc[df['SD_direct_insertion']==True, 'SD_DI_repeat_pat'].str.len()

col_to_save = ['CHR', 'POS', 'original_pos', 'variant_id', 'direction', #'REF','ALT',
            'ANC','DER', 'indel_type', 'indel_len',  
            # deletions
            'del_mmej', 'del_mmejl','del_mmej_cand', 'del_mmej_marked_on_ref', 'del_mmej_marked',
            'del_last_dimer','del_mmej_cand_len',
            #del_mmej motif
            'del_mmej_motif_pos','del_mmej_freq_small', 'del_mmej_freq_large',
                      
            #SD-direct_deletion
            'SD_direct_deletion','SD_DD_mutant_pattern','SD_DD_MH2','SD_DD_P2',
            'SD_DD_repeat_pat','SD_DD_repeat_pat_len','SD_DD_last_dimer','SD_DD_dist_between_reps',
            'SD_DD_motif_pos','SD_DD_motif_freq_small','SD_DD_motif_freq_large',
            
            # SD_inverted_deletions
            'SD_inverted_deletion','SD_ID_mutant_pattern','SD_ID_MHrevc','SD_ID_Prevc',
            'SD_ID_repeat_pat','SD_ID_repeat_pat_len','SD_ID_last_dimer','SD_ID_dist_between_reps',
            'SD_ID_motif_pos','SD_ID_motif_freq_small','SD_ID_motif_freq_large',

            # SD_direct_insertion
            'SD_direct_insertion','SD_DI_mutant_pattern', 'SD_DI_MH2',
            'SD_DI_P2','SD_DI_repeat_pat', 'SD_DI_repeat_pat_len',
            'SD_DI_last_dimer', 'SD_DI_dist_between_reps','SD_DI_motif_pos',
            'SD_DI_motif_freq_small','SD_DI_motif_freq_large',
            
            #SD_inverted_insertion
            'SD_inverted_insertion','SD_II_mutant_pattern','SD_II_MHrevc','SD_II_Prevc',
            'SD_II_repeat_pat','SD_II_repeat_pat_len','SD_II_last_dimer','SD_II_dist_between_reps',
            'SD_II_motif_pos','SD_II_motif_freq_small','SD_II_motif_freq_large',
            
            #SD_direct_substitution
            'SD_direct_substitution','SD_DS_mutant_pattern','SD_DS_P2','SD_DS_MH2',
            'SD_DS_repeat_pat','SD_DS_repeat_pat_len','SD_DS_dist_between_reps','SD_DS_last_dimer',
            'SD_DS_motif_pos','SD_DS_motif_freq_small','SD_DS_motif_freq_large',
            
            #SD_inverted_substitution
            'SD_inverted_substitution','SD_IS_mutant_pattern','SD_IS_Prevc','SD_IS_MHrevc',
            'SD_IS_repeat_pat','SD_IS_repeat_pat_len','SD_IS_dist_between_reps','SD_IS_last_dimer',
            'SD_IS_motif_pos','SD_IS_motif_freq_small','SD_IS_motif_freq_large',
                               
            # polymerase slippage
            'pol_slip', 'pol_slip_submotif', 'pol_slippage_repeatsIndel',
            'pol_slippage_repeatsDownstream', 'pol_slip_pos'] # 'pol_slippage_times']

if args["include_context"]: col_to_save = col_to_save + ['ref_genome_context', 'mutant_sequence']

for c in col_to_save: #se una colonna non viene creato, come nel caso di plo_slip per le sub allora è falsa
    if c not in df.columns:
        df[c] = False

df.loc[:,'del_mmej'].fillna(value=False, inplace=True)
df.loc[:,'SD_inverted_insertion'].fillna(value=False, inplace=True)
df.loc[:,'SD_direct_insertion'].fillna(value=False, inplace=True)
df.loc[:,'pol_slip'].fillna(value=False, inplace=True)
df.loc[:,'SD_direct_deletion'].fillna(value=False, inplace=True)
df.loc[:,'SD_inverted_deletion'].fillna(value=False, inplace=True)
df.loc[:,'SD_direct_substitution'].fillna(value=False, inplace=True)
df.loc[:,'SD_inverted_substitution'].fillna(value=False, inplace=True)
df.loc[:, col_to_save].to_csv(f"{args['outputfile']}", sep='\t', index=False)
