# This script detect MMEJ in a given dataset
# importing the relevant libraries/modules:
import datetime
import sys
from io import StringIO
import logging
sys.path.append('/home/labs/alevy/guyta/guy_master_project/src')
import argparse
import regex as re
import pandas as pd
import numpy as np

# from VcfProcess_module import VcfProcess
from data_exploration_util import *
from MicroHomology_module_v3 import MicroHomology
from MMEJ_2nd_order_MM_v2 import main_markovian_process
from MMEJ_2nd_order_MM_v2_old import  motif_probabily_calc
from ExpectationMaximization import ExpectationMaximization
from ExpectationMaximization_q import EMq
pd.options.display.max_colwidth = 3100
pd.set_option('display.max_columns', 500)


buf, buf2 = StringIO(), StringIO()
# Construct an argument parser
all_args = argparse.ArgumentParser()
# Add arguments to the parser
all_args.add_argument("-i", "--inputpath", required=True,
   help="input vcf")
all_args.add_argument("-o", "--outputfile", required=True,
      help="output file in vcf-like format")
all_args.add_argument("-w", "--windowsize", required=True,
      help="size (in bp) of the context window")
all_args.add_argument("-f", "--inputfilename", required=True,
      help="a bed file with the fasta context")
all_args.add_argument("-d", "--date", required=True,
      help="date and time of job submision")
all_args.add_argument("-sim", "--simulations", required=True,
      help="a flag to indicate wether the vcf file is originated using simulations or not")

args = vars(all_args.parse_args())
fname = args["inputfilename"].split('.')[0]
# setting logging config
logging.basicConfig(filename=f'{args["outputfile"]}/docs/{args["date"]}_{fname}_MMEJ_detection_pipline_log.txt',
                         level=logging.DEBUG, 
                        format='%(asctime)s:%(levelname)s%(message)s')

logging.info(f'############ FILE NAME: {args["inputfilename"]} ############')
logging.info('MMEJ detection pipline STARTS')
print('MMEJ detection pipline STARTS')
# logging.info(mem_reporter())
# Loading the data:
# Specify Dtypes
context_window = int(args["windowsize"])

# Loading the data
dtypes = {'CHR':str, 'POS':int, 'original_pos':str, 'REF':str, 
            'ALT':str, 'indel_type':str, 'indel_len':int, 
            f'ref_context_seq_{context_window}bp':str, 'accession_context':str}

# set data path
path_to_df = args["inputpath"]
# loading data
df = pd.read_csv(path_to_df, sep = '\t', dtype = dtypes)
# df = df.iloc[0:50, :]

"""
Detecting microhomology:
The data must have the following columns:
# accession_context -> A sequence with the indel and no gaps (type: str)
# derived_allele -> The indel itself (i.e. the sequence that got inserted/deleted)
# ancesstral_allele -> The reference sequence at the indel position 
    (in a case of insertion the value should be: "", in case of deletion 
    the value should be: the deleted sequence)
# reference_context_seq_120 -> the reference genome at a window of +-120bp from the indel position
    (i.e. -120bp ->INDEL->+ 120bp)
# indel_pos -> the indel position in the read
# indel_len -> the length of the indel
# indel_type -> the type of the indel (INS\DEL)
"""
# Creating MicroHomology objects (from MicroHomology_module_v3)
logging.debug(f'MMEJ detection starts for {args["inputfilename"]}')
indel_position = context_window # - 1

microhomology_object = df.apply( lambda x: MicroHomology(
                            accession_sequence=x['accession_context'],
                            indel_sequence=x['ALT'],
                            ancestral_sequence=x['REF'],
                            ref_sequence=x[f'ref_context_seq_{context_window}bp'],
                            indel_position=indel_position, # converting 1 based to 0 based
                            indel_type=x['indel_type']), axis = 1)

df2 = None
for i in microhomology_object.index:
    if (i == 0):
        df2 = microhomology_object[i].ex_data
    else:
        df2 = pd.concat([df2, microhomology_object[i].ex_data])
        
df2.reset_index(drop = True, inplace = True)

df = df.join(df2, how = 'right')

logging.debug(f'## MMEJ detection is done({args["inputfilename"]})')
df.info(buf=buf)
logging.info(buf.getvalue())

df['snap_mh2_len'] = df['snap_mh2'].str.len()
df['loop_mh2_len'] = df['loop_mh2'].str.len()

# 17.06.2022 -> trying to change the motif_d input for the Markov model for del MMEJ
df.loc[(df['del_mmej'] == True), 'del_motif_d'] = df.loc[(df['del_mmej'] == True), 'indel_len'] - df.loc[(df['del_mmej'] == True), 'del_mmej_cand_len']
"""
Calculating liklyhood and p-value of finding a motif in a
specific distance from the 1st motif.
"""
df['del_mmej_lk'] = np.nan
# df['trans_mmej_lk'] = np.nan
df['SD_snap_back_lk'] = np.nan
df['SD_loop_out_lk'] = np.nan
df['del_mmej_p_val'] = np.nan
# df['trans_mmej_p_val'] = np.nan
df['SD_snap_back_p_val'] = np.nan
df['SD_loop_out_p_val'] = np.nan
df['del_mmej_prob'] = np.nan
# df['trans_mmej_prob'] = np.nan
df['SD_snap_back_prob'] = np.nan
df['SD_loop_out_prob'] = np.nan
df['del_mmej_dist_CDF'] = np.nan
# df['trans_dist_CDF'] = np.nan
df['SD_snap_back_dist_CDF'] = np.nan
df['SD_loop_out_dist_CDF'] = np.nan
df['del_mmej_q'] = 0


"""
for insertions: the motif is the repeat ('tamplate_switch_repeat')
for  mmej snap back: the motif is the repeat ('S_tamplate_switch_repeat')
for deletions: the motif is the 'mmej_cand'
for all indels the last nucleotide is the one that comes right before the DSB on the reference seq
"""


logging.debug('DEL mmej Markov starts')


# ---------------------------------------------------------------------------------------------------
df['del_mmej_old_prob'] = np.nan
df['NHEJ_old_prob'] = np.nan


# creating a df with only deletions
df_with2K_del_mmej = df.loc[((df['del_mmej'] == True) & 
           (df['indel_type'] == 'DEL')), :].copy()


MM_del_df = df_with2K_del_mmej.apply(
        lambda x: motif_probabily_calc(
            _sequance = x[f'ref_context_seq_{context_window}bp'] , 
                                 _motif = x['del_mmej_cand'],
                                _indel_len = int(x['del_motif_d']), # 17.06.2022 motifs_d = int(x['indel_len']),
                    _memory_dimer = x['del_last_dimer']), axis = 1,result_type='expand')

df.loc[ df_with2K_del_mmej.index, 'del_mmej_old_prob'] = MM_del_df

df['del_mmej_old_prob'] = 1- df['del_mmej_old_prob']
df['del_mmej_old_prob'].fillna(value=0, inplace=True)

df.loc[:,'NHEJ_old_prob'] = 0
df.loc[(df['indel_type'] == 'DEL') ,'NHEJ_old_prob'] = 1-df.loc[(df['indel_type'] == 'DEL'),'del_mmej_old_prob']
# ---------------------------------------------------------------------------------------------------

# making the new version blind to 1bp MH length
# df_with2K_del_mmej = df.loc[((df['del_mmej'] == True) & (df['indel_type'] == 'DEL') & (df['del_mmej_cand_len'] != 1)), :].copy()

MM_del_df = df_with2K_del_mmej.apply(
        lambda x: main_markovian_process(
            sequance = x[f'ref_context_seq_{context_window}bp'] , 
                                 motif = x['del_mmej_cand'],
                                motifs_d = int(x['del_motif_d']), # 17.06.2022 motifs_d = int(x['indel_len']),
                    memory_dimer = x['del_last_dimer'],
                    early_stop=1), axis = 1,result_type='expand')
MM_del_df.rename(columns={0: 'del_mmej_prob', 1:'del_mmej_lk',
                2:'del_mmej_p_val', 3:'del_mmej_dist_CDF', 4: 'del_mmej_q'}, inplace=True)
df.loc[ df_with2K_del_mmej.index, 
    ['del_mmej_prob','del_mmej_lk', 'del_mmej_p_val', 'del_mmej_dist_CDF', 'del_mmej_q']] = MM_del_df.loc[:, 
                    ['del_mmej_prob','del_mmej_lk', 'del_mmej_p_val', 'del_mmej_dist_CDF', 'del_mmej_q']]
df.loc[:,'del_mmej_motif_count'] = np.nan
df.loc[df_with2K_del_mmej.index,'del_mmej_motif_count'] = df.loc[df_with2K_del_mmej.index,: ].apply(lambda x: 
    len([m.start() for m in re.finditer(x['del_mmej_cand'], x[f'ref_context_seq_{context_window}bp'], overlapped=True)]), axis = 1)

logging.debug('DEL mmej Markov ends')

# creating a df with only trans insertions
# df_with2K_ins_mmej = df[(df['SD_trans'] == True)].copy()
# df_with2K_ins_mmej.loc[:,'trans_mmej_cand_len'] = df_with2K_ins_mmej['trans_mmej_cand_len'].astype('int')
# logging.debug('Trans mmej likelihood (Markov) starts')
# MM_trans_df = df_with2K_ins_mmej.apply(
#         lambda x: main_markovian_process(sequance = x[f'ref_context_seq_{context_window}bp'] , 
#                                  motif = x['trans_TS_repeat'],
#                                 motifs_d = x['trans_mmej_cand_len'],
#                     memory_dimer = x['trans_last_dimer'],
#                     early_stop=1), axis = 1 ,result_type='expand')
# MM_trans_df.rename(columns={0: 'trans_mmej_prob', 1:'trans_mmej_lk',
#                 2:'trans_mmej_p_val', 3: 'trans_dist_CDF', 4: 'trans_q'}, inplace=True)

# df.loc[:,'trans_dist_CDF'] = np.nan
# df.loc[ df_with2K_ins_mmej.index, 
#             ['trans_mmej_prob','trans_mmej_lk','trans_mmej_p_val', 'trans_dist_CDF']] = MM_trans_df.loc[:, 
#                                             ['trans_mmej_prob','trans_mmej_lk','trans_mmej_p_val', 'trans_dist_CDF']]
# df.loc[:,'trans_motif_count'] = np.nan
# df.loc[df_with2K_ins_mmej.index,'trans_motif_count'] = df.loc[df_with2K_ins_mmej.index,: ].apply(lambda x: 
#     len([m.start() for m in re.finditer(x['trans_TS_repeat'], x[f'ref_context_seq_{context_window}bp'], overlapped=True)]), axis = 1)


# logging.debug('Trans mmej Markov ends')

logging.debug('Snap back mmej Markov starts')
# creating a df with only mmej snap back
df_with2K_mmej_snap_back = df[df['SD_snap_back'] == True].copy()
MM_snap_df = df_with2K_mmej_snap_back.apply(
    lambda x: main_markovian_process(
                sequance = x[f'ref_context_seq_{context_window}bp'] , 
                                 motif = x['snap_repeat_pat'],
                                motifs_d = x['snap_dist_between_reps'],
                    memory_dimer = x['snap_last_dimer'],
                    early_stop=1), axis = 1 ,result_type='expand')
MM_snap_df.rename(columns={0: 'SD_snap_back_prob',1: 'SD_snap_back_lk',
                2:'SD_snap_back_p_val', 3: 'SD_snap_back_dist_CDF', 4:'SD_snap_back_q'}, inplace=True)
df.loc[ MM_snap_df.index,
    ['SD_snap_back_prob','SD_snap_back_lk','SD_snap_back_p_val', 'SD_snap_back_dist_CDF']] = MM_snap_df.loc[
                            :, ['SD_snap_back_prob','SD_snap_back_lk','SD_snap_back_p_val', 'SD_snap_back_dist_CDF']]

df.loc[:,'snap_back_motif_count'] = np.nan
df.loc[MM_snap_df.index,'snap_back_motif_count'] = df.loc[MM_snap_df.index,: ].apply(lambda x: 
    len([m.start() for m in re.finditer(x['snap_repeat_pat'], x[f'ref_context_seq_{context_window}bp'], overlapped=True)]), axis = 1)
logging.debug('Snap back mmej Markov ends')

logging.debug('Loop out mmej Markov starts')
# creating a df with only mmej Loop out
df_with2K_mmej_loop_out = df[df['SD_loop_out'] == True].copy()
MM_loop_df = df_with2K_mmej_loop_out.apply(
    lambda x: main_markovian_process(sequance = x[f'ref_context_seq_{context_window}bp'] , 
                                 motif = x['loop_repeat_pat'],
                                motifs_d = x['loop_dist_between_reps'],
                    memory_dimer = x['loop_last_dimer'],
                    early_stop=1), axis = 1 ,result_type='expand')

MM_loop_df.rename(columns={0: 'SD_loop_out_prob',1: 'SD_loop_out_lk',
                2:'SD_loop_out_p_val', 3: 'SD_loop_out_dist_CDF', 4: 'SD_loop_out_q'}, inplace=True)
df.loc[ MM_loop_df.index,
    ['SD_loop_out_prob','SD_loop_out_lk', 'SD_loop_out_p_val','SD_loop_out_dist_CDF']] = MM_loop_df.loc[:, 
                            ['SD_loop_out_prob', 'SD_loop_out_lk', 'SD_loop_out_p_val', 'SD_loop_out_dist_CDF']]

df.loc[:,'loop_out_motif_count'] = np.nan
df.loc[MM_loop_df.index,'loop_out_motif_count'] = df.loc[MM_loop_df.index,: ].apply(lambda x: 
    len([m.start() for m in re.finditer(x['loop_repeat_pat'], x[f'ref_context_seq_{context_window}bp'], overlapped=True)]), axis = 1)

logging.debug('Loop out mmej Markov ends')

# logging.info(mem_reporter())

df['del_mmej_prob'] = 1 - df['del_mmej_prob']
# df['trans_mmej_prob'] = 1 - df['trans_mmej_prob']
df['SD_snap_back_prob'] = 1 - df['SD_snap_back_prob']
df['SD_loop_out_prob'] = 1 - df['SD_loop_out_prob']

df['del_mmej_prob'].fillna(value=0, inplace=True)
# df['trans_mmej_prob'].fillna(value=0, inplace=True)
df['SD_snap_back_prob'].fillna(value=0, inplace=True)
df['SD_loop_out_prob'].fillna(value=0, inplace=True)

df['del_mmej_lk'] = 1 - df['del_mmej_lk']
# df['trans_mmej_lk'] = 1 - df['trans_mmej_lk']
df['SD_snap_back_lk'] = 1 - df['SD_snap_back_lk']
df['SD_loop_out_lk'] = 1 - df['SD_loop_out_lk']

df['del_mmej_lk'].fillna(value=0, inplace=True)
# df['trans_mmej_lk'].fillna(value=0, inplace=True)
df['SD_snap_back_lk'].fillna(value=0, inplace=True)
df['SD_loop_out_lk'].fillna(value=0, inplace=True)

df['del_mmej_p_val'].fillna(value=1, inplace=True)
# df['trans_mmej_p_val'].fillna(value=1, inplace=True)
df['SD_snap_back_p_val'].fillna(value=1, inplace=True)
df['SD_loop_out_p_val'].fillna(value=1, inplace=True)

df['del_mmej_dist_CDF'].fillna(value=0, inplace=True)
# df['trans_dist_CDF'].fillna(value=0, inplace=True)
df['SD_snap_back_dist_CDF'].fillna(value=0, inplace=True)
df['SD_loop_out_dist_CDF'].fillna(value=0, inplace=True)



# logging.info(mem_reporter())

"""
Sice the Markov model calculates the probability to find 2 motifs in a
defined range and a given context, this probability is the probability
that something isn't MMEJ, therefore if one wants to have the probability
of something to be a MMEJ, he should do 1-probability
"""

df.loc[(df['indel_type'] == 'DEL'), 'indel_len'] = df.loc[(df['indel_type'] == 'DEL'), 'indel_len'] * (-1)

"""
Annotating NHEJ as 1-(all the other repair mechanisms)
Calculating NHEJ prob as 1-(max(all the other repair mechnisms prob))
"""

# insertions that are not any type of MMEJ will be called for now: unclassified_ins
df.loc[:,'unclassified_ins'] = False
df.loc[((df['indel_type'] == 'INS') &
    #    (df['SD_trans'] == False) &
       (df['SD_snap_back'] == False) &
       (df['SD_loop_out'] == False)),'unclassified_ins'] = True

df.loc[:,'NHEJ'] = False

df.loc[((df['indel_type'] == 'DEL') & (df['del_mmej'] == False)), 'NHEJ'] = True

df.loc[:,'NHEJ_prob'] = 0
df.loc[:,'NHEJ_lk'] = 0
df.loc[:,'NHEJ_p_val'] = 1
df.loc[:,'NHEJ_dist_CDF'] = 0
df.loc[:,'NHEJ_q'] = 0
df.loc[(df['indel_type'] == 'DEL'),'NHEJ_prob'] = df.apply(lambda x: (1-x['del_mmej_prob']), axis = 1)
df.loc[(df['indel_type'] == 'DEL'),'NHEJ_lk'] = df.apply(lambda x: (1-x['del_mmej_lk']), axis = 1)
df.loc[(df['indel_type'] == 'DEL'),'NHEJ_p_val'] = df.apply(lambda x: (1-x['del_mmej_p_val']), axis = 1)
df.loc[(df['indel_type'] == 'DEL'),'NHEJ_dist_CDF'] = df.apply(lambda x: (1-x['del_mmej_dist_CDF']), axis = 1)
df.loc[(df['indel_type'] == 'DEL'),'NHEJ_q'] = df.apply(lambda x: (1-x['del_mmej_q']), axis = 1)



df.loc[:,'unclassified_ins_prob'] = 0
df.loc[:,'unclassified_ins_lk'] = 0
df.loc[:,'unclassified_ins_p_val'] = 1
df.loc[:,'unclassified_ins_dist_CDF'] = 0
df.loc[(df['indel_type'] == 'INS'),'unclassified_ins_prob'] = df.apply(lambda x: (1-(x[[
                    #  'trans_mmej_prob',
                     'SD_snap_back_prob',
                     'SD_loop_out_prob']].max())), axis = 1)
df.loc[(df['indel_type'] == 'INS'),'unclassified_ins_p_val'] = df.apply(lambda x: (1-(x[[
                    #  'trans_mmej_p_val',
                     'SD_snap_back_p_val',
                     'SD_loop_out_p_val']].max())), axis = 1)

df.loc[(df['indel_type'] == 'INS'),'unclassified_ins_lk'] = df.apply(lambda x: (1-(x[[
                    #  'trans_mmej_lk',
                     'SD_snap_back_lk',
                     'SD_loop_out_lk']].max())), axis = 1)

df.loc[(df['indel_type'] == 'INS'),'unclassified_ins_dist_CDF'] = df.apply(lambda x: (1-(x[[
                    #  'trans_mmej_dist_CDF',
                     'SD_snap_back_dist_CDF',
                     'SD_loop_out_dist_CDF']].max())), axis = 1)


logging.debug('Old EM starts')
convergence_cutoff = 0.0000000045
df_del_mmej_lk = ExpectationMaximization(data=df,
                    mechanism_prob='del_mmej_lk',  
                    initial_theta=0.1,
                    convergence_threshold=convergence_cutoff)

df['del_mmej_post_decoded_lk'] = df_del_mmej_lk.posterior_decoding
logging.info(f'- df_del_mmej_lk theta = {df_del_mmej_lk.theta_a}')

df_del_mmej_q = ExpectationMaximization(data=df,
                    mechanism_prob='del_mmej_q',  
                    initial_theta=0.1,
                    convergence_threshold=convergence_cutoff)
df['del_mmej_post_decoded_q'] = df_del_mmej_q.posterior_decoding
logging.info(f'- df_del_mmej_q theta = {df_del_mmej_q.theta_a}')

df_del_mmej_prob = ExpectationMaximization(data=df,
                    mechanism_prob='del_mmej_prob',  
                    initial_theta=0.1,
                    convergence_threshold=convergence_cutoff)
df['del_mmej_post_decoded_prob'] = df_del_mmej_prob.posterior_decoding
logging.info(f'- df_del_mmej_prob theta = {df_del_mmej_prob.theta_a}')

df_del_mmej_p_val = ExpectationMaximization(data=df,
                    mechanism_prob='del_mmej_p_val',  
                    initial_theta=0.1,
                    convergence_threshold=convergence_cutoff)
df['del_mmej_post_decoded_p_val'] = df_del_mmej_p_val.posterior_decoding
logging.info(f'- df_del_mmej_p_val theta = {df_del_mmej_p_val.theta_a}')

df['NHEJ_post_decoded_lk'] = 1 - df['del_mmej_post_decoded_lk']
df['NHEJ_post_decoded_q'] = 1 - df['del_mmej_post_decoded_q']
df['NHEJ_post_decoded_prob'] = 1 - df['del_mmej_post_decoded_prob']
df['NHEJ_post_decoded_p_val'] = 1 - df['del_mmej_post_decoded_p_val']

logging.debug('Old EM ends')


# ------------------------- EMq --------------------------------------
def create_variant_id(chrom, original_pos):
    """Define name (id) for a set of potential realignments of the same original indel
    """
    return(chrom+"_"+str(original_pos).split(".")[0])

if df['original_pos'].isna().any():
    df['variant_id'] = df.apply(lambda x: create_variant_id(x['CHR'], (x['POS'])), axis=1)
else:
    df['variant_id'] = df.apply(lambda x: create_variant_id(x['CHR'], (x['original_pos'])), axis=1)

col_to_save = ['CHR', 'POS', 'original_pos', 'variant_id','REF', 
            'ALT', 'indel_type', 'indel_len', 
            f'ref_context_seq_{context_window}bp', 'accession_context']

df_for_EMq = df.loc[:,col_to_save]

EMq_obj = EMq(data=df,
                initial_theta=0.1,
                convergence_threshold=0.00000045,
                window_size=15)

print(f'EMq final theta MMEJ = {EMq_obj.theta_a}')

df['del_MMEJ_EMq_post_decoding'] = EMq_obj.del_MMEJ_post_decoding
df['del_NHEJ_EMq_post_decoding'] = EMq_obj.del_NHEJ_post_decoding

# ------------------------- EMq --------------------------------------



logging.info('############### OUTPUT DATA #################')
df.info(buf=buf2)
logging.info(buf2.getvalue())

col_to_save = ['CHR', 
        'POS', 'original_pos','variant_id', 'REF', 'ALT', 'del_mmej_prob','del_mmej_lk', 
        'del_mmej_p_val', 'del_mmej_dist_CDF', 'del_mmej_q','del_mmej_old_prob', 'del_mmej_motif_count',
        'del_mmej_post_decoded_prob', 'del_mmej_post_decoded_lk', 'del_mmej_post_decoded_p_val', 'del_mmej_post_decoded_q',
        'del_MMEJ_EMq_post_decoding', 
        'SD_snap_back_prob', 'SD_snap_back_lk','SD_snap_back_p_val', 'SD_snap_back_dist_CDF', 'snap_back_motif_count',
        'SD_loop_out_prob', 'SD_loop_out_lk', 'SD_loop_out_p_val', 'SD_loop_out_dist_CDF', 'loop_out_motif_count',
        'NHEJ_prob', 'NHEJ_lk','NHEJ_p_val', 'NHEJ_dist_CDF', 'NHEJ_q', 'NHEJ_old_prob',
        'NHEJ_post_decoded_lk', 'NHEJ_post_decoded_q', 'NHEJ_post_decoded_prob', 'NHEJ_post_decoded_p_val',
        'del_NHEJ_EMq_post_decoding', 
        'unclassified_ins_prob', 'unclassified_ins_lk', 'unclassified_ins_p_val', 'unclassified_ins_dist_CDF', 
        'del_mmej_cand_len', 'trans_mmej_cand_len', 'snap_mh2_len', 'loop_mh2_len']

df_mechanisms_prob = df.loc[:, col_to_save]
logging.debug({k:str(v[0]) for k,v in pd.DataFrame(df_mechanisms_prob.dtypes).T.to_dict('list').items()})

file_name = args["inputfilename"].split('.')[0]
df.to_csv(f'{args["outputfile"]}/{args["date"]}_{file_name}_classified_repair_mechanism.csv',
                                index=False, sep='\t')

df_mechanisms_prob.to_csv(f'{args["outputfile"]}/{args["date"]}_{file_name}_classified_repair_mechanism_prob.csv',
                                index=False, sep='\t')

print({k:str(v[0]) for k,v in pd.DataFrame(df_mechanisms_prob.dtypes).T.to_dict('list').items()})

logging.info('MMEJ detection pipline DONE')
print('MMEJ detection pipline DONE')