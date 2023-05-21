# 23.05.2022
"""
This script is a test for the function 
sd_snap_back_MMEJ (MicroHomology_module_v3)
to make sure it runs properly and return the same expected 
results as I develop it. 
"""
# importing the relevant libraries/modules:
import datetime
import timeit
from unicodedata import category, name
import pandas as pd
import numpy as np
import sys
print('Libreries imported')
sys.path.append('/home/labs/alevy/guyta/guy_master_project/src')
from MicroHomology_module_v3 import MicroHomology
from MMEJ_2nd_order_MM_v2 import motif_probabily_calc
from data_exploration_util import *
pd.options.display.max_colwidth = 1000
print('Modules imported')
mem_reporter()

Dtypes = {
    'CHR': str, 'POS': int, 'original_pos': float, 'REF': str,
    'ALT': str, 'indel_type': 'category', 'indel_len': int, 
    'ref_context_seq_1500bp': str, 'accession_context': str
}

path_to_data = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/20220517/1957_20220517_main_pipline_output_test3/1957_20220517_MMEJ_detection_formated_data.csv'

# df = pd.read_csv(path_to_data, sep='\t', dtype=Dtypes)
# df = df.iloc[6,:]
# # print(df)

# mmej_obj = MicroHomology(df['accession_context'], df['ALT'],
#                                         df['REF'],
#                                         df['ref_context_seq_1500bp'],
#                                         1499,
#                                         df['indel_len'],
#                                         df['indel_type'])
# print('mmej_obj assinged')
# mem_reporter()

###########################################################################################
########### THIS SECTION IS TO OPERATE WHEN DF HAS ONLY 1 ROW ########################

# df2 = mmej_obj.ex_data
# df = pd.DataFrame(df).T
# merged_df = pd.DataFrame(columns=(df.columns.to_list() + df2.columns.to_list()))
# merged_df.loc[0, df.columns] = df.loc[6,:]
# merged_df.loc[0, df2.columns] = df2.loc[0,:]
# print(merged_df.info())

###########################################################################################
########### THIS SECTION IS TO OPERATE WHEN DF HAS MORE THEN 1 ROW ########################
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
# indel_length -> the length of the indel
# indel_type -> the type of the indel (INS\DEL)
"""

df = pd.read_csv(path_to_data, sep='\t', dtype=Dtypes)
context_window = 1500
df = df.iloc[0:100, :]

microhomology_object = df.apply( lambda x: MicroHomology(x['accession_context'], x['ALT'],
                                        x['REF'],
                                        x['ref_context_seq_1500bp'],
                                        context_window,
                                        # x['indel_len'],
                                        x['indel_type']), axis = 1)

print('microhomology_object assinged')
mem_reporter()

df2 = None
for i in microhomology_object.index:
    if (i == 0):
        df2 = microhomology_object[i].ex_data
    else:
        df2 = pd.concat([df2, microhomology_object[i].ex_data])
        
df2.reset_index(drop = True, inplace = True)
df.reset_index(drop = True, inplace = True)
df = df.join(df2, how = 'right')
# df.drop(columns = ['index'], inplace = True)
# print(f'## MMEJ detection is done({args["inputfilename"]})')

"""
for insertions: the motif is the repeat ('tamplate_switch_repeat')
for  mmej snap back: the motif is the repeat ('S_tamplate_switch_repeat')
for deletions: the motif is the 'mmej_cand'
for all indels the last nucleotide is the one that comes right before the DSB on the reference seq
"""
# creating empty columns for the probabilities
df['del_mmej_prob'] = np.nan
df['trans_mmej_prob'] = np.nan
df['SD_snap_back_prob'] = np.nan
df['SD_loop_out_prob'] = np.nan

print('del_mmej: ',df['del_mmej'].value_counts(),'\n','###########\n', # 0
     'SD_trans: ',df['SD_trans'].value_counts(),'\n','###########\n', # 0
     'SD_snap_back: ', df['SD_snap_back'].value_counts(),'\n','###########\n', # 2
     'SD_loop_out: ', df['SD_loop_out'].value_counts(),'\n','###########\n') # 100


"""
The probability to find mmej will be calculate seperatly since each scenario has a different
motif. Other then the motif itself, the rest is similar.
"""
# creating a df with only deletions
df_with2K_del_mmej = df[((df['del_mmej'] == True) & 
           (df['indel_type'] == 'DEL'))].copy()

print('DEL mmej probability (Markov) starts', datetime.datetime.now())
df.loc[
    ((df['del_mmej'] == True) & (df['indel_type'] == 'DEL')),
    'del_mmej_prob'] = df_with2K_del_mmej.apply(lambda x: motif_probabily_calc(_sequance = x[f'ref_context_seq_{context_window}bp'] , 
                                 _motif = x['del_mmej_cand'],
                                _indel_len = int(x['indel_len']),
                    _memory_dimer = x['del_last_dimer']) , axis = 1 )
print('DEL mmej probability (Markov) ends', datetime.datetime.now())

# creating a df with only trans insertions
df_with2K_ins_mmej = df[(df['SD_trans'] == True)].copy()
df_with2K_ins_mmej.loc[:,'trans_mmej_cand_len'] = df_with2K_ins_mmej['trans_mmej_cand_len'].astype('int')
print('INS mmej probability (Markov) starts', datetime.datetime.now())
df.loc[
    (df['SD_trans'] == True),
    'trans_mmej_prob'] = df_with2K_ins_mmej.apply(lambda x: motif_probabily_calc(_sequance = x[f'ref_context_seq_{context_window}bp'] , 
                                 _motif = x['trans_TS_repeat'],
                                _indel_len = x['trans_mmej_cand_len'],
                    _memory_dimer = x['trans_last_dimer']), axis = 1 )
print('INS mmej probability (Markov) ends', datetime.datetime.now())


# creating a df with only mmej snap back
df_with2K_mmej_snap_back = df[df['SD_snap_back'] == True].copy()

print('snap back mmej probability (Markov) starts', datetime.datetime.now())
df.loc[
    (df['SD_snap_back'] == True),
    'SD_snap_back_prob'] = df_with2K_mmej_snap_back.apply(
    lambda x: motif_probabily_calc(_sequance = x[f'ref_context_seq_{context_window}bp'] , 
                                 _motif = x['snap_repeat_pat'],
                                _indel_len = x['snap_dist_between_reps'],
                    _memory_dimer = x['snap_last_dimer']), axis = 1 )

print('snap back mmej probability (Markov) ends', datetime.datetime.now())

# creating a df with only mmej Loop out
df_with2K_mmej_loop_out = df[df['SD_loop_out'] == True].copy()

print('Loop out mmej probability (Markov) starts', datetime.datetime.now())
df.loc[
    (df['SD_loop_out'] == True),
    'SD_loop_out_prob'] = df_with2K_mmej_loop_out.apply(
    lambda x: motif_probabily_calc(_sequance = x[f'ref_context_seq_{context_window}bp'] , 
                                 _motif = x['loop_repeat_pat'],
                                _indel_len = x['loop_dist_between_reps'],
                    _memory_dimer = x['loop_last_dimer']), axis = 1 )

print('Loop out mmej probability (Markov) ends', datetime.datetime.now())

mem_reporter()

"""
Sice the Markov model calculates the probability to find 2 motifs in a
defined range and a given context, this probability is the probability
that something isn't MMEJ, therefore if one wants to have the probability
of something to be a MMEJ, he should do 1-probability
"""
df['del_mmej_prob'] = 1 - df['del_mmej_prob']
df['del_mmej_prob'].fillna(value=0, inplace=True)
df['trans_mmej_prob'] = 1 - df['trans_mmej_prob']
df['trans_mmej_prob'].fillna(value=0, inplace=True)
df['SD_snap_back_prob'] = 1 - df['SD_snap_back_prob']
df['SD_snap_back_prob'].fillna(value=0, inplace=True)
df['SD_loop_out_prob'] = 1 - df['SD_loop_out_prob']
df['SD_loop_out_prob'].fillna(value=0, inplace=True)
pd.set_option("display.max_columns", None)

df.loc[(df['indel_type'] == 'DEL'), 'indel_len'] = df.loc[(df['indel_type'] == 'DEL'), 'indel_len'] * (-1)

"""
Annotating NHEJ as 1-(all the other repair mechanisms)
Calculating NHEJ prob as 1-(max(all the other repair mechnisms prob))
"""

# insertions that are not any type of MMEJ will be called for now: unclassified_ins
df.loc[:,'unclassified_ins'] = False
df.loc[((df['indel_type'] == 'INS') &
       (df['SD_trans'] == False) &
       (df['SD_snap_back'] == False) &
       (df['SD_loop_out'] == False)),'unclassified_ins'] = True

df.loc[:,'NHEJ'] = False

df.loc[((df['indel_type'] == 'DEL') & (df['del_mmej'] == False)), 'NHEJ'] = True

df.loc[:,'NHEJ_prob'] = 0
df.loc[(df['indel_type'] == 'DEL'),'NHEJ_prob'] = df.apply(lambda x: (1-x['del_mmej_prob']), axis = 1)


df.loc[:,'unclassified_ins_prob'] = 0
df.loc[(df['indel_type'] == 'INS'),'unclassified_ins_prob'] = df.apply(lambda x: (1-(x[[
                     'trans_mmej_prob',
                     'SD_snap_back_prob',
                     'SD_loop_out_prob']].max())), axis = 1)

print(df.info())
pd.set_option('display.max_rows', df.shape[0]+1)

###### WHEN UPDATING THE TEST PROCESS, UPDATE THE EXPECTED RESULTS FILE AS WELL ########
path_to_output='/home/labs/alevy/guyta/guy_master_project/data/repair_mechanisms_tests_expected_results/SD_loop_out'
filename='20220525_loop_out_test_correct_results.csv'
# df.to_csv(f'{path_to_output}/{filename}', sep='\t')


Dtypes = {'CHR' :object, 'POS' :int,'original_pos' :float,'REF':object,
        'ALT' :object,'indel_type':'category','indel_len':int,
        'ref_context_seq_1500bp':object,'accession_context':object,
        'del_mmej':object,'del_mmej_cand':object,'del_mmej_marked':object,
        'del_mmej_marked_on_ref':object,'del_last_dimer':object,'del_mmej_cand_len':object,
        'SD_trans':object,'trans_reps_pat':object,'trans_TS_repeat':object,
        'trans_last_dimer':object,'trans_mmej_cand_len':object,'trans_mmej_marked':object,
        'trans_mmej_marked_on_ref':object,'SD_snap_back':object,'snap_mmej_marked':object,
        'snap_P1':object,'snap_P2':object,'snap_mh1':object,'snap_mh2':object,'snap_repeat_pat':object,
        'snap_inv_comp_repeat_pat':object,'snap_last_dimer':object,'snap_dist_between_reps':object,
        'SD_loop_out':object,'loop_mmej_marked':object,'loop_P2':object,'loop_mh2':object,'loop_repeat_pat':object,
        'loop_last_dimer':object,'loop_dist_between_reps':object,'del_mmej_prob':float,'trans_mmej_prob':float,
        'SD_snap_back_prob':float,'SD_loop_out_prob':float,'unclassified_ins':object,'NHEJ':object,
        'NHEJ_prob':int,'unclassified_ins_prob':float}


#### TEST RESULTS REPOTING #### 
cols_to_compare = [#'SD_snap_back', 'SD_loop_out', 'del_mmej', 'SD_trans', 
                        #'NHEJ', 'unclassified_ins',
                    'del_mmej_prob', 'trans_mmej_prob',
                    'SD_snap_back_prob', 'SD_loop_out_prob',
                    'NHEJ_prob' ,
                    'unclassified_ins_prob', 'del_mmej_marked_on_ref',
                    'del_mmej_marked', 'trans_mmej_marked',
                    'trans_mmej_marked_on_ref', 'snap_mmej_marked',
                    'loop_mmej_marked']
results_df = pd.read_csv(f'{path_to_output}/{filename}', sep='\t', dtype=Dtypes)
results_df.drop(columns=['Unnamed: 0'], inplace=True)
print(results_df.info())

pd.testing.assert_frame_equal(df.loc[:, cols_to_compare], results_df.loc[:, cols_to_compare], 
        check_dtype=False, check_column_type=False)

print('#### TAST PASS ####')