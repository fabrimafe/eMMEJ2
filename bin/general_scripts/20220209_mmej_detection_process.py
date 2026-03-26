# 18.01.2022

"""
This notebook is basically how I put together all the different moduls
up to the point that we have detected the MMEJ and calculated the probability
that this is indeed an MMEJ using the 2nd order Markov model.
"""

# importing libreries
import datetime

import pandas as pd

import numpy as np
import sys
from time import localtime, strftime
import datetime
pd.options.display.max_colwidth = 2200
pd.set_option("display.max_columns", None)

sys.path.append('/home/labs/alevy/guyta/guy_master_project/src')
from VcfProcess_module import VcfProcess
from MicroHomology_module_v3 import MicroHomology
from MMEJ_2nd_order_MM_v2 import motif_probabily_calc

day = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
local_h = strftime("%H%M", localtime())
date = f'{local_h}_{day}'

print(date)
_chr = sys.argv[1]
_output_location = sys.argv[2] 

def mem_reporter():
    """
    Reports the RSS and VMS of a process
    No args.
    No return, just reporting VMS and RSS at a current time point
    """
    import os, psutil, datetime
    
    print(f'# Time: {datetime.datetime.now()}')
    print(f'# VMS(Mb):{psutil.Process(os.getpid()).memory_info().vms / 1024 ** 2}')
    print(f'# RSS(Mb):{psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2}')

pd.options.display.max_colwidth = 2400


print(f'Chr: {_chr},Date: {date}')

path_to_df = r'/home/labs/alevy/guyta/datasets/soybean/Liu_et.al.2020/analysis/Liu_data_exploration/20211006_Liu_annotated_chr_df.csv'
path_to_df_context = r'/home/labs/alevy/guyta/datasets/soybean/Liu_et.al.2020/analysis/Liu_getfasta/120bp_window/20211005_Liu_annotated_chr_df_getfasta_120window_intersected.fa'
path_to_df_context_2Kbp = r'/home/labs/alevy/guyta/datasets/soybean/Liu_et.al.2020/analysis/Liu_getfasta/20211216_2Kbp_window/20211216_Liu_annotated_chr_df_getfasta_2Kbp_window.fa'


df_columns = ['Chr', 'Start', 'end', 'Ref', 'Alt', 'ancesstral_allele',
        'derived_allele', 'derived_allele_freq', 'indel_type', 'indel_length',
        'indel_len_bin', 'Chr_num', 'Chr_intersection_format']
Dtypes = {
    'Chr' : 'str', 'Start' : 'int64', 'end': 'int64', 'Ref': 'str', 'Alt': 'str', 'ancesstral_allele': 'str',
        'derived_allele': 'str', 'derived_allele_freq': 'float64', 'indel_type': 'category', 'indel_length': 'int64',
        'indel_len_bin' : 'category', 'Chr_num' : 'int64', 'Chr_intersection_format': 'str'
}
df = VcfProcess(path_to_df, df_columns = df_columns, file_format = '.csv', df_dtypes = Dtypes)
df.data.drop(columns = ['Unnamed: 0'], inplace = True)

print('data loaded', df.data.shape, mem_reporter())

# creating a samll subset of the data
df.data = df.data[df.data['Chr'] == f'Chr{_chr}']
# df.data = df.data.iloc[:10000, :]

print('single chromosome data: ', df.data.shape, mem_reporter())

"""
generating the location of the fasta context in the format of the fasta intersected file
then using it to merge with the fasta context
"""
df.data.loc[:,'Start_minus120'] = df.data['Start'] - 120
df.data.loc[:,'end_plus120'] = df.data['end'] + 120
df.data.loc[:,'Start_minus1000'] = df.data['Start'] - 1000
df.data.loc[:,'end_plus1000'] = df.data['end'] + 1000

def Chr_intersection_fasta_format_genereter(Chr_intersection_format,
                                            Start_minus, end_plus ):
    return (f'>{Chr_intersection_format}:{Start_minus}-{end_plus}')

df.data.loc[:,'fasta_context_position'] = df.data.apply(lambda x:
                                Chr_intersection_fasta_format_genereter(x['Chr_intersection_format'],
                                                x['Start_minus120'], x['end_plus120']), axis=1)

df.data.loc[:,'fasta_context_position_2Kbp'] = df.data.apply(lambda x:
                                Chr_intersection_fasta_format_genereter(x['Chr_intersection_format'],
                                                x['Start_minus1000'], x['end_plus1000']), axis=1)

fasta_context_df = df.add_fasta_context_to_dataframe(path_to_df_context, 120)
df_context_2Kbp = df.add_fasta_context_to_dataframe(path_to_df_context_2Kbp, 1000)
df_context_2Kbp.rename(columns={"fasta_context_position": "fasta_context_position_2Kbp"}, inplace = True)

df = VcfProcess(path_to_df, df_columns = df_columns, file_format = '.csv', create_dataframe = False, df_dtypes = Dtypes,
    data = pd.merge(pd.DataFrame(df.data), pd.DataFrame(fasta_context_df), on="fasta_context_position"),
                    context_window = 120)

df.accession_context_generator('indel_type', 'derived_allele', 'ancesstral_allele')

# This section shows how I have checked wether the accession_context was calculeted properly or not:
df.data.loc[:,'accession_context_len'] = df.data['accession_context'].str.len()
df.data.loc[:,'fasta_context_seq_len'] = df.data['fasta_context_seq'].str.len()
df.data.loc[:,'REF_len'] = df.data['ancesstral_allele'].str.len()
df.data.loc[:,'ALT_len'] = df.data['derived_allele'].str.len()

df.data.loc[:,'insertion_check'] = (df.data['fasta_context_seq_len']
                                    + df.data['ALT_len']) ####
df.data.loc[:,'deletion_check'] = (df.data['fasta_context_seq_len']
                                    - df.data['REF_len'])

df.data.loc[:,'context_check'] = (((df.data['indel_type'] == 'INS')
                                    & (df.data['accession_context_len']
                                        == df.data['insertion_check']))
                                    | ((df.data['indel_type'] == 'DEL')
                                    & (df.data['accession_context_len']
                                        == df.data['deletion_check'])))

df.data.head()

"""
Looking at the problematic accesstions
"""
print('##problematic accession_context: \n',
    '##True\False counting: \n',
        df.data['context_check'].value_counts())
# print('\n##problematic accession_context: \n', df.data[df.data['context_check'] == False])


# looking for accetion in which the accession_context_generator() generated context that is'nt
# make sence in terms of length
problematic_accession_context = df.data['context_check'] != True
# # filtering out problematic_accession_context (3202 accesstions filtered here)
df.data = df.data[- problematic_accession_context]

# # adding the indel's relative position in the context
df.add_indel_pos_to_data('indel_type', 'ancesstral_allele')

df.data.drop(['accession_context_len', 'fasta_context_seq_len',
                            'insertion_check', 'deletion_check'], axis=1, inplace = True)


print(df.data.info())
print('mmej detection starts', mem_reporter())
"""
Detecting microhomology
"""
microhomology_object = df.data.apply( lambda x: MicroHomology(x['accession_context'],
                                                            x['derived_allele'],
                                        x['ancesstral_allele'],
                                        x['fasta_context_seq'],
                                        x['indel_pos'],
                                        x['indel_length'],
                                        x['indel_type']), axis = 1)

# extracting the data from the MicroHomology objects
print('extracting data from microhomology_object', datetime.datetime.now())

# memory allocating
df2 = pd.DataFrame(index=range(len(microhomology_object)),columns=microhomology_object[0].ex_data.columns)
# assigning the microhomology_objects to df2
for i in microhomology_object.index:
    df2.iloc[i,:] = microhomology_object[i].ex_data
        
print('extracting data from microhomology_object done', mem_reporter())
# print(df2)
df2.reset_index(drop = True, inplace = True)
# print(pd.DataFrame(df.data))
df = pd.DataFrame(df.data)
df = df.join(df2, how = 'right')
# df.drop(columns = ['index'], inplace = True)
# print(print(df.head()))
print('mmej detection ends', mem_reporter())

"""
Adding the long fasta context (2Kbp)
"""
df_with2K = pd.merge(pd.DataFrame(df), pd.DataFrame(df_context_2Kbp), on="fasta_context_position_2Kbp")

print(df_with2K.info())
# print(df_with2K.head())
df_with2K.rename(columns={"fasta_context_checker_y": "fasta_context_checker_2K",
                            'fasta_context_seq_y': 'reference_context_seq_2K',
                            "fasta_context_checker_x": "fasta_context_checker_120",
                            'fasta_context_seq_x': 'reference_context_seq_120'}, inplace = True)

# print(df_with2K.head())
print('mmej probability (Markov) starts', mem_reporter())

# for insertions: the motif is the repeat ('tamplate_switch_repeat')
# for  mmej snap back: the motif is the repeat ('S_tamplate_switch_repeat')
# for deletions: the motif is the 'mmej_cand'
# for all indels the last nucleotide is the one that comes right before the DSB on the reference seq

# creating empty columns for the probabilities
df_with2K.loc[:,'del_mmej_prob'] = np.nan
df_with2K.loc[:,'trans_mmej_prob'] = np.nan
df_with2K.loc[:,'SD_snap_back_prob'] = np.nan
df_with2K.loc[:,'SD_loop_out_prob'] = np.nan


df_with2K.loc[:,'indel_length'] = df_with2K['indel_length'].astype('int')

"""
The probability to find mmej will be calculate seperatly since each scenario has a different
motif. Other then the motif itself, the rest is similar.
"""
# creating a df with only deletions
df_with2K_del_mmej = df_with2K[((df_with2K['del_mmej'] == True) & 
           (df_with2K['indel_type'] == 'DEL'))].copy()

print('DEL mmej probability (Markov) starts', datetime.datetime.now())
    
df_with2K.loc[
    ((df_with2K['del_mmej'] == True) & (df_with2K['indel_type'] == 'DEL')),
    'del_mmej_prob'] = df_with2K_del_mmej.apply(lambda x: motif_probabily_calc(_sequance = x['reference_context_seq_2K'] , 
                                 _motif = x['del_mmej_cand'],
                                _indel_len = int(x['indel_length']),
                    _memory_dimer = x['del_last_dimer']) , axis = 1 )
print('DEL mmej probability (Markov) ends', datetime.datetime.now())
print(mem_reporter())

# creating a df with only trans insertions
df_with2K_ins_mmej = df_with2K[(df_with2K['SD_trans'] == True)].copy()

df_with2K_ins_mmej.loc[:,'trans_mmej_cand_len'] = df_with2K_ins_mmej['trans_mmej_cand_len'].astype('int')
print('INS mmej probability (Markov) starts', datetime.datetime.now())
df_with2K.loc[
    (df_with2K['SD_trans'] == True),
    'trans_mmej_prob'] = df_with2K_ins_mmej.apply(lambda x: motif_probabily_calc(_sequance = x['reference_context_seq_2K'] , 
                                 _motif = x['trans_TS_repeat'],
                                _indel_len = x['trans_mmej_cand_len'],
                    _memory_dimer = x['trans_last_dimer']), axis = 1 )
print('INS mmej probability (Markov) ends', datetime.datetime.now())
print(mem_reporter())

# creating a df with only mmej snap back
df_with2K_mmej_snap_back = df_with2K[df_with2K['SD_snap_back'] == True].copy()

print('snap back mmej probability (Markov) starts', datetime.datetime.now())
df_with2K.loc[
    (df_with2K['SD_snap_back'] == True),
    'SD_snap_back_prob'] = df_with2K_mmej_snap_back.apply(
    lambda x: motif_probabily_calc(_sequance = x['reference_context_seq_2K'] , 
                                 _motif = x['snap_repeat_pat'],
                                _indel_len = x['snap_dist_between_reps'],
                    _memory_dimer = x['snap_last_dimer']), axis = 1 )

print('snap back mmej probability (Markov) ends', datetime.datetime.now())
print(mem_reporter())
# creating a df with only mmej Loop out
df_with2K_mmej_loop_out = df_with2K[df_with2K['SD_loop_out'] == True].copy()

print('Loop out mmej probability (Markov) starts', datetime.datetime.now())
df_with2K.loc[
    (df_with2K['SD_loop_out'] == True),
    'SD_loop_out_prob'] = df_with2K_mmej_loop_out.apply(
    lambda x: motif_probabily_calc(_sequance = x['reference_context_seq_2K'] , 
                                 _motif = x['loop_repeat_pat'],
                                _indel_len = x['loop_dist_between_reps'],
                    _memory_dimer = x['loop_last_dimer']), axis = 1 )

print('Loop out mmej probability (Markov) ends', datetime.datetime.now())

# NHEJ: 
"""
Annotating NHEJ as 1-(all the other repair mechanisms)
Calculating NHEJ prob as 1-(max(all the other repair mechnisms prob))
"""

# insertions that are not any type of MMEJ will be called for now: NHEJ_ins
df_with2K.loc[:,'NHEJ_ins'] = False
df_with2K.loc[((df_with2K['indel_type'] == 'INS') &
       (df_with2K['SD_trans'] == False) &
       (df_with2K['SD_snap_back'] == False) &
       (df_with2K['SD_loop_out'] == False)),'NHEJ_ins'] = True

df_with2K.loc[:,'NHEJ'] = False
df_with2K.loc[((df_with2K['indel_type'] == 'DEL') & 
            (df_with2K['del_mmej'] == False)), 'NHEJ'] = True
# df_with2K.loc[:,'NHEJ'].value_counts()
df_with2K.loc[:,'NHEJ_prob'] = df_with2K.apply(lambda x: (1-(x[['del_mmej_prob',
                     'trans_mmej_prob',
                     'SD_snap_back_prob',
                     'SD_loop_out_prob']].max())), axis = 1)


df_with2K.loc[:,'NHEJ_ins_prob'] = 0
df_with2K.loc[(df_with2K['NHEJ_ins'] == True),'NHEJ_ins_prob'] = df_with2K.loc[(df_with2K['NHEJ_ins'] == True), 'NHEJ_prob']

df_with2K.loc[(df_with2K['NHEJ_ins'] == True), 'NHEJ_prob'] = 0

print(df_with2K.info())
# setting the columns that we want to keep
cols_to_save = ['Chr', 'Start', 'end', 'Ref', 'Alt', 'ancesstral_allele',
       'derived_allele', 'derived_allele_freq', 'indel_type', 'indel_pos', 'indel_length',
       'indel_len_bin','fasta_context_position', 'fasta_context_position_2Kbp',
       'accession_context', 'reference_context_seq_120', 'reference_context_seq_2K',
       'del_mmej',
       'del_mmej_cand', 'del_mmej_marked', 'del_mmej_marked_on_ref',
       'del_last_dimer', 'del_mmej_cand_len', 'SD_trans', 'trans_reps_pat',
       'trans_TS_repeat', 'trans_last_dimer', 'trans_mmej_cand_len',
       'trans_mmej_marked', 'SD_snap_back', 'snap_mmej_marked', 'snap_P1',
       'snap_P2', 'snap_mh1', 'snap_mh2', 'snap_repeat_pat',
       'snap_inv_comp_repeat_pat', 'snap_last_dimer', 'snap_dist_between_reps',
       'SD_loop_out', 'loop_mmej_marked', 'loop_P2', 'loop_mh2',
       'loop_repeat_pat', 'loop_last_dimer', 'loop_dist_between_reps',
       'del_mmej_prob', 'trans_mmej_prob', 'SD_snap_back_prob',
       'SD_loop_out_prob', 'NHEJ', 'NHEJ_prob', 'NHEJ_ins', 'NHEJ_ins_prob']


# creating the output dataframe
out_df = df_with2K.loc[:,cols_to_save]
# saving 
print('Saving starts', mem_reporter())
print('saving name and loc', f'{_output_location}/{date}_chr{_chr}_df_to_analyse.csv')
out_df.to_csv(f'{_output_location}/{date}_chr{_chr}_df_to_analyse.csv', sep = '\t')
print('Saving ends', mem_reporter())


