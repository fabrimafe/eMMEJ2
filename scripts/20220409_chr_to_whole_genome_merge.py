"""
This script takes the files from 20220207_mmej_detection_per_chr_submition.sh
and merge to one dataframe that contains all the chromosomes with the relevant
columns.

Example of how to operate with WEXAC
bsub -q new-short -R "rusage[mem=50000]" -e merge_err_try.txt -o merge_try.txt conda run -n guy_mmej_env python3 20220409_chr_to_whole_genome_merge.py '/home/labs/alevy/guyta/guy_master_project/results/soybean/Liu_et.al.2020/mmej_detection/1036_20220404'
"""

import sys
import os
import glob
sys.path.append('/home/labs/alevy/guyta/guy_master_project/src')
import datetime
from time import gmtime, strftime, localtime

import pandas as pd
import numpy as np

from data_exploration_util import *

# setting pandas display options
pd.options.display.max_colwidth = 2200
pd.set_option("display.max_columns", None)

mem_reporter()
date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
local_h = strftime("%H%M", localtime())

# geting file list
# wd = r'/home/labs/alevy/guyta/guy_master_project/results/soybean/Liu_et.al.2020/mmej_detection/1036_20220404'
wd = sys.argv[1]
output_folder = sys.argv[2]
print(f'### wd:\n {wd}')
print(f'### output folder:\n {output_folder}')

file_list = get_file_list(wd, 'csv')
print('file_list: ')
[print(file) for file in file_list]

# loading the annotations data
path_to_annotation_df = r'/home/labs/alevy/guyta/guy_master_project/results/soybean/Liu_et.al.2020/annotation_intersection/20220307_short_indels_grouped_annot_data.csv'
annotation_df_dtypes = {'Unnamed: 0' : str,'Chr': str,'start': int 
                ,'end': int ,'gene_ID' : str,'gene_annotation': 'category'
                ,'repeat_annotation': 'category','genic': bool  
                ,'transposon': bool,'Tandem_Repeat': bool  
                ,'LTR/Gypsy': bool ,'Unknown_rep_annot': bool  
                ,'LTR/Copia': bool ,'DNA/Mutator': bool  
                ,'NON_LTR/LINE': bool,'DNA/CACTA': bool  
                ,'DNA/stowaway': bool,'DNA/Helitron': bool  
                ,'DNA/tourist': bool ,'DNA/PIF': bool  
                ,'NON_LTR/SINE': bool,'LTR/Unclassified': bool  
                ,'DNA/PONG': bool,'tRNAScan-SE': bool  
                ,'DNA/Tc1': bool ,'DNA/hAT': bool  
                ,'nan_genic_annot': bool,'mRNA': bool  
                ,'gene': bool,'exon': bool  
                ,'CDS': bool ,'three_prime_UTR': bool  
                ,'five_prime_UTR': bool,'rRNA': bool  
                ,'miRNA_primary_transcript': bool,'snoRNA': bool  
                ,'tRNA': bool,'snRNA': bool,'miRNA' : bool,
                'max_gene_exp': float, 'cotyledon_exp_mean': float, 'leafbud_exp_mean': float,
           'leaf_exp_mean': float, 'stem_exp_mean': float, 'flo_exp_mean': float,
                    'pod_seed_exp_mean': float, 'pod_exp_mean': float, 'seed_exp_mean': float}

annotation_df = pd.read_csv(path_to_annotation_df, sep = '\t', dtype = annotation_df_dtypes)
# print(annotation_df.info())
#     mem_reporter()
#     annotation_df = annotation_df[annotation_df['Chr'] == 'Chr01']
annotation_df.drop(columns = ['Unnamed: 0'],inplace = True)
print(annotation_df.info())


columns_to_keep = ['Chr', 'Start','indel_length','indel_len_bin_pairs','indel_len_bin','del_mmej','SD_trans','SD_snap_back','SD_loop_out','NHEJ','del_mmej_prob', 'trans_mmej_prob', 'SD_snap_back_prob',
       'SD_loop_out_prob','NHEJ_prob','del_mmej_relative_prob', 'trans_mmej_relative_prob',
       'SD_snap_back_relative_prob', 'SD_loop_out_relative_prob',
       'NHEJ_relative_prob', 'repeat_annotation', 'genic', 'exon', 'loop_dist_between_reps', 'snap_dist_between_reps']

merged_chr_df = pd.DataFrame(columns=columns_to_keep)

dtypes = {'Unnamed: 0' : int, 'Chr' : str, 'Start' : int, 'end' : int ,
            'Ref' : str, 'Alt' : str, 'ancesstral_allele' : str, 'derived_allele' : str,
            'derived_allele_freq' : float, 'indel_type' : 'category', 'indel_pos' : int,
            'indel_length' : int, 'indel_len_bin' : 'category', 'fasta_context_position' : str,
            'fasta_context_position_2Kbp' : str, 'accession_context' : str,
            'reference_context_seq_120': str, 'reference_context_seq_2K': str, 
            'del_mmej_cand': str, 'del_mmej_marked': str, 'del_mmej_marked_on_ref': str,
            'del_last_dimer': str, 'del_mmej_cand_len': 'category', 'trans_reps_pat': str,
            'trans_TS_repeat': str, 'trans_last_dimer': str, 'trans_mmej_cand_len': 'category',
            'trans_mmej_marked': str,  'snap_mmej_marked': str, 'snap_P1': str,
            'snap_P2': str, 'snap_mh1': str, 'snap_mh2': str, 'snap_repeat_pat': str,
            'snap_inv_comp_repeat_pat': str, 'snap_last_dimer': str, 'snap_dist_between_reps': 'category',
            'loop_mmej_marked': str, 'loop_P2': str, 'loop_mh2': str,
            'loop_repeat_pat': str, 'loop_last_dimer': str, 'loop_dist_between_reps' : 'category',
            'del_mmej_prob': float, 'trans_mmej_prob': float, 'SD_snap_back_prob': float,
            'SD_loop_out_prob': float}

"""
Loading each file, process and concat to merged_chr_df
"""

for file in file_list:
#     print(merged_chr_df)
    print(f'##### Now working on:{file} #####')
    df = pd.read_csv(file, sep = '\t', dtype = dtypes)
    # df = df.iloc[0:100, :]
    mem_reporter()
    df.drop(columns = ['Alt', 'Ref',
                       'Unnamed: 0'], inplace = True)
    df.loc[:,['del_mmej', 'SD_trans', 'SD_snap_back', 'SD_loop_out']] = df.loc[:
                    ,['del_mmej', 'SD_trans', 'SD_snap_back', 'SD_loop_out']].fillna(value = False)

    df.loc[:,['del_mmej', 'SD_trans', 'SD_snap_back', 'SD_loop_out']].astype('bool')
    # df.info()
    mem_reporter()

    df.loc[(df['indel_type'] == 'DEL'), 'indel_length'] = df.loc[(df['indel_type'] == 'DEL'), 'indel_length'] * (-1)
    # Adding repeat and gene annotations
    df_without_annotations = df
    df_with_annotations = df.merge(annotation_df, on = ['Chr', 'Start', 'end'], how = 'left')
    df = df_with_annotations

    """
    Binning indel_length
    """
    bins = [i for i in range(-50,52,2)]
    bins = bins[0:round(len(bins)/2) -1] + [-1,1] + bins[round(len(bins)/2):]
    labels = [str(i) for i in bins]
    df.loc[:,'indel_len_bin_pairs'] = pd.cut(df.loc[:,'indel_length'], bins=bins, labels=labels[1:])
    df.loc[:,'indel_len_bin_pairs'] = df.loc[:,'indel_len_bin_pairs'].astype('float')

    conditions = [
        df['indel_len_bin_pairs'].isin(range(-10,11,1)),
        df['indel_len_bin_pairs'] > 10,
        df['indel_len_bin_pairs'] < (-10)
    ]

    choices = [
        df['indel_len_bin_pairs'],
        '10<indel_len',
        'indel_len<-10'
    ]
    df.loc[:,'indel_len_bin'] = np.select(conditions, choices)
    df.loc[:,'indel_len_bin'] = df.loc[:,'indel_len_bin'].astype('category')
    
    
    df.loc[ :,'SD_loop_out_prob'] = df.loc[ :,'SD_loop_out_prob'].fillna(value = 1)
    df.loc[ :,'SD_snap_back_prob'] = df.loc[ :,'SD_snap_back_prob'].fillna(value = 1)
    df.loc[ :,'del_mmej_prob'] = df.loc[ :,'del_mmej_prob'].fillna(value = 1)
    df.loc[ :,'trans_mmej_prob'] = df.loc[ :,'trans_mmej_prob'].fillna(value = 1)

    df.loc[ :,'SD_loop_out_prob'] = 1-df.loc[:,'SD_loop_out_prob']
    df.loc[ :,'SD_snap_back_prob'] = 1-df.loc[:,'SD_snap_back_prob']
    df.loc[ :,'del_mmej_prob'] = 1-df.loc[:,'del_mmej_prob']
    df.loc[ :,'trans_mmej_prob'] = 1-df.loc[:,'trans_mmej_prob']


    """
    Annotating NHEJ as 1-(all the other repair mechanisms)
    Calculating NHEJ prob as 1-(max(all the other repair mechnisms prob))
    """

    # insertions that are not any type of MMEJ will be called for now: NHEJ_ins
    df.loc[:,'NHEJ_ins'] = False
    df.loc[((df['indel_type'] == 'INS') &
           (df['SD_trans'] == False) &
           (df['SD_snap_back'] == False) &
           (df['SD_loop_out'] == False)),'NHEJ_ins'] = True

    df.loc[:,'NHEJ'] = False

    df.loc[((df['del_mmej'] == False) &
            (df['SD_trans'] == False) & 
            (df['SD_snap_back'] == False) & 
            (df['SD_loop_out'] == False) &
           (df['indel_length'] < 0)), 'NHEJ'] = True

    df.loc[:,'NHEJ_prob'] = df.apply(lambda x: (1-(x[['del_mmej_prob',
                         'trans_mmej_prob',
                         'SD_snap_back_prob',
                         'SD_loop_out_prob']].max())), axis = 1)


    df.loc[:,'NHEJ_ins_prob'] = 0
    df.loc[(df['NHEJ_ins'] == True),'NHEJ_ins_prob'] = df.loc[(df['NHEJ_ins'] == True), 'NHEJ_prob']
    df.loc[(df['NHEJ_ins'] == True), 'NHEJ_prob'] = 0
    # filtering out the insertion NHEJs:
    df = df.loc[~((df['NHEJ'] == True) & (df['indel_length'] > 0)), :].copy()
    df = df.loc[((df['indel_length'] != -50)), :].copy() # for some reason those are getting indel_len_bin = 0

    # Calculating the maximum relative probability between the mechanisms
    df['mechanisms_prob_sum'] = df.loc[:,['del_mmej_prob', 'trans_mmej_prob', 'SD_snap_back_prob',
           'SD_loop_out_prob', 'NHEJ_prob']].sum(axis=1)

    df['del_mmej_relative_prob'] = df['del_mmej_prob'] / df['mechanisms_prob_sum']
    df['trans_mmej_relative_prob'] = df['trans_mmej_prob'] / df['mechanisms_prob_sum']
    df['SD_snap_back_relative_prob'] = df['SD_snap_back_prob'] / df['mechanisms_prob_sum']
    df['SD_loop_out_relative_prob'] = df['SD_loop_out_prob'] / df['mechanisms_prob_sum']
    df['NHEJ_relative_prob'] = df['NHEJ_prob'] / df['mechanisms_prob_sum']

    no_repeat_df = df.loc[(df['repeat_annotation'] == '.'), :].copy()
    # no_repeat_df.loc[:,columns_to_keep].info()
    merged_chr_df = pd.concat([merged_chr_df, no_repeat_df.loc[:,columns_to_keep]], axis='rows')

merged_chr_df = merged_chr_df.sort_values(by='Chr').reset_index(drop=True)
print(merged_chr_df)
merged_chr_df.info()

print('START saving merged data')
merged_chr_df.to_csv(f'{output_folder}/{local_h}_{date}_merged_chr_df.csv', sep='\t')
print('DONE saving merged data')
