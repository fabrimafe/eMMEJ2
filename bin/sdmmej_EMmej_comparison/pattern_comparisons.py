# This script compare the actuall patterns that both algorithms find

import os

import pandas as pd
import numpy as np

pd.options.display.max_colwidth = 500

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3
def parallel(lst1, lst2):
    lst3 = [value for value in lst1 if value not in lst2]
    return lst3

res_pth = '/home/labs/alevy/guyta/guy_master_project/results/sdmmej_EMmej_comparison/plots/20221030'
WT_pth = f"{res_pth}/R0_WT_emmej_sdmmej_merged.tsv"
WT_ins_pth = f"{res_pth}/R0_WT_likeliest_ins.tsv"
POLQ_pth = f"{res_pth}/R0_POLQ_emmej_sdmmej_merged.tsv"
Lig4pth = f"{res_pth}/R0_Lig4_emmej_sdmmej_merged.tsv"

summary_table = pd.DataFrame(
    columns=['EMmej_only', 'sdmmej_only', 'both', 
    'MMEJ_matched_pattern',
    'Snap_matched_pattern',
    'Loop_matched_pattern', 'pattern_match_precentage'], index=['deletion_mmej', 'snap', 'loop'])

for pth in [WT_pth]: # , POLQ_pth, Lig4pth
    df = pd.read_csv(pth, sep='\t')
    # print(df.columns)
    df['sliced_variant_id'] = df.apply(lambda row: f"{str(row['variant_id'][:-2])}", axis=1)
    cols = ['CHR', 'POS_x', 'original_pos', 'variant_id', 'sliced_variant_id', 
        'direction','del_mmej', 'del_mmej_cand','MICROHOMOLOGY']
    EMmej_mmej = df.loc[df['del_mmej'], cols]
    sdmmej_mmej = df.loc[~df['MICROHOMOLOGY'].isna(), cols]

    # Deletions
    EMmej_only = parallel(EMmej_mmej['sliced_variant_id'].unique(), sdmmej_mmej['sliced_variant_id'].unique())
    sdmmej_only = parallel(sdmmej_mmej['sliced_variant_id'].unique(), EMmej_mmej['sliced_variant_id'].unique())
    EMmej_and_sdmmej = intersection(sdmmej_mmej['sliced_variant_id'].unique(), EMmej_mmej['sliced_variant_id'].unique())
    print(EMmej_only)
    summary_table.loc['deletion_mmej', 'EMmej_only'] = len(EMmej_only)
    summary_table.loc['deletion_mmej', 'sdmmej_only'] = len(sdmmej_only)
    summary_table.loc['deletion_mmej', 'both'] = len(EMmej_and_sdmmej)
    summary_table.loc['deletion_mmej', 'MMEJ_matched_pattern'] = EMmej_mmej.loc[
        EMmej_mmej['del_mmej_cand'] == EMmej_mmej['MICROHOMOLOGY'], 'sliced_variant_id'].unique().shape[0]
    summary_table.loc['deletion_mmej', 'pattern_match_precentage'] = summary_table.loc['deletion_mmej', 'MMEJ_matched_pattern']/summary_table.loc['deletion_mmej', 'both']
    
    
    # Insertions
    WT_ins_pth
    # df_ins = pd.read_csv(WT_ins_pth, sep='\t')
    df_ins = df.loc[df['indel_type'] == 'INS', :]
    cols = ['CHR', 'POS_x', 'original_pos', 'variant_id', 'sliced_variant_id', 'ancestral_indel', 'derived_indel',
        'direction','SD_snap_back', 'snap_mmej_marked','snap_repeat_pat','Snap-back',
        'SD_loop_out', 'loop_mmej_marked', 'loop_repeat_pat','Loop-out', 'RECONSTRUCTED_SEQ']
    
    df_ins['Loop-out'] = df_ins['Loop-out'].str.upper()
    df_ins['Snap-back'] = df_ins['Snap-back'].str.upper()

    # Snap-back
    df_ins.loc[((df_ins['direction'] == 1) & (df_ins['SD_snap_back'])), 'snap_repeat_pat'] = df_ins.loc[
        ((df_ins['direction'] == 1) & (df_ins['SD_snap_back'])), :].apply(lambda row: row['snap_repeat_pat'][::-1], axis=1)
        
    df_ins['snap_match'] = np.nan
    df_ins.loc[df_ins['SD_snap_back'], 'snap_match'] = df_ins.loc[df_ins['SD_snap_back'], :].apply(
        lambda row: str(row['snap_repeat_pat']) in str(row['Snap-back']), axis=1)

    EMmej_snap = df_ins.loc[(df_ins['SD_snap_back'] == True), 'sliced_variant_id'].unique()
    # print(EMmej_snap)
    # print('Iw7_CI_100_1082' in EMmej_snap)
    sdmmej_snap = df_ins.loc[((~df_ins['Snap-back'].isna()) & 
        (df_ins['Snap-back'] != "0") & 
        (df_ins['Snap-back'].str.contains('---')) &
        (df_ins['Snap-back'] != 0)), 'sliced_variant_id'].unique()
    EMmej_snap_only = parallel(EMmej_snap, sdmmej_snap)
    sdmmej_snap_only = parallel(sdmmej_snap, EMmej_snap)
    EMmej_and_sdmmej_snap = intersection(sdmmej_snap, EMmej_snap)
    
    summary_table.loc['snap', 'EMmej_only'] = len(EMmej_snap_only)
    summary_table.loc['snap', 'sdmmej_only'] = len(sdmmej_snap_only)
    summary_table.loc['snap', 'both'] = len(EMmej_and_sdmmej_snap)
    summary_table.loc['snap', 'Snap_matched_pattern'] = df_ins.loc[df_ins['snap_match'] == True, 'sliced_variant_id'].unique().shape[0]
    summary_table.loc['snap', 'pattern_match_precentage'] = summary_table.loc['snap', 'Snap_matched_pattern']/summary_table.loc['snap', 'both']

    # loop-out
    df_ins.loc[((df_ins['direction'] == 1) & (df_ins['SD_loop_out'])), 'loop_repeat_pat'] = df_ins.loc[
        ((df_ins['direction'] == 1) & (df_ins['SD_loop_out'])), :].apply(lambda row: row['loop_repeat_pat'][::-1], axis=1)
        
    df_ins['loop_match'] = np.nan
    df_ins.loc[df_ins['SD_loop_out'], 'loop_match'] = df_ins.loc[df_ins['SD_loop_out'], :].apply(
        lambda row: str(row['loop_repeat_pat']) in str(row['Loop-out']), axis=1)

    EMmej_loop = df_ins.loc[(df_ins['SD_loop_out'] == True), 'sliced_variant_id'].unique()
    # print(EMmej_loop)
    # print('Iw7_CI_100_1082' in EMmej_loop)
    sdmmej_loop = df_ins.loc[((~df_ins['Loop-out'].isna()) & (df_ins['Loop-out'] != '0')), 'sliced_variant_id'].unique()
    EMmej_loop_only = parallel(EMmej_loop, sdmmej_loop)
    sdmmej_loop_only = parallel(sdmmej_loop, EMmej_loop)
    EMmej_and_sdmmej_loop = intersection(sdmmej_loop, EMmej_loop)
    

    # print(df_ins.loc[df_ins['sliced_variant_id'].isin(sdmmej_loop_only), ['sliced_variant_id','loop_repeat_pat','Loop-out', 'RECONSTRUCTED_SEQ']])
    # print(df_ins.loc[df_ins['sliced_variant_id'].isin(sdmmej_loop_only), ['sliced_variant_id','loop_repeat_pat','Loop-out']])
    # print(df_ins.loc[df_ins['sliced_variant_id'].isin(sdmmej_loop_only), ['sliced_variant_id','Loop-out']])
    # print(df_ins.loc[df_ins['sliced_variant_id'].isin(sdmmej_loop_only), ['sliced_variant_id','loop_repeat_pat']])
    
    # print(df_ins.loc[df_ins['sliced_variant_id'] == 'Iw7_CI_100_1082', ['ancestral_indel', 'derived_indel','SD_loop_out']])
    # # print(df_ins.loc[df_ins['sliced_variant_id'] == 'Iw7_CI_100_1082', 'Loop-out'])
    # # print(df_ins.loc[df_ins['sliced_variant_id'] == 'Iw7_CI_100_1082', 'RECONSTRUCTED_SEQ'])
    # # print(df_ins.loc[df_ins['sliced_variant_id'] == 'Iw7_CI_53_1076', ['ancestral_indel', 'derived_indel',]])
    # # print(df_ins.loc[df_ins['sliced_variant_id'] == 'Iw7_CI_53_1076', 'Loop-out'])
    # # print(df_ins.loc[df_ins['sliced_variant_id'] == 'Iw7_CI_53_1076', 'RECONSTRUCTED_SEQ'])
    # print(df_ins.loc[df_ins['loop_match'] == True, ['loop_repeat_pat','Loop-out']])
    # print(df_ins.loc[df_ins['loop_match'] == True, 'sliced_variant_id'].unique().shape[0])
    summary_table.loc['loop', 'EMmej_only'] = len(EMmej_loop_only)
    summary_table.loc['loop', 'sdmmej_only'] = len(sdmmej_loop_only)
    summary_table.loc['loop', 'both'] = len(EMmej_and_sdmmej_loop)
    summary_table.loc['loop', 'Loop_matched_pattern'] = df_ins.loc[df_ins['loop_match'] == True, 'sliced_variant_id'].unique().shape[0]
    summary_table.loc['loop', 'pattern_match_precentage'] = summary_table.loc['loop', 'Loop_matched_pattern']/summary_table.loc['loop', 'both']

    print(summary_table)
    path_to_out = '/home/labs/alevy/guyta/guy_master_project/results/sdmmej_EMmej_comparison/plots/20221224'
    summary_table.to_csv(f"{path_to_out}/pattern_comparison_summary.tsv", sep='\t')