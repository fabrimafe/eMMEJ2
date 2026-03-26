"""
This script uses the simulations to generates
ROC curves based on which we decide what what
cutoff to use for classifying mehcanisms
"""



import os
import shutil
import glob
from subprocess import Popen, PIPE, STDOUT
import datetime
import logging
import regex as re
from operator import itemgetter

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.colors as mcolors
import matplotlib as mpl
import seaborn as sns

date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
organism = 'drosophila'
# out_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/20221102_sims/all_mechanisms_with_genomic_indel_length_dist/20221128_ROC'
# out_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/20221102_sims/all_mechanisms_with_genomic_indel_length_dist/20221129_ROC'
# out_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/20221102_sims/all_mechanisms_with_genomic_indel_length_dist/20221130_ROC'
if organism == 'drosophila':
    out_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/20221102_sims/all_mechanisms_with_genomic_indel_length_dist/20221204_ROC'
    out_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/20221102_sims/all_mechanisms_with_genomic_indel_length_dist/20230104_ROC_50bp_markov_window'
if organism == 'arabidopsis':
    out_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/arabidopsis_thaliana_sims/ROC_curves/20221204'
boot_n = 100
samp_size = 1000

markov_cols = ['CHR', 
        'POS', 'original_pos','variant_id', 'variant_id_N','direction','ancestral_indel', 'derived_indel', 
        'del_mmej_lk', 'del_mmej_p_val', 

        'SDMMEJ_lk','SD_snap_back_lk','SD_snap_back_lk_bc','SD_snap_back_p_val', 'pChance_snap',
        'snap_repeat_pat', 'snap_repeat_pat_len','snap_freq_motif_eq',
        'SD_loop_out_lk', 'SD_loop_out_lk_bc', 'SD_loop_out_p_val','pChance_loop', 'loop_freq_motif_eq',
        'loop_repeat_pat', 'loop_repeat_pat_len',
        'NHEJ_lk','NHEJ_p_val', 

        'unclassified_ins_lk', 'unclassified_ins_p_val',
        'indel_len','del_mmej_cand','del_mmej_cand_len','del_mmej_freq_motif_eq',
        'del_mmej_lk_MH2', 'del_mmej_p_val_MH2','pChance_MMEJ',
        'pChance_MMEJ_MH2', 'NHEJ_lk_MH2', 'NHEJ_p_val_MH2', 
        'del_mmej_motif_pos', 'del_mmej_freq_small_window', 'del_mmej_freq_large_window',
        'loop_motif_pos', 'loop_freq_small_window', 'loop_freq_large_window',
        'snap_motif_pos', 'snap_freq_small_window', 'snap_freq_large_window']

if organism == 'drosophila':
    sims_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/20221101_sims_from_fabrizio/sims_with_genomic_indel_lenght_dist/20221123_df_size_2000/20221125_EMmej_res/whole_sims_files'
    sims_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/20221101_sims_from_fabrizio/sims_with_genomic_indel_lenght_dist/20221221_df_size_2000/20230104_EMmej_res' # Drosophila melanogaster sims with 50bp marlov window
    # ins_1_and_4_bp_MH_SD_sims_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/20221101_sims_from_fabrizio/sims_with_genomic_indel_lenght_dist/20221204_df_size_2000/20221204_EMmej_res'
if organism == 'arabidopsis':
    sims_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/arabidopsis_taliana_sims_with_indel_len_dist/20221204_df_size_2000/20221204_EMmej_res'
    ins_1_and_4_bp_MH_SD_sims_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/arabidopsis_taliana_sims_with_indel_len_dist/20221204_df_size_2000/20221204_EMmej_res'

# fol = ['original_sims_EMmej_format', 'sims_in_repeats', 'sims_not_in_repeats']
# suf = ['_EMmej_format','_repeat', '_non_repeat']

fol = ['original_sims_EMmej_format']
suf = ['_EMmej_format']


shortest_indel = -2
MH_length = 2

for f,s in zip(fol, suf):
    print(f"=========={f}===========")
    if organism == 'drosophila':
        mmej_markov_df_pth = f'{sims_pth}/mmej_{f}_with_genomic_indel_len_dist_EMmej_output/mmej_{f}_with_genomic_indel_len_dist_markov_output.tsv'
        nhej_markov_df_pth = f'{sims_pth}/nhej_{f}_with_genomic_indel_len_dist_no_1bp_deletion_EMmej_output/nhej_{f}_with_genomic_indel_len_dist_no_1bp_deletion_markov_output.tsv'
        mmej_markov_df_pth_with_1bp_del = f'{sims_pth}/mmej_MH1_{f}_with_genomic_indel_len_dist_EMmej_output/mmej_MH1_{f}_with_genomic_indel_len_dist_markov_output.tsv'
        nhej_markov_df_pth_with_1bp_del = f'{sims_pth}/nhej_{f}_with_genomic_indel_len_dist_EMmej_output/nhej_{f}_with_genomic_indel_len_dist_markov_output.tsv'
        ins_markov_df_pth = f'{sims_pth}/ins_{f}_with_genomic_indel_len_dist_EMmej_output/ins_{f}_with_genomic_indel_len_dist_markov_output.tsv'
        snap_markov_df_pth = f'{sims_pth}/snap_{f}_with_genomic_indel_len_dist_EMmej_output/snap_{f}_with_genomic_indel_len_dist_markov_output.tsv'
        loop_markov_df_pth = f'{sims_pth}/loop_{f}_with_genomic_indel_len_dist_EMmej_output/loop_{f}_with_genomic_indel_len_dist_markov_output.tsv'
        # snap_markov_df_1bp_MH_SD_pth = f'{ins_1_and_4_bp_MH_SD_sims_pth}/snap_1bp_MH_SD_{f}_with_genomic_indel_len_dist_EMmej_output/snap_1bp_MH_SD_{f}_with_genomic_indel_len_dist_markov_output.tsv'
        # snap_markov_df_4bp_MH_SD_pth = f'{ins_1_and_4_bp_MH_SD_sims_pth}/snap_4bp_MH_SD_{f}_with_genomic_indel_len_dist_EMmej_output/snap_4bp_MH_SD_{f}_with_genomic_indel_len_dist_markov_output.tsv'
        # loop_markov_df_1bp_MH_SD_pth = f'{ins_1_and_4_bp_MH_SD_sims_pth}/loop_1bp_MH_SD_{f}_with_genomic_indel_len_dist_EMmej_output/loop_1bp_MH_SD_{f}_with_genomic_indel_len_dist_markov_output.tsv'
        # loop_markov_df_4bp_MH_SD_pth = f'{ins_1_and_4_bp_MH_SD_sims_pth}/loop_4bp_MH_SD_{f}_with_genomic_indel_len_dist_EMmej_output/loop_4bp_MH_SD_{f}_with_genomic_indel_len_dist_markov_output.tsv'
    
    if organism == 'arabidopsis':
        mmej_markov_df_pth = f'{sims_pth}/mmej_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/mmej_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
        nhej_markov_df_pth = f'{sims_pth}/nhej_{f}_with_genomic_indel_len_dist_no_1bp_deletion_EMmej_format_EMmej_output/nhej_{f}_with_genomic_indel_len_dist_no_1bp_deletion_EMmej_format_markov_output.tsv'
        mmej_markov_df_pth_with_1bp_del = f'{sims_pth}/mmej_MH1_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/mmej_MH1_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
        nhej_markov_df_pth_with_1bp_del = f'{sims_pth}/nhej_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/nhej_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
        ins_markov_df_pth = f'{sims_pth}/ins_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/ins_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
        snap_markov_df_pth = f'{sims_pth}/snap_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/snap_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
        loop_markov_df_pth = f'{sims_pth}/loop_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/loop_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'

        snap_markov_df_1bp_MH_SD_pth = f'{ins_1_and_4_bp_MH_SD_sims_pth}/snap_1bp_MH_SD_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/snap_1bp_MH_SD_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
        snap_markov_df_4bp_MH_SD_pth = f'{ins_1_and_4_bp_MH_SD_sims_pth}/snap_4bp_MH_SD_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/snap_4bp_MH_SD_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
        loop_markov_df_1bp_MH_SD_pth = f'{ins_1_and_4_bp_MH_SD_sims_pth}/loop_1bp_MH_SD_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/loop_1bp_MH_SD_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
        loop_markov_df_4bp_MH_SD_pth = f'{ins_1_and_4_bp_MH_SD_sims_pth}/loop_4bp_MH_SD_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/loop_4bp_MH_SD_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'

    # print(loop_markov_df_pth)
    outpth_parent = f"{out_pth}/"
    output = f'{date}_p_value'
    dir = os.path.join(f"{outpth_parent}/{output}")
    if not os.path.exists(dir):	
        if not os.path.exists(os.path.join(f"{outpth_parent}")):
            os.mkdir(os.path.join(f"{outpth_parent}"))
        os.mkdir(os.path.join(dir))

    boot_df = pd.DataFrame(columns=['bootN', 'ratio', 'cutoff',
                 'FPR', 'TPR']) # ,'FNR', 'TNR'
    
    mmej_markov_df = pd.read_csv(mmej_markov_df_pth, sep='\t')
    nhej_markov_df = pd.read_csv(nhej_markov_df_pth,sep='\t')
    mmej_markov_df_with_1bp = pd.read_csv(mmej_markov_df_pth_with_1bp_del, sep='\t')
    nhej_markov_df_with_1bp = pd.read_csv(nhej_markov_df_pth_with_1bp_del,sep='\t')
    snap_markov_df = pd.read_csv(snap_markov_df_pth,sep='\t')
    loop_markov_df = pd.read_csv(loop_markov_df_pth,sep='\t')
    ins_markov_df = pd.read_csv(ins_markov_df_pth,sep='\t')
    

    mmej_markov_df_no_1bp = mmej_markov_df.loc[mmej_markov_df['ancestral_indel'].str.len() > 2, :].copy()
    nhej_markov_df_no_1bp = nhej_markov_df.loc[nhej_markov_df['ancestral_indel'].str.len() > 2, :].copy()
    ins_markov_df_no_1bp = ins_markov_df.loc[ins_markov_df['derived_indel'].str.len() > 2, :].copy()
    loop_markov_df_no_1bp = loop_markov_df.loc[loop_markov_df['derived_indel'].str.len() > 2, :].copy()
    snap_markov_df_no_1bp = snap_markov_df.loc[snap_markov_df['derived_indel'].str.len() > 2, :].copy()

    # ins_markov_df_no_1bp = ins_markov_df.copy()
    # ins_markov_df_no_1bp.loc[ins_markov_df_no_1bp['loop_repeat_pat_len']<=5, 'SD_loop_out_p_val'] = 1
    # ins_markov_df_no_1bp.loc[ins_markov_df_no_1bp['snap_repeat_pat_len']<=5, 'SD_snap_back_p_val'] = 1
    
    # loop_markov_df_no_1bp = loop_markov_df.copy()
    # snap_markov_df_no_1bp = snap_markov_df.copy()
    # loop_markov_df_no_1bp.loc[loop_markov_df_no_1bp['loop_repeat_pat_len'] <= 5, 'SD_loop_out_p_val'] = 1
    # snap_markov_df_no_1bp.loc[snap_markov_df_no_1bp['snap_repeat_pat_len'] <= 5, 'SD_snap_back_p_val'] = 1

    # snap_markov_1MH_SD_df = pd.read_csv(snap_markov_df_1bp_MH_SD_pth,sep='\t')
    # snap_markov_4MH_SD_df = pd.read_csv(snap_markov_df_4bp_MH_SD_pth,sep='\t')
    # loop_markov_1MH_SD_df = pd.read_csv(loop_markov_df_1bp_MH_SD_pth,sep='\t')
    # loop_markov_4MH_SD_df = pd.read_csv(loop_markov_df_4bp_MH_SD_pth,sep='\t')
    # snap_markov_1MH_SD_df_no_1bp = snap_markov_1MH_SD_df.loc[snap_markov_1MH_SD_df['derived_indel'].str.len() > 2, :].copy()
    # snap_markov_4MH_SD_df_no_1bp = snap_markov_4MH_SD_df.loc[snap_markov_4MH_SD_df['derived_indel'].str.len() > 2, :].copy()
    # loop_markov_1MH_SD_df_no_1bp = loop_markov_1MH_SD_df.loc[loop_markov_1MH_SD_df['derived_indel'].str.len() > 2, :].copy()
    # loop_markov_4MH_SD_df_no_1bp = loop_markov_4MH_SD_df.loc[loop_markov_4MH_SD_df['derived_indel'].str.len() > 2, :].copy()
    # snap_markov_1MH_SD_df_no_1bp = snap_markov_1MH_SD_df.loc[snap_markov_1MH_SD_df['derived_indel'].str.len() > 2, :].copy()
    # snap_markov_4MH_SD_df_no_1bp = snap_markov_4MH_SD_df.loc[snap_markov_4MH_SD_df['derived_indel'].str.len() > 2, :].copy()
    # loop_markov_1MH_SD_df_no_1bp = loop_markov_1MH_SD_df.loc[loop_markov_1MH_SD_df['derived_indel'].str.len() > 2, :].copy()
    # loop_markov_4MH_SD_df_no_1bp = loop_markov_4MH_SD_df.loc[loop_markov_4MH_SD_df['derived_indel'].str.len() > 2, :].copy()
    
# MMEJ vs NHEJ
cuts = np.arange(0,0.005,0.0002).tolist() + np.arange(0.015,1.1,0.1).tolist()
cuts = cuts[:-1] + [1]
boot_df = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])

for b in range(boot_n):
    samp_mmej_var = pd.Series(mmej_markov_df_with_1bp['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_nhej_var = pd.Series(nhej_markov_df_with_1bp['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
        # for idx,cutoff in enumerate(np.arange(0,0.011,0.001)):
    for idx,cutoff in enumerate(cuts):
        mmej_markov_df_with_1bp['bonferroni_cutoff'] = cutoff / mmej_markov_df_with_1bp['variant_id_N']        
        nhej_markov_df_with_1bp['bonferroni_cutoff'] = cutoff / nhej_markov_df_with_1bp['variant_id_N']

        FP = nhej_markov_df_with_1bp.loc[
        (nhej_markov_df_with_1bp['variant_id'].isin(samp_nhej_var) & 
        (nhej_markov_df_with_1bp['del_mmej_p_val'] <= nhej_markov_df_with_1bp['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_nhej_var.shape[0] - FP

        TP = mmej_markov_df_with_1bp.loc[
            (mmej_markov_df_with_1bp['variant_id'].isin(samp_mmej_var) & 
            (mmej_markov_df_with_1bp['del_mmej_p_val'] <= mmej_markov_df_with_1bp['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_mmej_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)
        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df = pd.concat([boot_df, tmp_boot_df])
    boot_df.reset_index(inplace=True, drop=True)
    boot_df.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df = boot_df.groupby(by=['cutoff']).mean().reset_index()
boot_df.loc[(boot_df.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df.loc[(boot_df.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df['ins_len'] = 'all lengths'
chosen_cutoff = boot_df.loc[
    ((boot_df['FPR'] < 0.05) & 
    (boot_df['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()


boot_df_no_1bp = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])
for b in range(boot_n):
    samp_mmej_var = pd.Series(mmej_markov_df_no_1bp['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_nhej_var = pd.Series(nhej_markov_df_no_1bp['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
    for idx,cutoff in enumerate(cuts):
        mmej_markov_df_no_1bp['bonferroni_cutoff'] = cutoff / mmej_markov_df_no_1bp['variant_id_N']        
        nhej_markov_df_no_1bp['bonferroni_cutoff'] = cutoff / nhej_markov_df_no_1bp['variant_id_N']

        FP = nhej_markov_df_no_1bp.loc[
        (nhej_markov_df_no_1bp['variant_id'].isin(samp_nhej_var) & 
        (nhej_markov_df_no_1bp['del_mmej_p_val'] <= nhej_markov_df_no_1bp['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_nhej_var.shape[0] - FP

        TP = mmej_markov_df_no_1bp.loc[
            (mmej_markov_df_no_1bp['variant_id'].isin(samp_mmej_var) & 
            (mmej_markov_df_no_1bp['del_mmej_p_val'] <= mmej_markov_df_no_1bp['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_mmej_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)

        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df_no_1bp = pd.concat([boot_df_no_1bp, tmp_boot_df])
    boot_df_no_1bp.reset_index(inplace=True, drop=True)
    boot_df_no_1bp.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df_no_1bp = boot_df_no_1bp.groupby(by=['cutoff']).mean().reset_index()
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df_no_1bp['ins_len'] = 'Deletion length > 1'
chosen_cutoff_1bp = boot_df_no_1bp.loc[
    ((boot_df_no_1bp['FPR'] < 0.05) & 
    (boot_df_no_1bp['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()

boot_df = pd.concat([boot_df, boot_df_no_1bp])
boot_df.reset_index(inplace=True, drop=True)

ax = plt.figure(figsize = (7,5), dpi=80)
fig = ax.get_figure()
ax = sns.lineplot(data=boot_df, x='FPR', y='TPR', hue='ins_len')
plt.plot([0.05,0.05], [0,1.02], color='green', linestyle='dashed')
plt.plot([0, 1], [0, 1],ls="--", c=".3")

ax.legend(title='Deletion length')
ax.set(ylim=(0, 1.02), xlim=(0, 1), 
    title=f'MMEJ vs NHEJ')
fig.savefig(f"{dir}/MMEJ_NHEJ_ROC.svg", bbox_inches='tight')
boot_df.to_csv(f"{dir}/MMEJ_NHEJ.tsv", sep='\t', index=False)


# LOOP vs INS
cuts = np.arange(0,0.011,0.001).tolist() + np.arange(0.15,1.1,0.1).tolist()
cuts = cuts[:-1] + [1]
boot_df = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])

for b in range(boot_n):
    samp_loop_var = pd.Series(loop_markov_df['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_ins_var = pd.Series(ins_markov_df['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
        # for idx,cutoff in enumerate(np.arange(0,0.011,0.001)):
    for idx,cutoff in enumerate(cuts):
        loop_markov_df['bonferroni_cutoff'] = cutoff / loop_markov_df['variant_id_N']        
        ins_markov_df['bonferroni_cutoff'] = cutoff / ins_markov_df['variant_id_N']

        FP = ins_markov_df.loc[
        (ins_markov_df['variant_id'].isin(samp_ins_var) & 
        (ins_markov_df['SD_loop_out_p_val'] <= ins_markov_df['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_ins_var.shape[0] - FP

        TP = loop_markov_df.loc[
            (loop_markov_df['variant_id'].isin(samp_loop_var) & 
            (loop_markov_df['SD_loop_out_p_val'] <= loop_markov_df['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_loop_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)
        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df = pd.concat([boot_df, tmp_boot_df])
    boot_df.reset_index(inplace=True, drop=True)
    boot_df.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df = boot_df.groupby(by=['cutoff']).mean().reset_index()
boot_df.loc[(boot_df.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df.loc[(boot_df.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df['ins_len'] = 'all lengths'
chosen_cutoff = boot_df.loc[
    ((boot_df['FPR'] < 0.05) & 
    (boot_df['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()



boot_df_no_1bp = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])
for b in range(boot_n):
    samp_loop_var = pd.Series(loop_markov_df_no_1bp['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_ins_var = pd.Series(ins_markov_df_no_1bp['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
    for idx,cutoff in enumerate(cuts):
        loop_markov_df_no_1bp['bonferroni_cutoff'] = cutoff / loop_markov_df_no_1bp['variant_id_N']        
        ins_markov_df_no_1bp['bonferroni_cutoff'] = cutoff / ins_markov_df_no_1bp['variant_id_N']

        FP = ins_markov_df_no_1bp.loc[
        (ins_markov_df_no_1bp['variant_id'].isin(samp_ins_var) & 
        (ins_markov_df_no_1bp['SD_loop_out_p_val'] <= ins_markov_df_no_1bp['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_ins_var.shape[0] - FP

        TP = loop_markov_df_no_1bp.loc[
            (loop_markov_df_no_1bp['variant_id'].isin(samp_loop_var) & 
            (loop_markov_df_no_1bp['SD_loop_out_p_val'] <= loop_markov_df_no_1bp['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_loop_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)

        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df_no_1bp = pd.concat([boot_df_no_1bp, tmp_boot_df])
    boot_df_no_1bp.reset_index(inplace=True, drop=True)
    boot_df_no_1bp.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df_no_1bp = boot_df_no_1bp.groupby(by=['cutoff']).mean().reset_index()
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df_no_1bp['ins_len'] = 'insertions length > 1'
chosen_cutoff_1bp = boot_df_no_1bp.loc[
    ((boot_df_no_1bp['FPR'] < 0.05) & 
    (boot_df_no_1bp['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()

boot_df = pd.concat([boot_df, boot_df_no_1bp])
boot_df.reset_index(inplace=True, drop=True)

ax = plt.figure(figsize = (7,5), dpi=80)
fig = ax.get_figure()

ax = sns.lineplot(data=boot_df, x='FPR', y='TPR', hue='ins_len')
plt.plot([0.05,0.05], [0,1.02], color='green', linestyle='dashed')
plt.plot([0, 1], [0, 1],ls="--", c=".3")

plt.text(0.055,boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'TPR'],chosen_cutoff)
plt.plot(boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'FPR'],
        boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'TPR'], 
                 color='black', marker='o', linewidth=2, markersize=4,alpha=0.8)
plt.text(0.055,boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'TPR'],chosen_cutoff_1bp)
plt.plot(boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'FPR'],
        boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'TPR'], 
                 color='black', marker='o', linewidth=2, markersize=4,alpha=0.8)

ax.legend(title='Insertion length')
ax.set(ylim=(0, 1.02), xlim=(0, 1), 
    title=f'SD-Loop out vs NHEJ insertions')
fig.savefig(f"{dir}/LOOP_INS_ROC.svg", bbox_inches='tight')
boot_df.to_csv(f"{dir}/LOOP_INS.tsv", sep='\t', index=False)



# SNAP vs INS
cuts = np.arange(0,0.011,0.001).tolist() + np.arange(0.15,1.1,0.1).tolist()
cuts = cuts[:-1] + [1]
boot_df = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])

for b in range(boot_n):
    samp_snap_var = pd.Series(snap_markov_df['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_ins_var = pd.Series(ins_markov_df['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
        # for idx,cutoff in enumerate(np.arange(0,0.011,0.001)):
    for idx,cutoff in enumerate(cuts):
        snap_markov_df['bonferroni_cutoff'] = cutoff / snap_markov_df['variant_id_N']        
        ins_markov_df['bonferroni_cutoff'] = cutoff / ins_markov_df['variant_id_N']

        FP = ins_markov_df.loc[
        (ins_markov_df['variant_id'].isin(samp_ins_var) & 
        (ins_markov_df['SD_snap_back_p_val'] <= ins_markov_df['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_ins_var.shape[0] - FP

        TP = snap_markov_df.loc[
            (snap_markov_df['variant_id'].isin(samp_snap_var) & 
            (snap_markov_df['SD_snap_back_p_val'] <= snap_markov_df['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_snap_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)
        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df = pd.concat([boot_df, tmp_boot_df])
    boot_df.reset_index(inplace=True, drop=True)
    boot_df.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df = boot_df.groupby(by=['cutoff']).mean().reset_index()
boot_df.loc[(boot_df.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df.loc[(boot_df.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df['ins_len'] = 'all lengths'
chosen_cutoff = boot_df.loc[
    ((boot_df['FPR'] < 0.05) & 
    (boot_df['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()



boot_df_no_1bp = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])
for b in range(boot_n):
    samp_snap_var = pd.Series(snap_markov_df_no_1bp['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_ins_var = pd.Series(ins_markov_df_no_1bp['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
    for idx,cutoff in enumerate(cuts):
        snap_markov_df_no_1bp['bonferroni_cutoff'] = cutoff / snap_markov_df_no_1bp['variant_id_N']        
        ins_markov_df_no_1bp['bonferroni_cutoff'] = cutoff / ins_markov_df_no_1bp['variant_id_N']

        FP = ins_markov_df_no_1bp.loc[
        (ins_markov_df_no_1bp['variant_id'].isin(samp_ins_var) & 
        (ins_markov_df_no_1bp['SD_snap_back_p_val'] <= ins_markov_df_no_1bp['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_ins_var.shape[0] - FP

        TP = snap_markov_df_no_1bp.loc[
            (snap_markov_df_no_1bp['variant_id'].isin(samp_snap_var) & 
            (snap_markov_df_no_1bp['SD_snap_back_p_val'] <= snap_markov_df_no_1bp['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_snap_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)

        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df_no_1bp = pd.concat([boot_df_no_1bp, tmp_boot_df])
    boot_df_no_1bp.reset_index(inplace=True, drop=True)
    boot_df_no_1bp.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df_no_1bp = boot_df_no_1bp.groupby(by=['cutoff']).mean().reset_index()
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df_no_1bp['ins_len'] = 'insertions length > 1'
chosen_cutoff_1bp = boot_df_no_1bp.loc[
    ((boot_df_no_1bp['FPR'] < 0.05) & 
    (boot_df_no_1bp['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()

boot_df = pd.concat([boot_df, boot_df_no_1bp])
boot_df.reset_index(inplace=True, drop=True)

ax = plt.figure(figsize = (7,5), dpi=80)
fig = ax.get_figure()
ax = sns.lineplot(data=boot_df, x='FPR', y='TPR', hue='ins_len')
plt.plot([0.05,0.05], [0,1.02], color='green', linestyle='dashed')
plt.plot([0, 1], [0, 1],ls="--", c=".3")

plt.text(0.055,boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'TPR'],chosen_cutoff)
plt.plot(boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'FPR'],
        boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'TPR'], 
                 color='black', marker='o', linewidth=2, markersize=4,alpha=0.8)
plt.text(0.055,boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'TPR'],chosen_cutoff_1bp)
plt.plot(boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'FPR'],
        boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'TPR'], 
                 color='black', marker='o', linewidth=2, markersize=4,alpha=0.8)


ax.legend(title='Insertion length')
ax.set(ylim=(0, 1.02), xlim=(0, 1), 
    title=f'SD Snap back vs NHEJ insertions')
fig.savefig(f"{dir}/SNAP_INS_ROC.svg", bbox_inches='tight')
boot_df.to_csv(f"{dir}/SNAP_INS.tsv", sep='\t', index=False)





# ------ MH4 SD4 ----------

# LOOP vs INS
cuts = np.arange(0,0.011,0.001).tolist() + np.arange(0.15,1.1,0.1).tolist()
cuts = cuts[:-1] + [1]
boot_df = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])

for b in range(boot_n):
    samp_loop_var = pd.Series(loop_markov_4MH_SD_df['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_ins_var = pd.Series(ins_markov_df['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
        # for idx,cutoff in enumerate(np.arange(0,0.011,0.001)):
    for idx,cutoff in enumerate(cuts):
        loop_markov_4MH_SD_df['bonferroni_cutoff'] = cutoff / loop_markov_4MH_SD_df['variant_id_N']        
        ins_markov_df['bonferroni_cutoff'] = cutoff / ins_markov_df['variant_id_N']

        FP = ins_markov_df.loc[
        (ins_markov_df['variant_id'].isin(samp_ins_var) & 
        (ins_markov_df['SD_loop_out_p_val'] <= ins_markov_df['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_ins_var.shape[0] - FP

        TP = loop_markov_4MH_SD_df.loc[
            (loop_markov_4MH_SD_df['variant_id'].isin(samp_loop_var) & 
            (loop_markov_4MH_SD_df['SD_loop_out_p_val'] <= loop_markov_4MH_SD_df['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_loop_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)
        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df = pd.concat([boot_df, tmp_boot_df])
    boot_df.reset_index(inplace=True, drop=True)
    boot_df.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df = boot_df.groupby(by=['cutoff']).mean().reset_index()
boot_df.loc[(boot_df.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df.loc[(boot_df.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df['ins_len'] = 'all lengths'
chosen_cutoff = boot_df.loc[
    ((boot_df['FPR'] < 0.05) & 
    (boot_df['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()



boot_df_no_1bp = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])
for b in range(boot_n):
    samp_loop_var = pd.Series(loop_markov_4MH_SD_df_no_1bp['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_ins_var = pd.Series(ins_markov_df_no_1bp['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
    for idx,cutoff in enumerate(cuts):
        loop_markov_4MH_SD_df_no_1bp['bonferroni_cutoff'] = cutoff / loop_markov_4MH_SD_df_no_1bp['variant_id_N']        
        ins_markov_df_no_1bp['bonferroni_cutoff'] = cutoff / ins_markov_df_no_1bp['variant_id_N']

        FP = ins_markov_df_no_1bp.loc[
        (ins_markov_df_no_1bp['variant_id'].isin(samp_ins_var) & 
        (ins_markov_df_no_1bp['SD_loop_out_p_val'] <= ins_markov_df_no_1bp['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_ins_var.shape[0] - FP

        TP = loop_markov_4MH_SD_df_no_1bp.loc[
            (loop_markov_4MH_SD_df_no_1bp['variant_id'].isin(samp_loop_var) & 
            (loop_markov_4MH_SD_df_no_1bp['SD_loop_out_p_val'] <= loop_markov_4MH_SD_df_no_1bp['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_loop_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)

        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df_no_1bp = pd.concat([boot_df_no_1bp, tmp_boot_df])
    boot_df_no_1bp.reset_index(inplace=True, drop=True)
    boot_df_no_1bp.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df_no_1bp = boot_df_no_1bp.groupby(by=['cutoff']).mean().reset_index()
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df_no_1bp['ins_len'] = 'insertions length > 1'
chosen_cutoff_1bp = boot_df_no_1bp.loc[
    ((boot_df_no_1bp['FPR'] < 0.05) & 
    (boot_df_no_1bp['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()

boot_df = pd.concat([boot_df, boot_df_no_1bp])
boot_df.reset_index(inplace=True, drop=True)

ax = plt.figure(figsize = (7,5), dpi=80)
fig = ax.get_figure()

ax = sns.lineplot(data=boot_df, x='FPR', y='TPR', hue='ins_len')
plt.plot([0.05,0.05], [0,1.02], color='green', linestyle='dashed')
plt.plot([0, 1], [0, 1],ls="--", c=".3")

plt.text(0.055,boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'TPR'],chosen_cutoff)
plt.plot(boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'FPR'],
        boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'TPR'], 
                 color='black', marker='o', linewidth=2, markersize=4,alpha=0.8)
plt.text(0.055,boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'TPR'],chosen_cutoff_1bp)
plt.plot(boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'FPR'],
        boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'TPR'], 
                 color='black', marker='o', linewidth=2, markersize=4,alpha=0.8)

ax.legend(title='Insertion length')
ax.set(ylim=(0, 1.02), xlim=(0, 1), 
    title=f'SD-Loop out vs NHEJ insertions')
fig.savefig(f"{dir}/LOOP_MH4_SD4_INS_ROC.svg", bbox_inches='tight')
boot_df.to_csv(f"{dir}/LOOP_MH4_SD4_INS.tsv", sep='\t', index=False)



# SNAP vs INS
cuts = np.arange(0,0.011,0.001).tolist() + np.arange(0.15,1.1,0.1).tolist()
cuts = cuts[:-1] + [1]
boot_df = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])

for b in range(boot_n):
    samp_snap_var = pd.Series(snap_markov_4MH_SD_df['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_ins_var = pd.Series(ins_markov_df['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
        # for idx,cutoff in enumerate(np.arange(0,0.011,0.001)):
    for idx,cutoff in enumerate(cuts):
        snap_markov_4MH_SD_df['bonferroni_cutoff'] = cutoff / snap_markov_4MH_SD_df['variant_id_N']        
        ins_markov_df['bonferroni_cutoff'] = cutoff / ins_markov_df['variant_id_N']

        FP = ins_markov_df.loc[
        (ins_markov_df['variant_id'].isin(samp_ins_var) & 
        (ins_markov_df['SD_snap_back_p_val'] <= ins_markov_df['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_ins_var.shape[0] - FP

        TP = snap_markov_4MH_SD_df.loc[
            (snap_markov_4MH_SD_df['variant_id'].isin(samp_snap_var) & 
            (snap_markov_4MH_SD_df['SD_snap_back_p_val'] <= snap_markov_4MH_SD_df['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_snap_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)
        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df = pd.concat([boot_df, tmp_boot_df])
    boot_df.reset_index(inplace=True, drop=True)
    boot_df.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df = boot_df.groupby(by=['cutoff']).mean().reset_index()
boot_df.loc[(boot_df.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df.loc[(boot_df.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df['ins_len'] = 'all lengths'
chosen_cutoff = boot_df.loc[
    ((boot_df['FPR'] < 0.05) & 
    (boot_df['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()



boot_df_no_1bp = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])
for b in range(boot_n):
    samp_snap_var = pd.Series(snap_markov_4MH_SD_df_no_1bp['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_ins_var = pd.Series(ins_markov_df_no_1bp['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
    for idx,cutoff in enumerate(cuts):
        snap_markov_4MH_SD_df_no_1bp['bonferroni_cutoff'] = cutoff / snap_markov_4MH_SD_df_no_1bp['variant_id_N']        
        ins_markov_df_no_1bp['bonferroni_cutoff'] = cutoff / ins_markov_df_no_1bp['variant_id_N']

        FP = ins_markov_df_no_1bp.loc[
        (ins_markov_df_no_1bp['variant_id'].isin(samp_ins_var) & 
        (ins_markov_df_no_1bp['SD_snap_back_p_val'] <= ins_markov_df_no_1bp['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_ins_var.shape[0] - FP

        TP = snap_markov_4MH_SD_df_no_1bp.loc[
            (snap_markov_4MH_SD_df_no_1bp['variant_id'].isin(samp_snap_var) & 
            (snap_markov_4MH_SD_df_no_1bp['SD_snap_back_p_val'] <= snap_markov_4MH_SD_df_no_1bp['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_snap_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)

        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df_no_1bp = pd.concat([boot_df_no_1bp, tmp_boot_df])
    boot_df_no_1bp.reset_index(inplace=True, drop=True)
    boot_df_no_1bp.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df_no_1bp = boot_df_no_1bp.groupby(by=['cutoff']).mean().reset_index()
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df_no_1bp['ins_len'] = 'insertions length > 1'
chosen_cutoff_1bp = boot_df_no_1bp.loc[
    ((boot_df_no_1bp['FPR'] < 0.05) & 
    (boot_df_no_1bp['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()

boot_df = pd.concat([boot_df, boot_df_no_1bp])
boot_df.reset_index(inplace=True, drop=True)

ax = plt.figure(figsize = (7,5), dpi=80)
fig = ax.get_figure()
ax = sns.lineplot(data=boot_df, x='FPR', y='TPR', hue='ins_len')
plt.plot([0.05,0.05], [0,1.02], color='green', linestyle='dashed')
plt.plot([0, 1], [0, 1],ls="--", c=".3")

plt.text(0.055,boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'TPR'],chosen_cutoff)
plt.plot(boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'FPR'],
        boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'TPR'], 
                 color='black', marker='o', linewidth=2, markersize=4,alpha=0.8)
plt.text(0.055,boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'TPR'],chosen_cutoff_1bp)
plt.plot(boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'FPR'],
        boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'TPR'], 
                 color='black', marker='o', linewidth=2, markersize=4,alpha=0.8)


ax.legend(title='Insertion length')
ax.set(ylim=(0, 1.02), xlim=(0, 1), 
    title=f'SD Snap back vs NHEJ insertions')
fig.savefig(f"{dir}/SNAP_MH4_SD4_INS_ROC.svg", bbox_inches='tight')
boot_df.to_csv(f"{dir}/SNAP_MH4_SD4_INS.tsv", sep='\t', index=False)





# ------ MH1 SD1 ----------

# LOOP vs INS
cuts = np.arange(0,0.011,0.001).tolist() + np.arange(0.15,1.1,0.1).tolist()
cuts = cuts[:-1] + [1]
boot_df = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])

for b in range(boot_n):
    samp_loop_var = pd.Series(loop_markov_1MH_SD_df['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_ins_var = pd.Series(ins_markov_df['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
        # for idx,cutoff in enumerate(np.arange(0,0.011,0.001)):
    for idx,cutoff in enumerate(cuts):
        loop_markov_1MH_SD_df['bonferroni_cutoff'] = cutoff / loop_markov_1MH_SD_df['variant_id_N']        
        ins_markov_df['bonferroni_cutoff'] = cutoff / ins_markov_df['variant_id_N']

        FP = ins_markov_df.loc[
        (ins_markov_df['variant_id'].isin(samp_ins_var) & 
        (ins_markov_df['SD_loop_out_p_val'] <= ins_markov_df['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_ins_var.shape[0] - FP

        TP = loop_markov_1MH_SD_df.loc[
            (loop_markov_1MH_SD_df['variant_id'].isin(samp_loop_var) & 
            (loop_markov_1MH_SD_df['SD_loop_out_p_val'] <= loop_markov_1MH_SD_df['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_loop_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)
        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df = pd.concat([boot_df, tmp_boot_df])
    boot_df.reset_index(inplace=True, drop=True)
    boot_df.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df = boot_df.groupby(by=['cutoff']).mean().reset_index()
boot_df.loc[(boot_df.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df.loc[(boot_df.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df['ins_len'] = 'all lengths'
chosen_cutoff = boot_df.loc[
    ((boot_df['FPR'] < 0.05) & 
    (boot_df['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()



boot_df_no_1bp = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])
for b in range(boot_n):
    samp_loop_var = pd.Series(loop_markov_1MH_SD_df_no_1bp['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_ins_var = pd.Series(ins_markov_df_no_1bp['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
    for idx,cutoff in enumerate(cuts):
        loop_markov_1MH_SD_df_no_1bp['bonferroni_cutoff'] = cutoff / loop_markov_1MH_SD_df_no_1bp['variant_id_N']        
        ins_markov_df_no_1bp['bonferroni_cutoff'] = cutoff / ins_markov_df_no_1bp['variant_id_N']

        FP = ins_markov_df_no_1bp.loc[
        (ins_markov_df_no_1bp['variant_id'].isin(samp_ins_var) & 
        (ins_markov_df_no_1bp['SD_loop_out_p_val'] <= ins_markov_df_no_1bp['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_ins_var.shape[0] - FP

        TP = loop_markov_1MH_SD_df_no_1bp.loc[
            (loop_markov_1MH_SD_df_no_1bp['variant_id'].isin(samp_loop_var) & 
            (loop_markov_1MH_SD_df_no_1bp['SD_loop_out_p_val'] <= loop_markov_1MH_SD_df_no_1bp['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_loop_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)

        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df_no_1bp = pd.concat([boot_df_no_1bp, tmp_boot_df])
    boot_df_no_1bp.reset_index(inplace=True, drop=True)
    boot_df_no_1bp.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df_no_1bp = boot_df_no_1bp.groupby(by=['cutoff']).mean().reset_index()
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df_no_1bp['ins_len'] = 'insertions length > 1'
chosen_cutoff_1bp = boot_df_no_1bp.loc[
    ((boot_df_no_1bp['FPR'] < 0.05) & 
    (boot_df_no_1bp['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()

boot_df = pd.concat([boot_df, boot_df_no_1bp])
boot_df.reset_index(inplace=True, drop=True)

ax = plt.figure(figsize = (7,5), dpi=80)
fig = ax.get_figure()

ax = sns.lineplot(data=boot_df, x='FPR', y='TPR', hue='ins_len')
plt.plot([0.05,0.05], [0,1.02], color='green', linestyle='dashed')
plt.plot([0, 1], [0, 1],ls="--", c=".3")

plt.text(0.055,boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'TPR'],chosen_cutoff)
plt.plot(boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'FPR'],
        boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'TPR'], 
                 color='black', marker='o', linewidth=2, markersize=4,alpha=0.8)
plt.text(0.055,boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'TPR'],chosen_cutoff_1bp)
plt.plot(boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'FPR'],
        boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'TPR'], 
                 color='black', marker='o', linewidth=2, markersize=4,alpha=0.8)

ax.legend(title='Insertion length')
ax.set(ylim=(0, 1.02), xlim=(0, 1), 
    title=f'SD-Loop out vs NHEJ insertions')
fig.savefig(f"{dir}/LOOP_MH1_SD1_INS_ROC.svg", bbox_inches='tight')
boot_df.to_csv(f"{dir}/LOOP_MH1_SD1_INS.tsv", sep='\t', index=False)



# SNAP vs INS
cuts = np.arange(0,0.011,0.001).tolist() + np.arange(0.15,1.1,0.1).tolist()
cuts = cuts[:-1] + [1]
boot_df = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])

for b in range(boot_n):
    samp_snap_var = pd.Series(snap_markov_1MH_SD_df['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_ins_var = pd.Series(ins_markov_df['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
        # for idx,cutoff in enumerate(np.arange(0,0.011,0.001)):
    for idx,cutoff in enumerate(cuts):
        snap_markov_1MH_SD_df['bonferroni_cutoff'] = cutoff / snap_markov_1MH_SD_df['variant_id_N']        
        ins_markov_df['bonferroni_cutoff'] = cutoff / ins_markov_df['variant_id_N']

        FP = ins_markov_df.loc[
        (ins_markov_df['variant_id'].isin(samp_ins_var) & 
        (ins_markov_df['SD_snap_back_p_val'] <= ins_markov_df['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_ins_var.shape[0] - FP

        TP = snap_markov_1MH_SD_df.loc[
            (snap_markov_1MH_SD_df['variant_id'].isin(samp_snap_var) & 
            (snap_markov_1MH_SD_df['SD_snap_back_p_val'] <= snap_markov_1MH_SD_df['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_snap_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)
        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df = pd.concat([boot_df, tmp_boot_df])
    boot_df.reset_index(inplace=True, drop=True)
    boot_df.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df = boot_df.groupby(by=['cutoff']).mean().reset_index()
boot_df.loc[(boot_df.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df.loc[(boot_df.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df['ins_len'] = 'all lengths'
chosen_cutoff = boot_df.loc[
    ((boot_df['FPR'] < 0.05) & 
    (boot_df['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()



boot_df_no_1bp = pd.DataFrame(columns=['bootN', 'cutoff',
                    'FPR', 'TPR'])
for b in range(boot_n):
    samp_snap_var = pd.Series(snap_markov_1MH_SD_df_no_1bp['variant_id'].unique()).sample(
        int(round(samp_size*0.5, 0)), replace=False)
    samp_ins_var = pd.Series(ins_markov_df_no_1bp['variant_id'].unique()).sample(
                int(round(samp_size*(1-0.5), 0)), replace=False)
    tmp_boot_df = pd.DataFrame(columns=
        ['bootN', 'cutoff', 'FPR', 'TPR'], 
                        index=[i for i in range(len(cuts))])
    for idx,cutoff in enumerate(cuts):
        snap_markov_1MH_SD_df_no_1bp['bonferroni_cutoff'] = cutoff / snap_markov_1MH_SD_df_no_1bp['variant_id_N']        
        ins_markov_df_no_1bp['bonferroni_cutoff'] = cutoff / ins_markov_df_no_1bp['variant_id_N']

        FP = ins_markov_df_no_1bp.loc[
        (ins_markov_df_no_1bp['variant_id'].isin(samp_ins_var) & 
        (ins_markov_df_no_1bp['SD_snap_back_p_val'] <= ins_markov_df_no_1bp['bonferroni_cutoff'])), 
            'variant_id'].unique().shape[0]
        TN = samp_ins_var.shape[0] - FP

        TP = snap_markov_1MH_SD_df_no_1bp.loc[
            (snap_markov_1MH_SD_df_no_1bp['variant_id'].isin(samp_snap_var) & 
            (snap_markov_1MH_SD_df_no_1bp['SD_snap_back_p_val'] <= snap_markov_1MH_SD_df_no_1bp['bonferroni_cutoff'])), 
                'variant_id'].unique().shape[0]
        FN =  samp_snap_var.shape[0] - TP

        tmp_boot_df.loc[idx, 'FPR'] = round((FP/(FP+TN)), 5)
        tmp_boot_df.loc[idx, 'TPR'] = round((TP/(FN+TP)), 5)

        tmp_boot_df.loc[idx,'bootN'] = b
        tmp_boot_df.loc[idx,'cutoff'] = round(cutoff,4)
        
    boot_df_no_1bp = pd.concat([boot_df_no_1bp, tmp_boot_df])
    boot_df_no_1bp.reset_index(inplace=True, drop=True)
    boot_df_no_1bp.sort_values(by=['bootN', 'cutoff'], inplace=True)


boot_df_no_1bp = boot_df_no_1bp.groupby(by=['cutoff']).mean().reset_index()
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]), ['FPR', 'TPR']] = 0,0
boot_df_no_1bp.loc[(boot_df_no_1bp.shape[0]+1), ['FPR', 'TPR']] = 1,1
boot_df_no_1bp['ins_len'] = 'insertions length > 1'
chosen_cutoff_1bp = boot_df_no_1bp.loc[
    ((boot_df_no_1bp['FPR'] < 0.05) & 
    (boot_df_no_1bp['FPR'] != 0)), :].groupby(by=['bootN']).max()['cutoff'].min()

boot_df = pd.concat([boot_df, boot_df_no_1bp])
boot_df.reset_index(inplace=True, drop=True)

ax = plt.figure(figsize = (7,5), dpi=80)
fig = ax.get_figure()
ax = sns.lineplot(data=boot_df, x='FPR', y='TPR', hue='ins_len')
plt.plot([0.05,0.05], [0,1.02], color='green', linestyle='dashed')
plt.plot([0, 1], [0, 1],ls="--", c=".3")

plt.text(0.055,boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'TPR'],chosen_cutoff)
plt.plot(boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'FPR'],
        boot_df.loc[((boot_df['ins_len'] == 'all lengths') &
                 (boot_df['cutoff'] == chosen_cutoff)), 'TPR'], 
                 color='black', marker='o', linewidth=2, markersize=4,alpha=0.8)
plt.text(0.055,boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'TPR'],chosen_cutoff_1bp)
plt.plot(boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'FPR'],
        boot_df.loc[((boot_df['ins_len'] == 'insertions length > 1') &
                 (boot_df['cutoff'] == chosen_cutoff_1bp)), 'TPR'], 
                 color='black', marker='o', linewidth=2, markersize=4,alpha=0.8)


ax.legend(title='Insertion length')
ax.set(ylim=(0, 1.02), xlim=(0, 1), 
    title=f'SD Snap back vs NHEJ insertions')
fig.savefig(f"{dir}/SNAP_MH1_SD1_INS_ROC.svg", bbox_inches='tight')
boot_df.to_csv(f"{dir}/SNAP_MH1_SD1_INS.tsv", sep='\t', index=False)










