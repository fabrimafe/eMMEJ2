"""
This script evaluate EMmej preformance using sims

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

date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
organism = 'drosophila'
# organism = 'arabidopsis'
if organism == 'drosophila':
    reference = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/ref_genome/ISO1_GCF_000001215.4.fa'
if organism == 'arabidopsis':
    reference = '/home/labs/alevy/guyta/guy_master_project/data/arabidopsis/TAIR10_ref_genome/TAIR10_fasta/GCF_000001735.4_TAIR10.1_genomic.fna'

if organism == 'drosophila':
    # out_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/20221102_sims/all_mechanisms_with_genomic_indel_length_dist/20221201_evaluations/bootstrapN_100_sample_size_1000_p_val_e_12'
    out_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/20221102_sims/all_mechanisms_with_genomic_indel_length_dist/20221221_evaluations/bootstrapN_100_sample_size_1000_p_val_e_12'
    # out_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/20221102_sims/all_mechanisms_with_genomic_indel_length_dist/20230104_evaluations/bootstrapN_100_sample_size_1000_p_val_e_12_wm50bp'
if organism == 'arabidopsis':
    # out_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/arabidopsis_thaliana_sims/20221204/bootstrapN_100_sample_size_950_p_val_e_12'
    out_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/arabidopsis_thaliana_sims/20221221/bootstrapN_100_sample_size_950_p_val_e_12'

boot_n = 100
samp_size = 1000
# samp_size = 950

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

RM_cols = ['CHR', 'POS', 'original_pos', 'variant_id', 'direction', 'ancestral_indel', 
            'derived_indel', 'indel_type', 'indel_len',  
            # deletions
            'del_mmej', 'del_mmej_cand', 'del_mmej_marked' ,'del_mmej_marked_on_ref',
            'del_last_dimer','del_mmej_cand_len',
            # snap-back
            'SD_snap_back', 'snap_mmej_marked', 'snap_P1', 'snap_P2','snap_mh1', 
            'snap_mh2', 'snap_repeat_pat', 'snap_repeat_pat_len',
            'snap_inv_comp_repeat_pat', 'snap_last_dimer',
            'snap_dist_between_reps',
            # loop-out
            'SD_loop_out','loop_mmej_marked', 'loop_P2',
            'loop_mh2','loop_repeat_pat', 'loop_repeat_pat_len',
            'loop_last_dimer', 'loop_dist_between_reps', 
            'del_mmej_motif_pos', 'del_mmej_freq_small_window', 'del_mmej_freq_large_window',
            'loop_motif_pos', 'loop_freq_small_window', 'loop_freq_large_window',
            'snap_motif_pos', 'snap_freq_small_window', 'snap_freq_large_window']


if organism == 'drosophila':
    # sims_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/20221101_sims_from_fabrizio/sims_with_genomic_indel_lenght_dist/20221123_df_size_2000/20221125_EMmej_res/whole_sims_files' # Drosophila melanogaster sims
    sims_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/20221101_sims_from_fabrizio/sims_with_genomic_indel_lenght_dist/20221221_df_size_2000/20221221_EMmej_res' # Drosophila melanogaster sims
    # sims_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/20221101_sims_from_fabrizio/sims_with_genomic_indel_lenght_dist/20221221_df_size_2000/20230104_EMmej_res' # Drosophila melanogaster sims with 50bp marlov window
if organism == 'arabidopsis':
    # sims_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/arabidopsis_taliana_sims_with_indel_len_dist/20221204_df_size_2000/20221204_EMmej_res' # Arabidopsis thaliana sims
    sims_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/arabidopsis_taliana_sims_with_indel_len_dist/20221221_df_size_2000/20221221_EMmej_res' # Arabidopsis thaliana sims
fol = ['original_sims_EMmej_format', 'sims_in_repeats', 'sims_not_in_repeats']
suf = ['_EMmej_format','_repeat', '_non_repeat']

# fol = ['original_sims_EMmej_format']
# suf = ['_EMmej_format']

shortest_indel = -2
MH_length = 2

mh_cutoff = 0.001
snap_ins_cutoff = 0.005
loop_ins_cutoff = 0.003
snap_ins_cutoff_no_1bp_ins = 0.01
loop_ins_cutoff_no_1bp_ins = 0.01
snap_loop_cutoff = 0.001

for f,s in zip(fol, suf):
    print(f"=========={f}===========")
    if organism == 'drosophila':
        mmej_markov_df_pth = f'{sims_pth}/mmej_{f}_with_genomic_indel_len_dist_EMmej_output/mmej_{f}_with_genomic_indel_len_dist_markov_output.tsv'
        nhej_markov_df_pth = f'{sims_pth}/nhej_{f}_with_genomic_indel_len_dist_no_1bp_deletion_EMmej_output/nhej_{f}_with_genomic_indel_len_dist_no_1bp_deletion_markov_output.tsv'
        ins_markov_df_pth = f'{sims_pth}/ins_{f}_with_genomic_indel_len_dist_EMmej_output/ins_{f}_with_genomic_indel_len_dist_markov_output.tsv'
        snap_markov_df_pth = f'{sims_pth}/snap_{f}_with_genomic_indel_len_dist_EMmej_output/snap_{f}_with_genomic_indel_len_dist_markov_output.tsv'
        loop_markov_df_pth = f'{sims_pth}/loop_{f}_with_genomic_indel_len_dist_EMmej_output/loop_{f}_with_genomic_indel_len_dist_markov_output.tsv'
        
        mmej_RM_df_pth = f'{sims_pth}/mmej_{f}_with_genomic_indel_len_dist_EMmej_output/mmej_{f}_with_genomic_indel_len_dist_RMdetector_output.tsv'
        nhej_RM_df_pth = f'{sims_pth}/nhej_{f}_with_genomic_indel_len_dist_no_1bp_deletion_EMmej_output/nhej_{f}_with_genomic_indel_len_dist_no_1bp_deletion_RMdetector_output.tsv'
        ins_RM_df_pth = f'{sims_pth}/ins_{f}_with_genomic_indel_len_dist_EMmej_output/ins_{f}_with_genomic_indel_len_dist_RMdetector_output.tsv'
        snap_RM_df_pth = f'{sims_pth}/snap_{f}_with_genomic_indel_len_dist_EMmej_output/snap_{f}_with_genomic_indel_len_dist_RMdetector_output.tsv'
        loop_RM_df_pth = f'{sims_pth}/loop_{f}_with_genomic_indel_len_dist_EMmej_output/loop_{f}_with_genomic_indel_len_dist_RMdetector_output.tsv'
    if organism == 'arabidopsis':
        mmej_markov_df_pth = f'{sims_pth}/mmej_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/mmej_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
        nhej_markov_df_pth = f'{sims_pth}/nhej_{f}_with_genomic_indel_len_dist_no_1bp_deletion_EMmej_format_EMmej_output/nhej_{f}_with_genomic_indel_len_dist_no_1bp_deletion_EMmej_format_markov_output.tsv'
        ins_markov_df_pth = f'{sims_pth}/ins_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/ins_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
        snap_markov_df_pth = f'{sims_pth}/snap_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/snap_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
        loop_markov_df_pth = f'{sims_pth}/loop_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/loop_{f}_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
        
        mmej_RM_df_pth = f'{sims_pth}/mmej_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/mmej_{f}_with_genomic_indel_len_dist_EMmej_format_RMdetector_output.tsv'
        nhej_RM_df_pth = f'{sims_pth}/nhej_{f}_with_genomic_indel_len_dist_no_1bp_deletion_EMmej_format_EMmej_output/nhej_{f}_with_genomic_indel_len_dist_no_1bp_deletion_EMmej_format_RMdetector_output.tsv'
        ins_RM_df_pth = f'{sims_pth}/ins_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/ins_{f}_with_genomic_indel_len_dist_EMmej_format_RMdetector_output.tsv'
        snap_RM_df_pth = f'{sims_pth}/snap_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/snap_{f}_with_genomic_indel_len_dist_EMmej_format_RMdetector_output.tsv'
        loop_RM_df_pth = f'{sims_pth}/loop_{f}_with_genomic_indel_len_dist_EMmej_format_EMmej_output/loop_{f}_with_genomic_indel_len_dist_EMmej_format_RMdetector_output.tsv'
    """
	Deletion EM section:
	This section uses the files with the indel length dist as genomic data to 
	create dataframes with different ratios of MMEJ:NHEJ and runs the deletion
	EM over it
	"""
    outpth_parent = f"{out_pth}/{date}_EM_res"
    output = f'MH_{f}_sample_size_{samp_size}'

    dir = os.path.join(f"{outpth_parent}/{output}")
    if not os.path.exists(dir):	
        if not os.path.exists(os.path.join(f"{outpth_parent}")):
            os.mkdir(os.path.join(f"{outpth_parent}"))
        os.mkdir(os.path.join(dir))

    for ratio in np.arange(0,1.1,0.1):
        samp_mmej = int(samp_size*round(ratio, 2))
        samp_nhej = int(samp_size - samp_mmej)
        EM_operator = ['bsub', '-q', 'new-short', 
        '-e', '/home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt',
        '-o', '/home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt',
        '-m', 'public_hosts', '-R', "rusage[mem=2000]", '-J', f'{samp_size}_{ratio}',
        "python3","/home/labs/alevy/guyta/guy_master_project/scripts/EMmej/EM_operator.py", 
        '-v',f"{mmej_markov_df_pth}", 
        '-u', f"{nhej_markov_df_pth}",
        '-o',f"{dir}/{round(ratio, 2)}_MMEJ_l{MH_length}_{samp_size}_shortest_deletion_{str(shortest_indel)}.tsv", 
        '-r',f"{reference}", 
        # '-w', "150" ,
        '-bs', f"{samp_mmej}", '-bs2', f"{samp_nhej}",
        '-b', f"{boot_n-1}", '-e', '12'] 
        Popen(EM_operator)

    mmej_RM_df = pd.read_csv(mmej_RM_df_pth, sep='\t',  header=0)
    nhej_RM_df = pd.read_csv(nhej_RM_df_pth, sep='\t',  header=0)
    mmej_markov_df = pd.read_csv(mmej_markov_df_pth, sep='\t')
    nhej_markov_df = pd.read_csv(nhej_markov_df_pth,sep='\t')

    ratios = np.arange(0,1.1,0.1)
    boot_df = pd.DataFrame(columns=['bootN', 'ratio','MMEJ_cut_prop',
        'NHEJ_cut_prop', 'MMEJ_count_prop',
                        'NHEJ_count_prop', 'mmej_FPR', 'mmej_FNR'])
    
    cutoff = mh_cutoff
    for b in range(boot_n):
        tmp_boot_df = pd.DataFrame(columns=['bootN', 'ratio', 'MMEJ_cut_prop', 
                        'NHEJ_cut_prop', 'MMEJ_count_prop',
                        'NHEJ_count_prop', 'mmej_FPR', 'mmej_FNR', 'ACC'], 
                    index=[i for i in range(len(ratios))])
        for idx,ratio in enumerate(ratios):
            # print(ratio)
            samp_mmej_var = pd.Series(mmej_RM_df['variant_id'].unique()).sample(
            int(round(samp_size*ratio, 0)), replace=False)
            samp_nhej_var = pd.Series(nhej_RM_df['variant_id'].unique()).sample(
                    int(round(samp_size*(1-ratio), 0)), replace=False)

            mmej_RM_c = mmej_RM_df.loc[
                (mmej_RM_df['variant_id'].isin(samp_mmej_var) & 
                (mmej_RM_df['del_mmej'] == True)), 
                    'variant_id'].unique().shape[0]
            mmej_markov_df['bonferroni_cutoff'] = cutoff / mmej_markov_df['variant_id_N']
            mmej_markov_cutoff = mmej_markov_df.loc[
                (mmej_markov_df['variant_id'].isin(samp_mmej_var) & 
                (mmej_markov_df['del_mmej_p_val'] <= mmej_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]

            nhej_RM_c = nhej_RM_df.loc[
                (nhej_RM_df['variant_id'].isin(samp_nhej_var) & 
                (nhej_RM_df['del_mmej'] == True)), 
                    'variant_id'].unique().shape[0]
            nhej_markov_df['bonferroni_cutoff'] = cutoff / nhej_markov_df['variant_id_N']
            nhej_markov_cutoff = nhej_markov_df.loc[
                (nhej_markov_df['variant_id'].isin(samp_nhej_var) & 
                (nhej_markov_df['del_mmej_p_val'] <= nhej_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
                        
            prop_cutoff = (mmej_markov_cutoff+nhej_markov_cutoff)/samp_size
            prop_count = (mmej_RM_c+nhej_RM_c)/samp_size 

            tmp_boot_df.loc[idx, 'ratio'] = round(ratio, 1)
            tmp_boot_df.loc[idx, 'bootN'] = b
            tmp_boot_df.loc[idx, 'MMEJ_cut_prop'] = round(prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'NHEJ_cut_prop'] = round(1-prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'MMEJ_count_prop'] = round(prop_count, 6)
            tmp_boot_df.loc[idx, 'NHEJ_count_prop'] = round(1-prop_count, 6)

            FP = nhej_markov_df.loc[
                (nhej_markov_df['variant_id'].isin(samp_nhej_var) & 
                (nhej_markov_df['del_mmej_p_val'] <= nhej_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
            TN = samp_nhej_var.shape[0] - FP
            
            TP = mmej_markov_df.loc[
                (mmej_markov_df['variant_id'].isin(samp_mmej_var) & 
                (mmej_markov_df['del_mmej_p_val'] <= mmej_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
            FN =  samp_mmej_var.shape[0] - TP

            if ratio == 1.0: tmp_boot_df.loc[idx, 'mmej_FPR'] = 0
            else: tmp_boot_df.loc[idx, 'mmej_FPR'] = round((FP/(FP+TN)), 4)
            if ratio == 0.0: tmp_boot_df.loc[idx, 'mmej_FNR'] = 0
            else: tmp_boot_df.loc[idx, 'mmej_FNR'] = round((FN/(FN+TP)), 4)
            tmp_boot_df.loc[idx, 'ACC'] = round(((TP+TN)/(TP+FP+TN+FN)), 4)

        boot_df = pd.concat([boot_df, tmp_boot_df])
        boot_df.reset_index(inplace=True, drop=True)

    boot_df.to_csv(
            f"{out_pth}/{f}_MMEJ_NHEJ_proportions_in_using_cutoff_{cutoff}_bootN_{boot_n}_sampsize_{samp_size}_l{MH_length}_shortest_{str(shortest_indel*(-1))}.tsv",
            sep='\t', index=False)


    # Evaluating EMmej performance in insertions
    # Snap back vs unclassified ins
    snap_RM_df = pd.read_csv(snap_RM_df_pth, sep='\t',  header=0)
    ins_RM_df = pd.read_csv(ins_RM_df_pth, sep='\t',  header=0)
    snap_markov_df = pd.read_csv(snap_markov_df_pth, sep='\t')
    ins_markov_df = pd.read_csv(ins_markov_df_pth,sep='\t')
    
  

    boot_df = pd.DataFrame(columns=['bootN', 'ratio', 'snap_cut_prop', 
        'unclassified_ins_cut_prop', 'snap_count_prop',
                        'unclassified_ins_count_prop', 'snap_FPR', 'snap_FNR', 'ACC'])
    cutoff = snap_ins_cutoff
    # b = 1
    for b in range(boot_n):
        tmp_boot_df = pd.DataFrame(columns=boot_df.columns, 
                    index=[i for i in range(len(ratios))])
        for idx,ratio in enumerate(ratios):
            samp_snap_var = pd.Series(snap_RM_df['variant_id'].unique()).sample(
            int(round(samp_size*ratio, 0)), replace=False)
            samp_ins_var = pd.Series(ins_RM_df['variant_id'].unique()).sample(
                    int(round(samp_size*(1-ratio), 0)), replace=False)
            
            snap_RM_c = snap_RM_df.loc[
                (snap_RM_df['variant_id'].isin(samp_snap_var) & 
                (snap_RM_df['SD_snap_back'] == True)), 
                    'variant_id'].unique().shape[0]
            snap_markov_df['bonferroni_cutoff'] = cutoff / snap_markov_df['variant_id_N']
            
            snap_markov_cutoff = snap_markov_df.loc[
                (snap_markov_df['variant_id'].isin(samp_snap_var) & 
                (snap_markov_df['SD_snap_back_p_val'] <= snap_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]

            ins_RM_c = ins_RM_df.loc[
                (ins_RM_df['variant_id'].isin(samp_ins_var) & 
                (ins_RM_df['SD_snap_back'] == True)), 
                    'variant_id'].unique().shape[0]
            ins_markov_df['bonferroni_cutoff'] = cutoff / ins_markov_df['variant_id_N']
                                   
            
            ins_markov_cutoff = ins_markov_df.loc[
                (ins_markov_df['variant_id'].isin(samp_ins_var) & 
                (ins_markov_df['SD_snap_back_p_val'] <= ins_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
                        
            prop_cutoff = (snap_markov_cutoff+ins_markov_cutoff)/samp_size
            prop_count = (snap_RM_c+ins_RM_c)/samp_size 

            tmp_boot_df.loc[idx, 'ratio'] = round(ratio, 1)
            tmp_boot_df.loc[idx, 'bootN'] = b
            tmp_boot_df.loc[idx, 'snap_cut_prop'] = round(prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'unclassified_ins_cut_prop'] = round(1-prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'snap_count_prop'] = round(prop_count, 6)
            tmp_boot_df.loc[idx, 'unclassified_ins_count_prop'] = round(1-prop_count, 6)

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

            if ratio == 1.0: tmp_boot_df.loc[idx, 'snap_FPR'] = 0
            else: tmp_boot_df.loc[idx, 'snap_FPR'] = round((FP/(FP+TN)), 4)
            if ratio == 0.0: tmp_boot_df.loc[idx, 'snap_FNR'] = 0
            else: tmp_boot_df.loc[idx, 'snap_FNR'] = round((FN/(FN+TP)), 4)
            tmp_boot_df.loc[idx, 'ACC'] = round(((TP+TN)/(TP+FP+TN+FN)), 4)

        boot_df = pd.concat([boot_df, tmp_boot_df])
        boot_df.reset_index(inplace=True, drop=True)

    boot_df.to_csv(
            f"{out_pth}/{f}_snap_ins_proportions_in_using_cutoff_{cutoff}_bootN_{boot_n}_sampsize_{samp_size}_l{MH_length}.tsv",
            sep='\t', index=False)
    

    # Snap back vs unclassified ins with 1bp insertions excluded
    snap_RM_df_no_1bp_ins = snap_RM_df.loc[snap_RM_df['derived_indel'].str.len()>2,:].copy()
    ins_RM_df_no_1bp_ins = ins_RM_df.loc[ins_RM_df['derived_indel'].str.len()>2,:].copy()
    snap_markov_df_no_1bp_ins = snap_markov_df.loc[snap_markov_df['derived_indel'].str.len()>2,:].copy()
    ins_markov_df_no_1bp_ins = ins_markov_df.loc[ins_markov_df['derived_indel'].str.len()>2,:].copy()
    # snap_RM_df_no_1bp_ins = snap_RM_df.loc[snap_RM_df['snap_repeat_pat_len'].str.len()>2,:].copy()
    # ins_RM_df_no_1bp_ins = ins_RM_df.loc[ins_RM_df['snap_repeat_pat_len'].str.len()>2,:].copy()
    # snap_markov_df_no_1bp_ins = snap_markov_df.loc[snap_markov_df['snap_repeat_pat_len'].str.len()>2,:].copy()
    # ins_markov_df_no_1bp_ins = ins_markov_df.loc[ins_markov_df['snap_repeat_pat_len'].str.len()>2,:].copy()
    
    boot_df = pd.DataFrame(columns=['bootN', 'ratio', 'snap_cut_prop', 
        'unclassified_ins_cut_prop', 'snap_count_prop',
                        'unclassified_ins_count_prop', 'snap_FPR', 'snap_FNR', 'ACC'])
    cutoff = snap_ins_cutoff_no_1bp_ins
    # b = 1
    for b in range(boot_n):
        tmp_boot_df = pd.DataFrame(columns=boot_df.columns, 
                    index=[i for i in range(len(ratios))])
        for idx,ratio in enumerate(ratios):
            samp_snap_var = pd.Series(snap_RM_df_no_1bp_ins['variant_id'].unique()).sample(
            int(round(samp_size*ratio, 0)), replace=False)
            samp_ins_var = pd.Series(ins_RM_df_no_1bp_ins['variant_id'].unique()).sample(
                    int(round(samp_size*(1-ratio), 0)), replace=False)
            
            snap_RM_c = snap_RM_df_no_1bp_ins.loc[
                (snap_RM_df_no_1bp_ins['variant_id'].isin(samp_snap_var) & 
                (snap_RM_df_no_1bp_ins['SD_snap_back'] == True)), 
                    'variant_id'].unique().shape[0]
            snap_markov_df_no_1bp_ins['bonferroni_cutoff'] = cutoff / snap_markov_df_no_1bp_ins['variant_id_N']
            
            snap_markov_cutoff = snap_markov_df_no_1bp_ins.loc[
                (snap_markov_df_no_1bp_ins['variant_id'].isin(samp_snap_var) & 
                (snap_markov_df_no_1bp_ins['SD_snap_back_p_val'] <= snap_markov_df_no_1bp_ins['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]

            ins_RM_c = ins_RM_df_no_1bp_ins.loc[
                (ins_RM_df_no_1bp_ins['variant_id'].isin(samp_ins_var) & 
                (ins_RM_df_no_1bp_ins['SD_snap_back'] == True)), 
                    'variant_id'].unique().shape[0]
            ins_markov_df_no_1bp_ins['bonferroni_cutoff'] = cutoff / ins_markov_df_no_1bp_ins['variant_id_N']
                                   
            
            ins_markov_cutoff = ins_markov_df_no_1bp_ins.loc[
                (ins_markov_df_no_1bp_ins['variant_id'].isin(samp_ins_var) & 
                (ins_markov_df_no_1bp_ins['SD_snap_back_p_val'] <= ins_markov_df_no_1bp_ins['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
                        
            prop_cutoff = (snap_markov_cutoff+ins_markov_cutoff)/samp_size
            prop_count = (snap_RM_c+ins_RM_c)/samp_size 

            tmp_boot_df.loc[idx, 'ratio'] = round(ratio, 1)
            tmp_boot_df.loc[idx, 'bootN'] = b
            tmp_boot_df.loc[idx, 'snap_cut_prop'] = round(prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'unclassified_ins_cut_prop'] = round(1-prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'snap_count_prop'] = round(prop_count, 6)
            tmp_boot_df.loc[idx, 'unclassified_ins_count_prop'] = round(1-prop_count, 6)

            FP = ins_markov_df_no_1bp_ins.loc[
                (ins_markov_df_no_1bp_ins['variant_id'].isin(samp_ins_var) & 
                (ins_markov_df_no_1bp_ins['SD_snap_back_p_val'] <= ins_markov_df_no_1bp_ins['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
            TN = samp_ins_var.shape[0] - FP
            
            TP = snap_markov_df_no_1bp_ins.loc[
                (snap_markov_df_no_1bp_ins['variant_id'].isin(samp_snap_var) & 
                (snap_markov_df_no_1bp_ins['SD_snap_back_p_val'] <= snap_markov_df_no_1bp_ins['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
            FN =  samp_snap_var.shape[0] - TP

            if ratio == 1.0: tmp_boot_df.loc[idx, 'snap_FPR'] = 0
            else: tmp_boot_df.loc[idx, 'snap_FPR'] = round((FP/(FP+TN)), 4)
            if ratio == 0.0: tmp_boot_df.loc[idx, 'snap_FNR'] = 0
            else: tmp_boot_df.loc[idx, 'snap_FNR'] = round((FN/(FN+TP)), 4)
            tmp_boot_df.loc[idx, 'ACC'] = round(((TP+TN)/(TP+FP+TN+FN)), 4)

        boot_df = pd.concat([boot_df, tmp_boot_df])
        boot_df.reset_index(inplace=True, drop=True)

    boot_df.to_csv(
            f"{out_pth}/{f}_snap_ins_proportions_in_using_cutoff_{cutoff}_bootN_{boot_n}_sampsize_{samp_size}_l{MH_length}_no_1bp.tsv",
            sep='\t', index=False)


    # Loop out vs unclassified
    loop_RM_df = pd.read_csv(loop_RM_df_pth, sep='\t',  header=0)
    loop_markov_df = pd.read_csv(loop_markov_df_pth, sep='\t')
    cutoff = loop_ins_cutoff
    boot_df = pd.DataFrame(columns=['bootN', 'ratio', 'loop_cut_prop', 
        'unclassified_ins_cut_prop', 'loop_count_prop',
                        'unclassified_ins_count_prop', 
                        'loop_FPR', 'loop_FNR', 'ACC'])
    for b in range(boot_n):
        tmp_boot_df = pd.DataFrame(columns=boot_df.columns, 
                    index=[i for i in range(len(ratios))])
        for idx,ratio in enumerate(ratios):
            samp_loop_var = pd.Series(loop_RM_df['variant_id'].unique()).sample(
            int(round(samp_size*ratio, 0)), replace=False)
            samp_ins_var = pd.Series(ins_RM_df['variant_id'].unique()).sample(
                    int(round(samp_size*(1-ratio), 0)), replace=False)

            loop_RM_c = loop_RM_df.loc[
                (loop_RM_df['variant_id'].isin(samp_loop_var) & 
                (loop_RM_df['SD_loop_out'] == True)), 
                    'variant_id'].unique().shape[0]
            loop_markov_df['bonferroni_cutoff'] = cutoff / loop_markov_df['variant_id_N']
            loop_markov_cutoff = loop_markov_df.loc[
                (loop_markov_df['variant_id'].isin(samp_loop_var) & 
                (loop_markov_df['SD_loop_out_p_val'] <= loop_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]

            ins_RM_c = ins_RM_df.loc[
                (ins_RM_df['variant_id'].isin(samp_ins_var) & 
                (ins_RM_df['SD_loop_out'] == True)), 
                    'variant_id'].unique().shape[0]
            ins_markov_df['bonferroni_cutoff'] = cutoff / ins_markov_df['variant_id_N']
            ins_markov_cutoff = ins_markov_df.loc[
                (ins_markov_df['variant_id'].isin(samp_ins_var) & 
                (ins_markov_df['SD_loop_out_p_val'] <= ins_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
                        
            prop_cutoff = (loop_markov_cutoff+ins_markov_cutoff)/samp_size
            prop_count = (loop_RM_c+ins_RM_c)/samp_size 

            tmp_boot_df.loc[idx, 'ratio'] = round(ratio, 1)
            tmp_boot_df.loc[idx, 'bootN'] = b
            tmp_boot_df.loc[idx, 'loop_cut_prop'] = round(prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'unclassified_ins_cut_prop'] = round(1-prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'loop_count_prop'] = round(prop_count, 6)
            tmp_boot_df.loc[idx, 'unclassified_ins_count_prop'] = round(1-prop_count, 6)

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

            if ratio == 1.0: tmp_boot_df.loc[idx, 'loop_FPR'] = 0
            else: tmp_boot_df.loc[idx, 'loop_FPR'] = round((FP/(FP+TN)), 4)
            if ratio == 0.0: tmp_boot_df.loc[idx, 'loop_FNR'] = 0
            else: tmp_boot_df.loc[idx, 'loop_FNR'] = round((FN/(FN+TP)), 4)
            tmp_boot_df.loc[idx, 'ACC'] = round(((TP+TN)/(TP+FP+TN+FN)), 4)

        boot_df = pd.concat([boot_df, tmp_boot_df])
        boot_df.reset_index(inplace=True, drop=True)

    boot_df.to_csv(
            f"{out_pth}/{f}_loop_ins_proportions_in_using_cutoff_{cutoff}_bootN_{boot_n}_sampsize_{samp_size}_l{MH_length}.tsv",
            sep='\t', index=False)


    # Loop out vs unclassified # Without 1bp insertions
    loop_RM_df_no_1bp_ins = loop_RM_df.loc[(loop_RM_df['derived_indel'].str.len() > 2), :].copy()
    loop_markov_df_no_1bp_ins = loop_markov_df.loc[(loop_markov_df['derived_indel'].str.len() > 2), :].copy()
    ins_RM_df_no_1bp_ins = ins_RM_df.loc[ins_RM_df['derived_indel'].str.len()>2,:].copy()
    ins_markov_df_no_1bp_ins = ins_markov_df.loc[ins_markov_df['derived_indel'].str.len()>2,:].copy()

    # loop_RM_df_no_1bp_ins = loop_RM_df.loc[(loop_RM_df['loop_repeat_pat_len'] > 5), :].copy()
    # loop_markov_df_no_1bp_ins = loop_markov_df.loc[(loop_markov_df['loop_repeat_pat_len'] > 5), :].copy()
    # ins_RM_df_no_1bp_ins = ins_RM_df
    # ins_markov_df_no_1bp_ins = ins_markov_df
    # ins_RM_df_no_1bp_ins = ins_RM_df_no_1bp_ins.loc[ins_RM_df_no_1bp_ins['loop_repeat_pat_len']<5,'SD_loop_out_p_val']=1
    # ins_markov_df_no_1bp_ins = ins_markov_df_no_1bp_ins.loc[ins_markov_df_no_1bp_ins['loop_repeat_pat_len']<5,'SD_loop_out_p_val']=1
    cutoff = loop_ins_cutoff_no_1bp_ins
    boot_df = pd.DataFrame(columns=['bootN', 'ratio', 'loop_cut_prop', 
        'unclassified_ins_cut_prop', 'loop_count_prop',
                        'unclassified_ins_count_prop', 
                        'loop_FPR', 'loop_FNR', 'ACC'])
    for b in range(boot_n):
        tmp_boot_df = pd.DataFrame(columns=boot_df.columns, 
                    index=[i for i in range(len(ratios))])
        for idx,ratio in enumerate(ratios):
            samp_loop_var = pd.Series(loop_RM_df_no_1bp_ins['variant_id'].unique()).sample(
            int(round(samp_size*ratio, 0)), replace=False)
            samp_ins_var = pd.Series(ins_RM_df_no_1bp_ins['variant_id'].unique()).sample(
                    int(round(samp_size*(1-ratio), 0)), replace=False)

            loop_RM_c = loop_RM_df_no_1bp_ins.loc[
                (loop_RM_df_no_1bp_ins['variant_id'].isin(samp_loop_var) & 
                (loop_RM_df_no_1bp_ins['SD_loop_out'] == True)), 
                    'variant_id'].unique().shape[0]
            loop_markov_df_no_1bp_ins['bonferroni_cutoff'] = cutoff / loop_markov_df_no_1bp_ins['variant_id_N']
            loop_markov_cutoff = loop_markov_df_no_1bp_ins.loc[
                (loop_markov_df_no_1bp_ins['variant_id'].isin(samp_loop_var) & 
                (loop_markov_df_no_1bp_ins['SD_loop_out_p_val'] <= loop_markov_df_no_1bp_ins['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]

            ins_RM_c = ins_RM_df_no_1bp_ins.loc[
                (ins_RM_df_no_1bp_ins['variant_id'].isin(samp_ins_var) & 
                (ins_RM_df_no_1bp_ins['SD_loop_out'] == True)), 
                    'variant_id'].unique().shape[0]
            ins_markov_df_no_1bp_ins['bonferroni_cutoff'] = cutoff / ins_markov_df_no_1bp_ins['variant_id_N']
            ins_markov_cutoff = ins_markov_df_no_1bp_ins.loc[
                (ins_markov_df_no_1bp_ins['variant_id'].isin(samp_ins_var) & 
                (ins_markov_df_no_1bp_ins['SD_loop_out_p_val'] <= ins_markov_df_no_1bp_ins['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
                        
            prop_cutoff = (loop_markov_cutoff+ins_markov_cutoff)/samp_size
            prop_count = (loop_RM_c+ins_RM_c)/samp_size 

            tmp_boot_df.loc[idx, 'ratio'] = round(ratio, 1)
            tmp_boot_df.loc[idx, 'bootN'] = b
            tmp_boot_df.loc[idx, 'loop_cut_prop'] = round(prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'unclassified_ins_cut_prop'] = round(1-prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'loop_count_prop'] = round(prop_count, 6)
            tmp_boot_df.loc[idx, 'unclassified_ins_count_prop'] = round(1-prop_count, 6)

            FP = ins_markov_df_no_1bp_ins.loc[
                (ins_markov_df_no_1bp_ins['variant_id'].isin(samp_ins_var) & 
                (ins_markov_df_no_1bp_ins['SD_loop_out_p_val'] <= ins_markov_df_no_1bp_ins['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
            TN = samp_ins_var.shape[0] - FP
            
            TP = loop_markov_df_no_1bp_ins.loc[
                (loop_markov_df_no_1bp_ins['variant_id'].isin(samp_loop_var) & 
                (loop_markov_df_no_1bp_ins['SD_loop_out_p_val'] <= loop_markov_df_no_1bp_ins['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
            FN =  samp_loop_var.shape[0] - TP

            if ratio == 1.0: tmp_boot_df.loc[idx, 'loop_FPR'] = 0
            else: tmp_boot_df.loc[idx, 'loop_FPR'] = round((FP/(FP+TN)), 4)
            if ratio == 0.0: tmp_boot_df.loc[idx, 'loop_FNR'] = 0
            else: tmp_boot_df.loc[idx, 'loop_FNR'] = round((FN/(FN+TP)), 4)
            tmp_boot_df.loc[idx, 'ACC'] = round(((TP+TN)/(TP+FP+TN+FN)), 4)

        boot_df = pd.concat([boot_df, tmp_boot_df])
        boot_df.reset_index(inplace=True, drop=True)

    boot_df.to_csv(
            f"{out_pth}/{f}_loop_ins_proportions_in_using_cutoff_{cutoff}_bootN_{boot_n}_sampsize_{samp_size}_l{MH_length}_witn_no_1bp.tsv",
            sep='\t', index=False)

    # Snap back vs Loop out
    boot_df = pd.DataFrame(columns=['bootN', 'ratio', 'loop_cut_prop', 
        'snap_cut_prop', 'loop_count_prop',
                        'snap_count_prop', 
                        'snap_FPR', 'snap_FNR', 'ACC'])
    cutoff = snap_loop_cutoff
    for b in range(boot_n):
        tmp_boot_df = pd.DataFrame(columns=boot_df.columns, 
                    index=[i for i in range(len(ratios))])
        for idx,ratio in enumerate(ratios):
            samp_snap_var = pd.Series(snap_RM_df['variant_id'].unique()).sample(
                    int(round(samp_size*(ratio), 0)), replace=False)
            samp_loop_var = pd.Series(loop_RM_df['variant_id'].unique()).sample(
            int(round(samp_size*(1-ratio), 0)), replace=False)
            
            snap_RM_c = snap_RM_df.loc[
                (snap_RM_df['variant_id'].isin(samp_snap_var) & 
                (snap_RM_df['SD_snap_back'] == True)), 
                    'variant_id'].unique().shape[0]
            snap_markov_df['bonferroni_cutoff'] = cutoff / snap_markov_df['variant_id_N']
            snap_markov_cutoff = snap_markov_df.loc[
                (snap_markov_df['variant_id'].isin(samp_snap_var) & 
                (snap_markov_df['SD_snap_back_p_val'] <= snap_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]

            loop_RM_c = loop_RM_df.loc[
                (loop_RM_df['variant_id'].isin(samp_loop_var) & 
                (loop_RM_df['SD_snap_back'] == True)), 
                    'variant_id'].unique().shape[0]
            loop_markov_df['bonferroni_cutoff'] = cutoff / loop_markov_df['variant_id_N']
            loop_markov_cutoff = loop_markov_df.loc[
                (loop_markov_df['variant_id'].isin(samp_loop_var) & 
                (loop_markov_df['SD_snap_back_p_val'] <= loop_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
                        
            prop_cutoff = (snap_markov_cutoff+loop_markov_cutoff)/samp_size
            prop_count = (snap_RM_c+loop_RM_c)/samp_size 

            tmp_boot_df.loc[idx, 'ratio'] = round(ratio, 1)
            tmp_boot_df.loc[idx, 'bootN'] = b
            tmp_boot_df.loc[idx, 'loop_cut_prop'] = round(prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'snap_cut_prop'] = round(1-prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'loop_count_prop'] = round(prop_count, 6)
            tmp_boot_df.loc[idx, 'snap_count_prop'] = round(1-prop_count, 6)

            FP = loop_markov_df.loc[
                (loop_markov_df['variant_id'].isin(samp_loop_var) & 
                (loop_markov_df['SD_snap_back_p_val'] <= loop_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
            TN = samp_loop_var.shape[0] - FP
            
            TP = snap_markov_df.loc[
                (snap_markov_df['variant_id'].isin(samp_snap_var) & 
                (snap_markov_df['SD_snap_back_p_val'] <= snap_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
            FN =  samp_snap_var.shape[0] - TP

            if ratio == 1.0: tmp_boot_df.loc[idx, 'snap_FPR'] = 0
            else: tmp_boot_df.loc[idx, 'snap_FPR'] = round((FP/(FP+TN)), 4)
            if ratio == 0.0: tmp_boot_df.loc[idx, 'snap_FNR'] = 0
            else: tmp_boot_df.loc[idx, 'snap_FNR'] = round((FN/(FN+TP)), 4)
            tmp_boot_df.loc[idx, 'ACC'] = round(((TP+TN)/(TP+FP+TN+FN)), 4)

        boot_df = pd.concat([boot_df, tmp_boot_df])
        boot_df.reset_index(inplace=True, drop=True)

    boot_df.to_csv(
            f"{out_pth}/{f}_snap_loop_proportions_in_using_cutoff_{cutoff}_bootN_{boot_n}_sampsize_{samp_size}_l{MH_length}.tsv",
            sep='\t', index=False)


    # Loop out vs Snap back
    boot_df = pd.DataFrame(columns=['bootN', 'ratio', 'loop_cut_prop', 
        'snap_cut_prop', 'loop_count_prop',
                        'snap_count_prop', 
                        'loop_FPR', 'loop_FNR', 'ACC'])
    cutoff = snap_loop_cutoff
    for b in range(boot_n):
        tmp_boot_df = pd.DataFrame(columns=boot_df.columns, 
                    index=[i for i in range(len(ratios))])
        for idx,ratio in enumerate(ratios):
            samp_loop_var = pd.Series(loop_RM_df['variant_id'].unique()).sample(
                    int(round(samp_size*ratio, 0)), replace=False)
            samp_snap_var = pd.Series(snap_RM_df['variant_id'].unique()).sample(
                    int(round(samp_size*(1-ratio), 0)), replace=False)

            snap_RM_c = snap_RM_df.loc[
                (snap_RM_df['variant_id'].isin(samp_snap_var) & 
                (snap_RM_df['SD_loop_out'] == True)), 
                    'variant_id'].unique().shape[0]
            snap_markov_df['bonferroni_cutoff'] = cutoff / snap_markov_df['variant_id_N']
            snap_markov_cutoff = snap_markov_df.loc[
                (snap_markov_df['variant_id'].isin(samp_snap_var) & 
                (snap_markov_df['SD_loop_out_p_val'] <= snap_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]

            loop_RM_c = loop_RM_df.loc[
                (loop_RM_df['variant_id'].isin(samp_loop_var) & 
                (loop_RM_df['SD_loop_out'] == True)), 
                    'variant_id'].unique().shape[0]
            loop_markov_df['bonferroni_cutoff'] = cutoff / loop_markov_df['variant_id_N']
            loop_markov_cutoff = loop_markov_df.loc[
                (loop_markov_df['variant_id'].isin(samp_loop_var) & 
                (loop_markov_df['SD_loop_out_p_val'] <= loop_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
                        
            prop_cutoff = (snap_markov_cutoff+loop_markov_cutoff)/samp_size
            prop_count = (snap_RM_c+loop_RM_c)/samp_size 

            TP = loop_markov_df.loc[
                (loop_markov_df['variant_id'].isin(samp_loop_var) & 
                (loop_markov_df['SD_loop_out_p_val'] <= loop_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
            FN = samp_loop_var.shape[0] - TP
            
            FP = snap_markov_df.loc[
                (snap_markov_df['variant_id'].isin(samp_snap_var) & 
                (snap_markov_df['SD_loop_out_p_val'] <= snap_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
            TN =  samp_snap_var.shape[0] - FP

            if ratio == 1.0: tmp_boot_df.loc[idx, 'loop_FPR'] = 0
            else: tmp_boot_df.loc[idx, 'loop_FPR'] = round((FP/(FP+TN)), 4)
            if ratio == 0.0: tmp_boot_df.loc[idx, 'loop_FNR'] = 0
            else: tmp_boot_df.loc[idx, 'loop_FNR'] = round((FN/(FN+TP)), 4)
            tmp_boot_df.loc[idx, 'ACC'] = round(((TP+TN)/(TP+FP+TN+FN)), 4)


            tmp_boot_df.loc[idx, 'ratio'] = round(ratio, 1)
            tmp_boot_df.loc[idx, 'bootN'] = b
            tmp_boot_df.loc[idx, 'loop_cut_prop'] = round(prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'snap_cut_prop'] = round(1-prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'loop_count_prop'] = round(prop_count, 6)
            tmp_boot_df.loc[idx, 'snap_count_prop'] = round(1-prop_count, 6)

        boot_df = pd.concat([boot_df, tmp_boot_df])
        boot_df.reset_index(inplace=True, drop=True)
    boot_df.to_csv(
            f"{out_pth}/{f}_loop_snap_proportions_in_using_cutoff_{cutoff}_bootN_{boot_n}_sampsize_{samp_size}_l{MH_length}.tsv",
            sep='\t', index=False)


# out_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/20221102_sims/all_mechanisms_with_genomic_indel_length_dist/20221123_evaluations/bootstrapN_100_sample_size_500'
# Evaluating sims with MH=1,4 seperatly:
f = 'original_sims_EMmej_format'
for i in [1,4]:
    if organism == 'arabidopsis':
        mmej_markov_df_pth = f'{sims_pth}/mmej_MH{i}_original_sims_EMmej_format_with_genomic_indel_len_dist_EMmej_format_EMmej_output/mmej_MH{i}_original_sims_EMmej_format_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
        nhej_markov_df_pth = f'{sims_pth}/nhej_original_sims_EMmej_format_with_genomic_indel_len_dist_EMmej_format_EMmej_output/nhej_original_sims_EMmej_format_with_genomic_indel_len_dist_EMmej_format_markov_output.tsv'
            
        mmej_RM_df_pth = f'{sims_pth}/mmej_MH{i}_original_sims_EMmej_format_with_genomic_indel_len_dist_EMmej_format_EMmej_output/mmej_MH{i}_original_sims_EMmej_format_with_genomic_indel_len_dist_EMmej_format_RMdetector_output.tsv'
        nhej_RM_df_pth = f'{sims_pth}/nhej_original_sims_EMmej_format_with_genomic_indel_len_dist_EMmej_format_EMmej_output/nhej_original_sims_EMmej_format_with_genomic_indel_len_dist_EMmej_format_RMdetector_output.tsv'
    if organism == 'drosophila':
        mmej_markov_df_pth = f'{sims_pth}/mmej_MH{i}_original_sims_EMmej_format_with_genomic_indel_len_dist_EMmej_output/mmej_MH{i}_original_sims_EMmej_format_with_genomic_indel_len_dist_markov_output.tsv'
        nhej_markov_df_pth = f'{sims_pth}/nhej_original_sims_EMmej_format_with_genomic_indel_len_dist_EMmej_output/nhej_original_sims_EMmej_format_with_genomic_indel_len_dist_markov_output.tsv'
            
        mmej_RM_df_pth = f'{sims_pth}/mmej_MH{i}_original_sims_EMmej_format_with_genomic_indel_len_dist_EMmej_output/mmej_MH{i}_original_sims_EMmej_format_with_genomic_indel_len_dist_RMdetector_output.tsv'
        nhej_RM_df_pth = f'{sims_pth}/nhej_original_sims_EMmej_format_with_genomic_indel_len_dist_EMmej_output/nhej_original_sims_EMmej_format_with_genomic_indel_len_dist_RMdetector_output.tsv'
        
    """
	Deletion EM section:
	This section uses the files with the indel length dist as genomic data to 
	create dataframes with different ratios of MMEJ:NHEJ and runs the deletion
	EM over it
	"""
    outpth_parent = f"{out_pth}/{date}_EM_res"
    output = f'MH{i}_original_sims_EMmej_format_sample_size_{samp_size}'
    dir = os.path.join(f"{outpth_parent}/{output}")
    if not os.path.exists(dir):	
        if not os.path.exists(os.path.join(f"{outpth_parent}")):
            os.mkdir(os.path.join(f"{outpth_parent}"))
        os.mkdir(os.path.join(dir))
    for ratio in np.arange(0,1.1,0.1):
        samp_mmej = int(samp_size*round(ratio, 2))
        samp_nhej = int(samp_size - samp_mmej)
        EM_operator = ['bsub', '-q', 'new-short', 
        '-e', '/home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt',
        '-o', '/home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt',
        '-m', 'public_hosts', '-R', "rusage[mem=2000]", '-J', f'{samp_size}_{ratio}',
        "python3","/home/labs/alevy/guyta/guy_master_project/scripts/EMmej/EM_operator.py", 
        '-v',f"{mmej_markov_df_pth}", 
        '-u', f"{nhej_markov_df_pth}",
        '-o',f"{dir}/{round(ratio, 2)}_MMEJ_l{i}_{samp_size}_shortest_deletion_{str(((-1)*i))}.tsv", 
        '-r',f"{reference}", 
        # '-w', "150" ,
        '-bs', f"{samp_mmej}", '-bs2', f"{samp_nhej}",
        '-b', f"{boot_n-1}", '-e', '12'] 
        Popen(EM_operator)

    mmej_RM_df = pd.read_csv(mmej_RM_df_pth, sep='\t',  header=0)
    nhej_RM_df = pd.read_csv(nhej_RM_df_pth, sep='\t',  header=0)
    mmej_markov_df = pd.read_csv(mmej_markov_df_pth, sep='\t')
    nhej_markov_df = pd.read_csv(nhej_markov_df_pth,sep='\t')

    ratios = np.arange(0,1.1,0.1)
    boot_df = pd.DataFrame(columns=['bootN', 'ratio', 
        'MMEJ_cut_prop', 'NHEJ_cut_prop', 'MMEJ_count_prop',
        'NHEJ_count_prop', 'mmej_FPR', 'mmej_FNR', 'ACC'])
    cutoff = mh_cutoff
    for b in range(boot_n):
        tmp_boot_df = pd.DataFrame(columns=['bootN', 'ratio', 
        'MMEJ_cut_prop', 'NHEJ_cut_prop', 'MMEJ_count_prop',
        'NHEJ_count_prop', 'mmej_FPR', 'mmej_FNR', 'ACC'], 
                    index=[i for i in range(len(ratios))])
        for idx,ratio in enumerate(ratios):
            samp_mmej_var = pd.Series(mmej_RM_df['variant_id'].unique()).sample(
                    int(round(samp_size*ratio, 0)), replace=False)
            samp_nhej_var = pd.Series(nhej_RM_df['variant_id'].unique()).sample(
                    int(round(samp_size*(1-ratio), 0)), replace=False)
            
            mmej_RM_c = mmej_RM_df.loc[
                (mmej_RM_df['variant_id'].isin(samp_mmej_var) & 
                (mmej_RM_df['del_mmej'] == True)), 
                    'variant_id'].unique().shape[0]
            mmej_markov_df['bonferroni_cutoff'] = cutoff / mmej_markov_df['variant_id_N']
            mmej_markov_cutoff = mmej_markov_df.loc[
                (mmej_markov_df['variant_id'].isin(samp_mmej_var) & 
                (mmej_markov_df['del_mmej_p_val']<= mmej_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]

            nhej_RM_c = nhej_RM_df.loc[
                (nhej_RM_df['variant_id'].isin(samp_nhej_var) & 
                (nhej_RM_df['del_mmej'] == True)), 
                    'variant_id'].unique().shape[0]
            nhej_markov_df['bonferroni_cutoff'] = cutoff / nhej_markov_df['variant_id_N']
            nhej_markov_cutoff = nhej_markov_df.loc[
                (nhej_markov_df['variant_id'].isin(samp_nhej_var) & 
                (nhej_markov_df['del_mmej_p_val']<= nhej_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
                        
            prop_cutoff = (mmej_markov_cutoff+nhej_markov_cutoff)/samp_size
            prop_count = (mmej_RM_c+nhej_RM_c)/samp_size 

            tmp_boot_df.loc[idx, 'ratio'] = round(ratio, 1)
            tmp_boot_df.loc[idx, 'bootN'] = b
            tmp_boot_df.loc[idx, 'MMEJ_cut_prop'] = round(prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'NHEJ_cut_prop'] = round(1-prop_cutoff, 6)
            tmp_boot_df.loc[idx, 'MMEJ_count_prop'] = round(prop_count, 6)
            tmp_boot_df.loc[idx, 'NHEJ_count_prop'] = round(1-prop_count, 6)

            FP = nhej_markov_df.loc[
                (nhej_markov_df['variant_id'].isin(samp_nhej_var) & 
                (nhej_markov_df['del_mmej_p_val']<= nhej_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
            TN = samp_nhej_var.shape[0] - FP
            
            TP = mmej_markov_df.loc[
                (mmej_markov_df['variant_id'].isin(samp_mmej_var) & 
                (mmej_markov_df['del_mmej_p_val']<= mmej_markov_df['bonferroni_cutoff'])), 
                    'variant_id'].unique().shape[0]
            FN =  samp_mmej_var.shape[0] - TP

            if ratio == 1.0: tmp_boot_df.loc[idx, 'mmej_FPR'] = 0
            else: tmp_boot_df.loc[idx, 'mmej_FPR'] = round((FP/(FP+TN)), 4)
            if ratio == 0.0: tmp_boot_df.loc[idx, 'mmej_FNR'] = 0
            else: tmp_boot_df.loc[idx, 'mmej_FNR'] = round((FN/(FN+TP)), 4)
            tmp_boot_df.loc[idx, 'ACC'] = round(((TP+TN)/(TP+FP+TN+FN)), 4)
            
        boot_df = pd.concat([boot_df, tmp_boot_df])
        boot_df.reset_index(inplace=True, drop=True)


    boot_df.to_csv(
            f"{out_pth}/{f}_MMEJ_NHEJ_proportions_in_using_cutoff_{cutoff}_bootN_{boot_n}_sampsize_{samp_size}_l{i}_shortest_{str(i*(-1))}.tsv",
            sep='\t', index=False)











