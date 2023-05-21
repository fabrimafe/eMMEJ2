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


# os.chdir('/home/labs/alevy/guyta/guy_master_project/scripts/EMmej')
# set constants
date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
df_size = 2000
shortest_indel = -2
reference = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/ref_genome/ISO1_GCF_000001215.4.fa'
out_pth = f'/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/20221101_sims_from_fabrizio/sims_with_genomic_indel_lenght_dist/{date}_df_size_{df_size}'
print(out_pth)
dir = os.path.join(f"{out_pth}")
if not os.path.exists(dir):
    os.mkdir(os.path.join(dir))
# get indel length distributions based on 
# genomic data (Clark's data)
whole_df_pth = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/construct_ancestral_state/GDL_indel_EMmej_format.vcf.gz'
rep_df_pth = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/stratification_by_annotation/GDL_indel_EMmej_format_repeated_regions.vcf.gz'
non_rep_df_pth = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/stratification_by_annotation/GDL_indel_EMmej_format_non_repeated_regions.vcf.gz'

def get_indel_length_dist(df :pd.DataFrame) -> np.array:
    """
    A function that calculate indel length distribution in
    a given df
    """
    # calc indel length
    df['indel_len'] = df['ALT'].str.len() - df['REF'].str.len()
    return {i : (df.loc[(df['indel_len'] == i), :].shape[0] /
                df.shape[0]) for i in range(df['indel_len'].min(),df['indel_len'].max())}




def generate_df_with_indel_length_dist(
        df: pd.DataFrame, genomic_dist: dict,
        shortest_indel: int, indel_type:str):
        df = df.loc[(df['indel_len'].isin(genomic_dist.keys())), :]
        ind_lenghts = pd.Series(df['indel_len'].unique().tolist()).sort_values()
        if indel_type == 'DEL': ind_lenghts = ind_lenghts[ind_lenghts <= shortest_indel]
        norm_len_freq = 0
        for ind in ind_lenghts: norm_len_freq = norm_len_freq + genomic_dist[ind]
        out_df = pd.DataFrame(columns=['#CHR', 'POS', 'REF', 'ALT'])
        tot_mmej_samp_var_id = []
        ind_samp_size_count = 0
        for ind in ind_lenghts:
            # sample variant ids
            ind_samp_size = int(round((df_size*genomic_dist[ind]/norm_len_freq), 0))
            # print(ind , ind_samp_size)
            ind_samp_size_count = ind_samp_size_count + ind_samp_size
            samp_var_id = pd.Series(
                    df.loc[((df['indel_len'] == ind) & 
                    (~df['variant_id'].isin(tot_mmej_samp_var_id))),
                        'variant_id'].unique()).sample(ind_samp_size)
            tot_mmej_samp_var_id = tot_mmej_samp_var_id+samp_var_id.to_list()
            
            out_df = pd.concat([out_df, 
                    df.loc[df['variant_id'].isin(samp_var_id), :]])
        return out_df



cols = ['CHR', 'POS','POS_e', 'REF', 'ALT', 'ANC']
whole_df = pd.read_csv(whole_df_pth, sep='\t', 
        names=['CHR', 'POS', 'REF', 'ALT', 'ANC'], index_col=None)
non_rep_df = pd.read_csv(non_rep_df_pth, sep='\t', names=cols, index_col=None)
non_rep_df.drop(columns=['POS_e'], inplace=True)
rep_df = pd.read_csv(rep_df_pth, sep='\t', names=cols, index_col=None)
rep_df.drop(columns=['POS_e'], inplace=True)

# Get indel length distribution from genomic datasets
whole_df_indel_len_dist = get_indel_length_dist(df=whole_df)    
non_rep_indel_len_dist = get_indel_length_dist(df=non_rep_df)
rep_indel_len_dist = get_indel_length_dist(df=rep_df)

sims_pth = '/home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/20221101_sims_from_fabrizio'
cols = ['#CHR', 'POS', 'REF', 'ALT']

dists = [whole_df_indel_len_dist, rep_indel_len_dist, non_rep_indel_len_dist]
fol = ['original_sims_EMmej_format', 'sims_in_repeats', 'sims_not_in_repeats']
suf = ['_EMmej_format','_repeat', '_non_repeat']

# dists = [whole_df_indel_len_dist]
# fol = ['original_sims_EMmej_format']
# suf = ['_EMmej_format']

for f,d,s in zip(fol, dists, suf):
    print(f)
    # NHEJ
    nhej_sims = pd.DataFrame(columns=cols)
    for i in range(1,10):
        nhej_sims = pd.concat([nhej_sims, 
            pd.read_csv(f'{sims_pth}/{f}/deletion_NHEJ_l{i}_d{i}{s}.vcf', sep='\t', comment='#', names=cols)])
    nhej_sims = pd.concat([nhej_sims, 
            pd.read_csv(f'{sims_pth}/{f}/deletion_NHEJ_l10_d50{s}.vcf', sep='\t', comment='#', names=cols)])
    nhej_sims.reset_index(drop=True, inplace=True)
    nhej_sims['variant_id'] = nhej_sims.apply(lambda row: f"{row['#CHR']}_{row['POS']}", axis=1)
    nhej_sims['indel_len'] = nhej_sims.apply(lambda row: abs(len(row['REF']) - len(row['ALT'])), axis=1)
    nhej_sims['indel_len'] = nhej_sims['indel_len']*(-1)
    nhej_sims['indel_len'] = nhej_sims['indel_len'].astype(int)
    nhej_df = generate_df_with_indel_length_dist(df=nhej_sims, 
            genomic_dist=d, shortest_indel=shortest_indel, indel_type='DEL')
    nhej_df_all_indel_lengths = generate_df_with_indel_length_dist(df=nhej_sims, 
            genomic_dist=d, shortest_indel=-1, indel_type='DEL')

    # MMEJ
    mmej_sims = pd.DataFrame(columns=cols)
    mmej_sims = pd.concat([mmej_sims, 
            pd.read_csv(f'{sims_pth}/{f}/deletion_MMEJ_l2_d10{s}.vcf', sep='\t', comment='#', names=cols)])
    mmej_sims = pd.concat([mmej_sims, 
            pd.read_csv(f'{sims_pth}/{f}/deletion_MMEJ_l2_d200{s}.vcf', sep='\t', comment='#', names=cols)])
    mmej_sims.reset_index(drop=True, inplace=True)
    mmej_sims['variant_id'] = mmej_sims.apply(lambda row: f"{row['#CHR']}_{row['POS']}", axis=1)
    mmej_sims['indel_len'] = mmej_sims.apply(lambda row: abs(len(row['REF']) - len(row['ALT'])), axis=1)
    mmej_sims['indel_len'] = mmej_sims['indel_len']*(-1)
    mmej_df = generate_df_with_indel_length_dist(df=mmej_sims, 
            genomic_dist=d, shortest_indel=shortest_indel, indel_type='DEL')
    
    # Unclassified ins
    ins_sims = pd.DataFrame(columns=cols)
    for i in range(1,10):
        ins_sims = pd.concat([ins_sims, 
            pd.read_csv(f'{sims_pth}/{f}/insertion_l{i}_d{i}{s}.vcf', sep='\t', comment='#', names=cols)])
    ins_sims = pd.concat([ins_sims, 
            pd.read_csv(f'{sims_pth}/{f}/insertion_l10_d50{s}.vcf', sep='\t', comment='#', names=cols)])
    ins_sims.reset_index(drop=True, inplace=True)
    ins_sims['variant_id'] = ins_sims.apply(lambda row: f"{row['#CHR']}_{row['POS']}", axis=1)
    ins_sims['indel_len'] = ins_sims.apply(lambda row: abs(len(row['REF']) - len(row['ALT'])), axis=1)

    ins_df = generate_df_with_indel_length_dist(df=ins_sims, 
            genomic_dist=d, shortest_indel=shortest_indel, indel_type='INS')
    
    # Snap-Back
    snap_sims = pd.DataFrame(columns=cols)
    snap_sims = pd.concat([snap_sims, 
            pd.read_csv(f'{sims_pth}/{f}/ins_d140/SDsnapback_l2_d140_s2{s}.vcf', sep='\t', comment='#', names=cols)])
    # snap_sims = pd.concat([snap_sims, 
    #         pd.read_csv(f'{sims_pth}/{f}/ins_140/SDsnapback_l2_d140{s}.vcf', sep='\t', comment='#', names=cols)])
    snap_sims.reset_index(drop=True, inplace=True)
    snap_sims['variant_id'] = snap_sims.apply(lambda row: f"{row['#CHR']}_{row['POS']}", axis=1)
    snap_sims['indel_len'] = snap_sims.apply(lambda row: abs(len(row['REF']) - len(row['ALT'])), axis=1)
    snap_df = generate_df_with_indel_length_dist(df=snap_sims, 
            genomic_dist=d, shortest_indel=shortest_indel, indel_type='INS')

    # Loop-out
    loop_sims = pd.DataFrame(columns=cols)
    loop_sims = pd.concat([loop_sims, 
            pd.read_csv(f'{sims_pth}/{f}/ins_d140/SDloopout_l2_d140_s2{s}.vcf', sep='\t', comment='#', names=cols)])
    # loop_sims = pd.concat([loop_sims, 
    #         pd.read_csv(f'{sims_pth}/{f}/ins_140/SDloopout_l2_d140{s}.vcf', sep='\t', comment='#', names=cols)])
    loop_sims.reset_index(drop=True, inplace=True)
    loop_sims['variant_id'] = loop_sims.apply(lambda row: f"{row['#CHR']}_{row['POS']}", axis=1)
    loop_sims['indel_len'] = loop_sims.apply(lambda row: abs(len(row['REF']) - len(row['ALT'])), axis=1)
    loop_df = generate_df_with_indel_length_dist(df=loop_sims, 
            genomic_dist=d, shortest_indel=shortest_indel, indel_type='INS')
    
    if f == 'original_sims_EMmej_format':
        # MMEJ with MH = 1
        shortest_indel = -1
        mmej_MH1_sims = pd.DataFrame(columns=cols)
        mmej_MH1_sims = pd.concat([mmej_MH1_sims, 
                pd.read_csv(f'{sims_pth}/{f}/deletion_MMEJ_l1_d10{s}.vcf', sep='\t', comment='#', names=cols)])
        mmej_MH1_sims = pd.concat([mmej_MH1_sims, 
                pd.read_csv(f'{sims_pth}/{f}/deletion_MMEJ_l1_d200{s}.vcf', sep='\t', comment='#', names=cols)])
        mmej_MH1_sims.reset_index(drop=True, inplace=True)
        mmej_MH1_sims['variant_id'] = mmej_MH1_sims.apply(lambda row: f"{row['#CHR']}_{row['POS']}", axis=1)
        mmej_MH1_sims['indel_len'] = mmej_MH1_sims.apply(lambda row: abs(len(row['REF']) - len(row['ALT'])), axis=1)
        mmej_MH1_sims['indel_len'] = mmej_MH1_sims['indel_len']*(-1)

        mmej_MH1_df = generate_df_with_indel_length_dist(df=mmej_MH1_sims, 
                genomic_dist=d, shortest_indel=shortest_indel, indel_type='DEL')

        
        # MMEJ with MH = 4
        shortest_indel = -4
        mmej_MH4_sims = pd.DataFrame(columns=cols)
        mmej_MH4_sims = pd.concat([mmej_MH4_sims, 
                pd.read_csv(f'{sims_pth}/{f}/deletion_MMEJ_l4_d10{s}.vcf', sep='\t', comment='#', names=cols)])
        mmej_MH4_sims = pd.concat([mmej_MH4_sims, 
                pd.read_csv(f'{sims_pth}/{f}/deletion_MMEJ_l4_d200{s}.vcf', sep='\t', comment='#', names=cols)])
        mmej_MH4_sims.reset_index(drop=True, inplace=True)
        mmej_MH4_sims['variant_id'] = mmej_MH4_sims.apply(lambda row: f"{row['#CHR']}_{row['POS']}", axis=1)
        mmej_MH4_sims['indel_len'] = mmej_MH4_sims.apply(lambda row: abs(len(row['REF']) - len(row['ALT'])), axis=1)
        mmej_MH4_sims['indel_len'] = mmej_MH4_sims['indel_len']*(-1)

        mmej_MH4_df = generate_df_with_indel_length_dist(df=mmej_MH4_sims, 
                genomic_dist=d, shortest_indel=shortest_indel, indel_type='DEL')


        # Snap-Back with 1bp MH and SD
        snap_1bp_MH_SD_sims = pd.DataFrame(columns=cols)
        snap_1bp_MH_SD_sims = pd.concat([snap_1bp_MH_SD_sims, 
                pd.read_csv(f'{sims_pth}/{f}/ins_d140/SDsnapback_l1_d140_s1{s}.vcf', sep='\t', comment='#', names=cols)])
        # snap_1bp_MH_SD_sims = pd.concat([snap_1bp_MH_SD_sims, 
        #         pd.read_csv(f'{sims_pth}/{f}/ins_140/SDsnapback_l2_d140{s}.vcf', sep='\t', comment='#', names=cols)])
        snap_1bp_MH_SD_sims.reset_index(drop=True, inplace=True)
        snap_1bp_MH_SD_sims['variant_id'] = snap_1bp_MH_SD_sims.apply(lambda row: f"{row['#CHR']}_{row['POS']}", axis=1)
        snap_1bp_MH_SD_sims['indel_len'] = snap_1bp_MH_SD_sims.apply(lambda row: abs(len(row['REF']) - len(row['ALT'])), axis=1)
        snap_1bp_MH_SD_df = generate_df_with_indel_length_dist(df=snap_1bp_MH_SD_sims, 
                genomic_dist=d, shortest_indel=shortest_indel, indel_type='INS')
        
        # Snap-Back with 4bp MH and SD
        snap_4bp_MH_SD_sims = pd.DataFrame(columns=cols)
        snap_4bp_MH_SD_sims = pd.concat([snap_4bp_MH_SD_sims, 
                pd.read_csv(f'{sims_pth}/{f}/ins_d140/SDsnapback_l4_d140_s4{s}.vcf', sep='\t', comment='#', names=cols)])
        # snap_4bp_MH_SD_sims = pd.concat([snap_4bp_MH_SD_sims, 
        #         pd.read_csv(f'{sims_pth}/{f}/ins_140/SDsnapback_l2_d140{s}.vcf', sep='\t', comment='#', names=cols)])
        snap_4bp_MH_SD_sims.reset_index(drop=True, inplace=True)
        snap_4bp_MH_SD_sims['variant_id'] = snap_4bp_MH_SD_sims.apply(lambda row: f"{row['#CHR']}_{row['POS']}", axis=1)
        snap_4bp_MH_SD_sims['indel_len'] = snap_4bp_MH_SD_sims.apply(lambda row: abs(len(row['REF']) - len(row['ALT'])), axis=1)
        snap_4bp_MH_SD_df = generate_df_with_indel_length_dist(df=snap_4bp_MH_SD_sims, 
                genomic_dist=d, shortest_indel=shortest_indel, indel_type='INS')

        # Loop-out with 1bp MH and SD
        loop_1bp_MH_SD_sims = pd.DataFrame(columns=cols)
        loop_1bp_MH_SD_sims = pd.concat([loop_1bp_MH_SD_sims, 
                pd.read_csv(f'{sims_pth}/{f}/ins_d140/SDloopout_l1_d140_s1{s}.vcf', sep='\t', comment='#', names=cols)])
        # loop_1bp_MH_SD_sims = pd.concat([loop_1bp_MH_SD_sims, 
        #         pd.read_csv(f'{sims_pth}/{f}/ins_140/SDloopout_l2_d140{s}.vcf', sep='\t', comment='#', names=cols)])
        loop_1bp_MH_SD_sims.reset_index(drop=True, inplace=True)
        loop_1bp_MH_SD_sims['variant_id'] = loop_1bp_MH_SD_sims.apply(lambda row: f"{row['#CHR']}_{row['POS']}", axis=1)
        loop_1bp_MH_SD_sims['indel_len'] = loop_1bp_MH_SD_sims.apply(lambda row: abs(len(row['REF']) - len(row['ALT'])), axis=1)
        loop_1bp_MH_SD_df = generate_df_with_indel_length_dist(df=loop_1bp_MH_SD_sims, 
                genomic_dist=d, shortest_indel=shortest_indel, indel_type='INS')
        
        # Loop-out with 4bp MH and SD
        loop_4bp_MH_SD_sims = pd.DataFrame(columns=cols)
        loop_4bp_MH_SD_sims = pd.concat([loop_4bp_MH_SD_sims, 
                pd.read_csv(f'{sims_pth}/{f}/ins_d140/SDloopout_l4_d140_s4{s}.vcf', sep='\t', comment='#', names=cols)])
        # loop_4bp_MH_SD_sims = pd.concat([loop_4bp_MH_SD_sims, 
        #         pd.read_csv(f'{sims_pth}/{f}/ins_140/SDloopout_l2_d140{s}.vcf', sep='\t', comment='#', names=cols)])
        loop_4bp_MH_SD_sims.reset_index(drop=True, inplace=True)
        loop_4bp_MH_SD_sims['variant_id'] = loop_4bp_MH_SD_sims.apply(lambda row: f"{row['#CHR']}_{row['POS']}", axis=1)
        loop_4bp_MH_SD_sims['indel_len'] = loop_4bp_MH_SD_sims.apply(lambda row: abs(len(row['REF']) - len(row['ALT'])), axis=1)
        loop_4bp_MH_SD_df = generate_df_with_indel_length_dist(df=loop_4bp_MH_SD_sims, 
                genomic_dist=d, shortest_indel=shortest_indel, indel_type='INS')

    
    
    
#     Exporting data
    nhej_df.loc[:,cols].to_csv(f'{out_pth}/nhej_{f}_with_genomic_indel_len_dist_no_1bp_deletion.vcf', sep='\t', index=False)
    nhej_df_all_indel_lengths.loc[:,cols].to_csv(f'{out_pth}/nhej_{f}_with_genomic_indel_len_dist.vcf', sep='\t', index=False)
    mmej_df.loc[:,cols].to_csv(f'{out_pth}/mmej_{f}_with_genomic_indel_len_dist.vcf', sep='\t', index=False)
    mmej_MH1_df.loc[:,cols].to_csv(f'{out_pth}/mmej_MH1_{f}_with_genomic_indel_len_dist.vcf', sep='\t', index=False)
    mmej_MH4_df.loc[:,cols].to_csv(f'{out_pth}/mmej_MH4_{f}_with_genomic_indel_len_dist.vcf', sep='\t', index=False)
    ins_df.loc[:,cols].to_csv(f'{out_pth}/ins_{f}_with_genomic_indel_len_dist.vcf', sep='\t', index=False)
    snap_df.loc[:,cols].to_csv(f'{out_pth}/snap_{f}_with_genomic_indel_len_dist.vcf', sep='\t', index=False)
    loop_df.loc[:,cols].to_csv(f'{out_pth}/loop_{f}_with_genomic_indel_len_dist.vcf', sep='\t', index=False)
    shortest_indel = -2

#     snap_1bp_MH_SD_df.loc[:,cols].to_csv(f'{out_pth}/snap_1bp_MH_SD_{f}_with_genomic_indel_len_dist.vcf', sep='\t', index=False)
#     snap_4bp_MH_SD_df.loc[:,cols].to_csv(f'{out_pth}/snap_4bp_MH_SD_{f}_with_genomic_indel_len_dist.vcf', sep='\t', index=False)
#     loop_1bp_MH_SD_df.loc[:,cols].to_csv(f'{out_pth}/loop_1bp_MH_SD_{f}_with_genomic_indel_len_dist.vcf', sep='\t', index=False)
#     loop_4bp_MH_SD_df.loc[:,cols].to_csv(f'{out_pth}/loop_4bp_MH_SD_{f}_with_genomic_indel_len_dist.vcf', sep='\t', index=False)



