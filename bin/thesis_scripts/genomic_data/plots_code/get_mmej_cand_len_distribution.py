
import os
import glob
import re 
import argparse
import datetime

import pandas as pd
import numpy as np



pth='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/20221126/Clark_GDL_EMmej_formated_non_repeated_regions_genic/20221126_2000000bp_blocks_EMmej_res'
file_ls = glob.glob(f"{pth}/*/*_markov_output.tsv")
print(len(file_ls))
df_i = pd.read_csv(file_ls[0], sep='\t')

genic_df = pd.DataFrame(columns=df_i.columns.to_list())

for d in file_ls:
    genic_df = pd.concat([genic_df, pd.read_csv(d, sep='\t')])
    genic_df.reset_index(drop=True, inplace=True)

genic_df = genic_df.loc[(~genic_df['indel_len'].isin([1,-1])),:]

pth='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/20221126/Clark_GDL_EMmej_formated_non_repeated_regions_non_genic/20221126_2000000bp_blocks_EMmej_res'
file_ls = glob.glob(f"{pth}/*/*_markov_output.tsv")
print(len(file_ls))
df_i = pd.read_csv(file_ls[0], sep='\t')

non_genic_df = pd.DataFrame(columns=df_i.columns.to_list())

for d in file_ls:
    non_genic_df = pd.concat([non_genic_df, pd.read_csv(d, sep='\t')])
    non_genic_df.reset_index(drop=True, inplace=True)

non_genic_df = non_genic_df.loc[(~non_genic_df['indel_len'].isin([1,-1])),:]

whole_data = pd.concat([genic_df, non_genic_df])


whole_data['del_mmej_cand']
print(whole_data['del_mmej_cand_len'].value_counts()/whole_data['del_mmej_cand_len'].value_counts().sum())







pth='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/20221126/Clark_GDL_EMmej_formated_non_repeated_regions/20221126_2000000bp_blocks_EMmej_res'
file_ls = glob.glob(f"{pth}/*/*_markov_output.tsv")
print(len(file_ls))
df_i = pd.read_csv(file_ls[0], sep='\t')

non_repeated_regions = pd.DataFrame(columns=df_i.columns.to_list())

for d in file_ls:
    non_repeated_regions = pd.concat([non_repeated_regions, pd.read_csv(d, sep='\t')])
    genic_df.reset_index(drop=True, inplace=True)

non_repeated_regions = non_repeated_regions.loc[(~non_repeated_regions['indel_len'].isin([1,-1])),:]

pth='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/20221126/Clark_GDL_EMmej_formated_repeated_regions/20221126_2000000bp_blocks_EMmej_res'
file_ls = glob.glob(f"{pth}/*/*_markov_output.tsv")
print(len(file_ls))
df_i = pd.read_csv(file_ls[0], sep='\t')

repeated_regions = pd.DataFrame(columns=df_i.columns.to_list())

for d in file_ls:
    repeated_regions = pd.concat([repeated_regions, pd.read_csv(d, sep='\t')])
    repeated_regions.reset_index(drop=True, inplace=True)

repeated_regions = repeated_regions.loc[(~repeated_regions['indel_len'].isin([1,-1])),:]

whole_data = pd.concat([non_repeated_regions, repeated_regions])


# print('repeated')
# print(repeated_regions['del_mmej_cand_len'].value_counts()/repeated_regions['del_mmej_cand_len'].value_counts().sum())
# print('non repeated')
# print(non_repeated_regions['del_mmej_cand_len'].value_counts()/non_repeated_regions['del_mmej_cand_len'].value_counts().sum())
# print('Whole data')
# print(whole_data['del_mmej_cand_len'].value_counts()/whole_data['del_mmej_cand_len'].value_counts().sum())
wd = whole_data['del_mmej_cand_len'].value_counts()/whole_data.shape[0]
print(wd,wd.sum())
print(wd[wd.index>=4].sum(), wd[wd.index<4].sum())

print('DGRP2')



pth='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/EMmej_results/20221126/DGRP_EMmej_formated_non_repeated_regions/20221126_2000000bp_blocks_EMmej_res'
file_ls = glob.glob(f"{pth}/*/*_markov_output.tsv")
print(len(file_ls))
df_i = pd.read_csv(file_ls[0], sep='\t')

non_repeated_regions = pd.DataFrame(columns=df_i.columns.to_list())

for d in file_ls:
    non_repeated_regions = pd.concat([non_repeated_regions, pd.read_csv(d, sep='\t')])
    genic_df.reset_index(drop=True, inplace=True)

non_repeated_regions = non_repeated_regions.loc[(~non_repeated_regions['indel_len'].isin([1,-1])),:]

pth='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/EMmej_results/20221126/DGRP_EMmej_formated_repeated_regions/20221126_2000000bp_blocks_EMmej_res'
file_ls = glob.glob(f"{pth}/*/*_markov_output.tsv")
print(len(file_ls))
df_i = pd.read_csv(file_ls[0], sep='\t')

repeated_regions = pd.DataFrame(columns=df_i.columns.to_list())

for d in file_ls:
    repeated_regions = pd.concat([repeated_regions, pd.read_csv(d, sep='\t')])
    repeated_regions.reset_index(drop=True, inplace=True)

repeated_regions = repeated_regions.loc[(~repeated_regions['indel_len'].isin([1,-1])),:]

whole_data = pd.concat([non_repeated_regions, repeated_regions])

# print('repeated')
# print(repeated_regions['del_mmej_cand_len'].value_counts()/repeated_regions.shape[0], )
# print('non repeated')
# print(non_repeated_regions['del_mmej_cand_len'].value_counts()/non_repeated_regions.shape[0])
print('Whole data')
wd = whole_data['del_mmej_cand_len'].value_counts()/whole_data.shape[0]
print(wd,wd.sum())
print(wd[wd.index>=4].sum(), wd[wd.index<4].sum())
