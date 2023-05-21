import os
import glob
import re 
import argparse
import datetime

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.colors as mcolors
import matplotlib as mpl
# mpl.style.use('seaborn-deep')
import seaborn as sns


out_pth = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/plots/20221227_plots/exon_non_genic_indel_len_dist'
# Exon
pth='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/EMmej_results/20221126/DGRP_EMmej_formated_non_repeated_regions_exons/20221126_2000000bp_blocks_EMmej_res'
file_ls = glob.glob(f"{pth}/*/*_markov_output.tsv")
print(len(file_ls))
df_i = pd.read_csv(file_ls[0], sep='\t')

exon_df = pd.DataFrame(columns=df_i.columns.to_list())

for d in file_ls:
    exon_df = pd.concat([exon_df, pd.read_csv(d, sep='\t')])
    exon_df.reset_index(drop=True, inplace=True)

exon_df = exon_df.loc[(~exon_df['indel_len'].isin([1,-1])),:]

# non-genic
pth='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/EMmej_results/20221126/DGRP_EMmej_formated_non_repeated_regions_non_genic/20221126_2000000bp_blocks_EMmej_res'
file_ls = glob.glob(f"{pth}/*/*_markov_output.tsv")
print(len(file_ls))
df_i = pd.read_csv(file_ls[0], sep='\t')

non_genic_df = pd.DataFrame(columns=df_i.columns.to_list())

for d in file_ls:
    non_genic_df = pd.concat([non_genic_df, pd.read_csv(d, sep='\t')])
    non_genic_df.reset_index(drop=True, inplace=True)

non_genic_df = non_genic_df.loc[(~non_genic_df['indel_len'].isin([1,-1])),:]
# indel length dist of Loop out exons
snap_ins_cutoff = 0.005
loop_ins_cutoff = 0.003
exon_df['bonferroni_cutoff_snap'] = snap_ins_cutoff / exon_df['variant_id_N']
exon_df['bonferroni_cutoff_loop'] = loop_ins_cutoff / exon_df['variant_id_N']
non_genic_df['bonferroni_cutoff_snap'] = snap_ins_cutoff / non_genic_df['variant_id_N']
non_genic_df['bonferroni_cutoff_loop'] = loop_ins_cutoff / non_genic_df['variant_id_N']


ins_len_ls = [i for i in range(2, 115)]
del_len_ls = [i for i in range(-115, -1)]#[::-1]

# Loop-Out Exon
loop_indel_lengths = pd.DataFrame(columns=['indel_len', 'count'])
loop_indel_lengths['indel_len'] = ins_len_ls

for l in loop_indel_lengths['indel_len']:
    tmp_count_loop = exon_df.loc[((exon_df['indel_len'] == l)
    & (exon_df['SD_loop_out_p_val'] <= exon_df['bonferroni_cutoff_loop'])), 
                                    'variant_id'].unique().shape[0]
    loop_indel_lengths.loc[loop_indel_lengths['indel_len'] == l, 'count'] = tmp_count_loop

loop_indel_lengths['prop'] = loop_indel_lengths['count']/loop_indel_lengths['count'].sum()
loop_indel_lengths['cat'] = 'Exon'

# indel length dist of Loop out non_genic
loop_indel_lengths_non_genic = pd.DataFrame(columns=['indel_len', 'count'])
loop_indel_lengths_non_genic['indel_len'] = ins_len_ls

for l in loop_indel_lengths_non_genic['indel_len']:
    tmp_count_loop = non_genic_df.loc[((non_genic_df['indel_len'] == l)
    & (non_genic_df['SD_loop_out_p_val'] <= non_genic_df['bonferroni_cutoff_loop'])), 
                                    'variant_id'].unique().shape[0]
    loop_indel_lengths_non_genic.loc[loop_indel_lengths_non_genic['indel_len'] == l, 'count'] = tmp_count_loop

loop_indel_lengths_non_genic['prop'] = loop_indel_lengths_non_genic['count']/loop_indel_lengths_non_genic['count'].sum()
loop_indel_lengths_non_genic['cat'] = 'Non_genic'
loop_indel_lengths = pd.concat([loop_indel_lengths, loop_indel_lengths_non_genic])
loop_indel_lengths.reset_index(drop=True, inplace=True)
loop_indel_lengths.to_csv(f'{out_pth}/exon_non_genic_indel_length_loop_out.tsv', sep='\t', index=False)

# Insertions only
indel_lengths_dist = pd.DataFrame(columns=['indel_len', 'count'])
indel_lengths_dist['indel_len'] = ins_len_ls

for l in indel_lengths_dist['indel_len']:
    tmp_count_loop = exon_df.loc[(exon_df['indel_len'] == l), 
                                    'variant_id'].unique().shape[0]
    indel_lengths_dist.loc[indel_lengths_dist['indel_len'] == l, 'count'] = tmp_count_loop

indel_lengths_dist['prop'] = indel_lengths_dist['count']/indel_lengths_dist['count'].sum()
indel_lengths_dist['cat'] = 'Exon'
# general indel length dist non genic
indel_lengths_dist_non_genic = pd.DataFrame(columns=['indel_len', 'count'])
indel_lengths_dist_non_genic['indel_len'] = ins_len_ls

for l in indel_lengths_dist_non_genic['indel_len']:
    tmp_count_loop = non_genic_df.loc[(non_genic_df['indel_len'] == l), 
                                    'variant_id'].unique().shape[0]
    indel_lengths_dist_non_genic.loc[indel_lengths_dist_non_genic['indel_len'] == l, 'count'] = tmp_count_loop

indel_lengths_dist_non_genic['prop'] = indel_lengths_dist_non_genic['count']/indel_lengths_dist_non_genic['count'].sum()
indel_lengths_dist_non_genic['cat'] = 'Non_genic'
# print(indel_lengths_dist_non_genic)
indel_len_dist = pd.concat([indel_lengths_dist, indel_lengths_dist_non_genic])
indel_len_dist.reset_index(drop=True, inplace=True)
indel_len_dist.to_csv(f'{out_pth}/insertions_exon_non_genic_indel_length.tsv', sep='\t', index=False)

# Deletions only
indel_lengths_dist = pd.DataFrame(columns=['indel_len', 'count'])
indel_lengths_dist['indel_len'] = del_len_ls[::-1]

for l in indel_lengths_dist['indel_len']:
    tmp_count_loop = exon_df.loc[(exon_df['indel_len'] == l), 
                                    'variant_id'].unique().shape[0]
    indel_lengths_dist.loc[indel_lengths_dist['indel_len'] == l, 'count'] = tmp_count_loop

indel_lengths_dist['prop'] = indel_lengths_dist['count']/indel_lengths_dist['count'].sum()
indel_lengths_dist['cat'] = 'Exon'
# general indel length dist non genic
indel_lengths_dist_non_genic = pd.DataFrame(columns=['indel_len', 'count'])
indel_lengths_dist_non_genic['indel_len'] = del_len_ls[::-1]

for l in indel_lengths_dist_non_genic['indel_len']:
    tmp_count_loop = non_genic_df.loc[(non_genic_df['indel_len'] == l), 
                                    'variant_id'].unique().shape[0]
    indel_lengths_dist_non_genic.loc[indel_lengths_dist_non_genic['indel_len'] == l, 'count'] = tmp_count_loop

indel_lengths_dist_non_genic['prop'] = indel_lengths_dist_non_genic['count']/indel_lengths_dist_non_genic['count'].sum()
indel_lengths_dist_non_genic['cat'] = 'Non_genic'
# print(indel_lengths_dist_non_genic)
indel_len_dist = pd.concat([indel_lengths_dist, indel_lengths_dist_non_genic])
indel_len_dist.reset_index(drop=True, inplace=True)
indel_len_dist.to_csv(f'{out_pth}/deletions_exon_non_genic_indel_length.tsv', sep='\t', index=False)



# Deletions MMEJ only
indel_lengths_dist = pd.DataFrame(columns=['indel_len', 'count'])
indel_lengths_dist['indel_len'] = del_len_ls[::-1]

for l in indel_lengths_dist['indel_len']:
    tmp_count_loop = exon_df.loc[((exon_df['indel_len'] == l) & (~exon_df['del_mmej_cand'].isna())), 
                                    'variant_id'].unique().shape[0]
    indel_lengths_dist.loc[indel_lengths_dist['indel_len'] == l, 'count'] = tmp_count_loop

indel_lengths_dist['prop'] = indel_lengths_dist['count']/indel_lengths_dist['count'].sum()
indel_lengths_dist['cat'] = 'Exon'
# general indel length dist non genic
indel_lengths_dist_non_genic = pd.DataFrame(columns=['indel_len', 'count'])
indel_lengths_dist_non_genic['indel_len'] = del_len_ls[::-1]

for l in indel_lengths_dist_non_genic['indel_len']:
    tmp_count_loop = non_genic_df.loc[((non_genic_df['indel_len'] == l)& (~non_genic_df['del_mmej_cand'].isna())), 
                                    'variant_id'].unique().shape[0]
    indel_lengths_dist_non_genic.loc[indel_lengths_dist_non_genic['indel_len'] == l, 'count'] = tmp_count_loop

indel_lengths_dist_non_genic['prop'] = indel_lengths_dist_non_genic['count']/indel_lengths_dist_non_genic['count'].sum()
indel_lengths_dist_non_genic['cat'] = 'Non_genic'
# print(indel_lengths_dist_non_genic)
indel_len_dist = pd.concat([indel_lengths_dist, indel_lengths_dist_non_genic])
indel_len_dist.reset_index(drop=True, inplace=True)
indel_len_dist.to_csv(f'{out_pth}/deletions_MMEJ_exon_non_genic_indel_length.tsv', sep='\t', index=False)



# All indels
# general indel length dist exion
indel_lengths_dist = pd.DataFrame(columns=['indel_len', 'count'])
indel_lengths_dist['indel_len'] = del_len_ls + ins_len_ls

for l in indel_lengths_dist['indel_len']:
    tmp_count_loop = exon_df.loc[(exon_df['indel_len'] == l), 
                                    'variant_id'].unique().shape[0]
    indel_lengths_dist.loc[indel_lengths_dist['indel_len'] == l, 'count'] = tmp_count_loop

indel_lengths_dist['prop'] = indel_lengths_dist['count']/indel_lengths_dist['count'].sum()
indel_lengths_dist['cat'] = 'Exon'
# general indel length dist non genic
indel_lengths_dist_non_genic = pd.DataFrame(columns=['indel_len', 'count'])
indel_lengths_dist_non_genic['indel_len'] = del_len_ls + ins_len_ls

for l in indel_lengths_dist_non_genic['indel_len']:
    tmp_count_loop = non_genic_df.loc[(non_genic_df['indel_len'] == l), 
                                    'variant_id'].unique().shape[0]
    indel_lengths_dist_non_genic.loc[indel_lengths_dist_non_genic['indel_len'] == l, 'count'] = tmp_count_loop

indel_lengths_dist_non_genic['prop'] = indel_lengths_dist_non_genic['count']/indel_lengths_dist_non_genic['count'].sum()
indel_lengths_dist_non_genic['cat'] = 'Non_genic'
# print(indel_lengths_dist_non_genic)
indel_len_dist = pd.concat([indel_lengths_dist, indel_lengths_dist_non_genic])
indel_len_dist.reset_index(drop=True, inplace=True)
indel_len_dist.to_csv(f'{out_pth}/all_indels_exon_non_genic_indel_length.tsv', sep='\t', index=False)



"""
Ploting
"""
# Insertions
color_dict = {
    'Exon': '#ef8a62', 'Non-genic': '#67a9cf'}

pth = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/plots/20221227_plots/exon_non_genic_indel_len_dist/exon_non_genic_indel_length_loop_out.tsv'
loop_df = pd.read_csv(pth, sep='\t')
plt_df = pd.melt(loop_df, id_vars=['cat','indel_len'], value_vars=['prop'])
plt_df.loc[plt_df['cat'] == 'Non_genic', 'cat'] = 'Non-genic'
plt_df = plt_df.loc[plt_df['indel_len'] < 40, :]

ax = plt.figure(figsize = (8.5,5), dpi=80)
ax = sns.barplot(x = plt_df['indel_len'], y=plt_df['value'], hue=plt_df['cat'],palette=color_dict)

ax.set_ylabel('Density',fontsize=12), ax.set_xlabel('Insertion length',fontsize=12)
ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation=45)

fig = ax.get_figure()
fig.savefig(f"{out_pth}/Insertions_exon_non_genic_loop_indel_length_dist.svg", bbox_inches='tight')


pth = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/plots/20221227_plots/exon_non_genic_indel_len_dist/insertions_exon_non_genic_indel_length.tsv'
df = pd.read_csv(pth, sep='\t')
plt_df = pd.melt(df, id_vars=['cat','indel_len'], value_vars=['prop'])
plt_df.loc[plt_df['cat'] == 'Non_genic', 'cat'] = 'Non-genic'
plt_df = plt_df.loc[((plt_df['indel_len'] < 40)), :]
print(plt_df['indel_len'])
ax = plt.figure(figsize = (8.5,5), dpi=80)
ax = sns.barplot(x = plt_df['indel_len'], y=plt_df['value'], hue=plt_df['cat'], palette=color_dict)

ax.set_ylabel('Density',fontsize=12), ax.set_xlabel('Insertion length',fontsize=12)
ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation=45)

fig = ax.get_figure()
fig.savefig(f"{out_pth}/Insertions_Exon_non_genic_indel_length_dist.svg", bbox_inches='tight')


# Deletions only
pth = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/plots/20221227_plots/exon_non_genic_indel_len_dist/deletions_exon_non_genic_indel_length.tsv'
df = pd.read_csv(pth, sep='\t')
plt_df = pd.melt(df, id_vars=['cat','indel_len'], value_vars=['prop'])
plt_df.loc[plt_df['cat'] == 'Non_genic', 'cat'] = 'Non-genic'
plt_df = plt_df.loc[((plt_df['indel_len'] > -40)), :]
print(plt_df['indel_len'])
ax = plt.figure(figsize = (8.5,5), dpi=80)
ax = sns.barplot(x = plt_df['indel_len'], y=plt_df['value'], hue=plt_df['cat'], palette=color_dict)
ax.invert_xaxis()
ax.set_ylabel('Density',fontsize=12), ax.set_xlabel('Deletion length',fontsize=12)
ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation=45)

fig = ax.get_figure()
fig.savefig(f"{out_pth}/deletions_Exon_non_genic_indel_length_dist.svg", bbox_inches='tight')


# Deletions MMEJ only
pth = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/plots/20221227_plots/exon_non_genic_indel_len_dist/deletions_MMEJ_exon_non_genic_indel_length.tsv'
df = pd.read_csv(pth, sep='\t')
plt_df = pd.melt(df, id_vars=['cat','indel_len'], value_vars=['prop'])
plt_df.loc[plt_df['cat'] == 'Non_genic', 'cat'] = 'Non-genic'
plt_df = plt_df.loc[((plt_df['indel_len'] > -40)), :]
print(plt_df['indel_len'])
ax = plt.figure(figsize = (8.5,5), dpi=80)
ax = sns.barplot(x = plt_df['indel_len'], y=plt_df['value'], hue=plt_df['cat'], palette=color_dict)
ax.invert_xaxis()
ax.set_ylabel('Density',fontsize=12), ax.set_xlabel('Deletion length',fontsize=12)
ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation=45)

fig = ax.get_figure()
fig.savefig(f"{out_pth}/deletions_MMEJ_Exon_non_genic_indel_length_dist.svg", bbox_inches='tight')

# All indels
pth = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/plots/20221227_plots/exon_non_genic_indel_len_dist/all_indels_exon_non_genic_indel_length.tsv'
df = pd.read_csv(pth, sep='\t')
plt_df = pd.melt(df, id_vars=['cat','indel_len'], value_vars=['prop'])
plt_df.loc[plt_df['cat'] == 'Non_genic', 'cat'] = 'Non-genic'
plt_df = plt_df.loc[((plt_df['indel_len'] < 40) & (plt_df['indel_len'] > -40)), :]
print(plt_df['indel_len'])
ax = plt.figure(figsize = (16,5), dpi=80)
ax = sns.barplot(x = plt_df['indel_len'], y=plt_df['value'], hue=plt_df['cat'], palette=color_dict)

ax.set_ylabel('Density',fontsize=12), ax.set_xlabel('Indel length',fontsize=12)
ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation=45)

fig = ax.get_figure()
fig.savefig(f"{out_pth}/all_indels_Exon_non_genic_indel_length_dist.svg", bbox_inches='tight')