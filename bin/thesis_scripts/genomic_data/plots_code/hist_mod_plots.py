

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



def calc_p_val(vec_a: pd.Series, vec_b: pd.Series,
            bonferroni_correction: int):
    """
    A function that calculates the p-value over
    the fake genomes.
    formula:
    # fake genomes of vec_a > # fake genomes of vec_b / # fake genomes
    """
    p_val = (pd.Series(vec_a[vec_a>vec_b]).shape[0]/vec_a.shape[0])*bonferroni_correction
    if p_val < 0.01: p,a = '0.01>p-value', '***'
    elif p_val < 0.05: p,a = '0.05>p-value', '**'
    elif p_val < 0.1: p,a = '0.1>p-value','*'
    else: p,a = 'not significant', ''
    return (p_val, p, a)
fs = 12

        
categories = ['H3K27Ac_chip_seq', 
    '8wg16_ChIP_chip', 'CBP_ChIP_chip', 'CBP_chip_seq', 
    'H3K27Ac_chip_seq', 'H3K4Me1_chip_seq', 'H3K4Me3_chip_seq', 'H3K9AC', 'H3K9Ac_chip_seq']
cn = [ 
    '8wg16 (ChIP chip)', 'CBP (ChIP chip)', 'CBP (chip seq)', 
    'H3K27Ac (chip seq)', 'H3K4Me1 (chip seq)', 'H3K4Me3 (chip seq)', 'H3K9AC', 'H3K9Ac (chip seq)']


datasets = ['Clark_AG_et_al_G3_2015', 'DGRP'] # 

datasets_pth_prefix=['Clark_GDL', 'DGRP']# 
for d,ds in zip(datasets,datasets_pth_prefix): 
    path_to_output = f'/home/labs/alevy/guyta/guy_master_project/results/Drosophila/{d}/plots'
    date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
    outpth_parent = f"{path_to_output}/{date}_plots"
    output = f"hist_mot_plots"
    dir = os.path.join(f"{outpth_parent}/{output}")
    if not os.path.exists(dir):	
        if not os.path.exists(os.path.join(f"{outpth_parent}")):
            os.mkdir(os.path.join(f"{outpth_parent}"))
        os.mkdir(os.path.join(dir))
    path_to_output = dir
    print(path_to_output)
    pth = f'/home/labs/alevy/guyta/guy_master_project/results/Drosophila/{d}/EMmej_results/20221126/{ds}_EMmej_formated_non_repeated_regions_AdultFemale_no_histon_modifications/20221201_blockbootstrap_output/EMmej_overall_estimations.tsv'
    df = pd.read_csv(pth, sep='\t')
    df['cat'] = 'No histone modifications'
    for cat,name in zip(categories, cn):
        pth = f'/home/labs/alevy/guyta/guy_master_project/results/Drosophila/{d}/EMmej_results/20221126/{ds}_EMmej_formated_non_repeated_regions_{cat}/20221201_blockbootstrap_output/EMmej_overall_estimations.tsv'
        tmp_df = pd.read_csv(pth, sep='\t')
        tmp_df['cat'] = name
        df = pd.concat([df,tmp_df])
        # cat = 'H3K27Ac_chip_seq'
        
    plt_df = pd.melt(df, id_vars=['cat','boot_n'], value_vars=['MMEJ',
                'NHEJ', 'snap', 'loop', 'SDMMEJ','unclassified_ins'])
    print(plt_df)
    mechanisms = ['MMEJ',
        'snap', 'loop', 'SDMMEJ']
    mechanisms_name=['MMEJ',
        'SD-Snap-Back', 'SD-Loop-Out', 'Snap-Back and Loop-Out']

    color_dict = {'MMEJ':'#fb8072', 'NHEJ':'#8dd3c7', 
        'snap':'#fdb462','loop':'#b3de69', 'unclassified_ins': '#bc80bd'}

    ann_df = pd.DataFrame(columns=['mechanism','p_val','p','ann', 'group1', 'group2'])
    for mech,mech_name in zip(mechanisms,mechanisms_name):
        ax = plt.figure(figsize = (7,5), dpi=80) # , layout='tight'
        ax = sns.boxplot(data=plt_df.loc[plt_df['variable'].isin([mech])],
                        y='value', x='variable', 
                        linewidth=1, hue='cat')  # , order= mechanisms, palette=color_dict
        ax.set_ylabel(f'Proportion of {mech_name}', fontsize=fs), ax.set_xlabel('')
        sns.move_legend(
            ax, loc="center left", ncol=1,
            bbox_to_anchor=(1.04, 0.5), 
            title=None, frameon=False,
        )
        plt.ylim(bottom=0, top=min(plt_df.loc[plt_df['variable'].isin([mech]), 'value'].max()*1.1, 1.05))
        plt.show()
        for cat in cn:
            if cat != 'No histone modifications':
                p_val,p,ann = calc_p_val(vec_b=df.loc[(df['cat'] == 'No histone modifications'), mech],
                vec_a=df.loc[(df['cat'] == cat), mech], bonferroni_correction=len(cn))

                tmp_ann_df = pd.DataFrame(
                            columns=['mechanism','p_val','p','ann', 'group1', 'group2'])
                tmp_ann_df.loc[0, ['p_val','p','ann']] = p_val,p,ann
                tmp_ann_df.loc[0, ['mechanism','group1', 'group2']] = mech, 'No histone modifications', cat
                # print(tmp_ann_df)
                ann_df = pd.concat([ann_df, tmp_ann_df])
        fig = ax.get_figure()
        fig.savefig(f"{path_to_output}/{ds}_hist_mod_{mech}_box.svg", bbox_inches='tight')
    ann_df.reset_index(inplace=True, drop=True)
    ann_df.to_csv(f"{path_to_output}/{ds}_stat_df.tsv", sep='\t')