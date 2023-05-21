"""
This script generates the plots for this dataset,
takes only summerized EMmej outputs as its input
"""

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
import scipy
import scipy.stats as stats
from scipy.stats import wilcoxon, pearsonr

sns.set_style("whitegrid")

pd.options.display.max_colwidth = 3500
pd.set_option('display.max_rows', 1000)
pd.set_option("display.max_columns", None)
boot_n = 100


# Construct an argument parser
all_args = argparse.ArgumentParser()
# Add arguments to the parser
all_args.add_argument("-cat1", "--cat1", required=True,
   help="name of first category (string)")
all_args.add_argument("-cat1p", "--catpath1", required=True,
   help="path to input data for first category (string)")
all_args.add_argument("-cat2", "--cat2", required=True,
   help="name of first category (string)")
all_args.add_argument("-cat2p", "--catpath2", required=True,
   help="path to input data for second category (string)")

all_args.add_argument("-o", "--outputpath", required=True,
      help="path to output file (string)")
all_args.add_argument("-c", "--comment", required=False, default='',
      help="comment to add to location (string)")

args = vars(all_args.parse_args())

path_to_output = args['outputpath']

date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
outpth_parent = f"{path_to_output}/{date}_plots"
output = f"{args['cat1']}_vs_{args['cat2']}_{args['comment']}"
dir = os.path.join(f"{outpth_parent}/{output}")
if not os.path.exists(dir):	
    if not os.path.exists(os.path.join(f"{outpth_parent}")):
        os.mkdir(os.path.join(f"{outpth_parent}"))
    os.mkdir(os.path.join(dir))
path_to_output = dir

def calc_p_val(vec_a: pd.Series, vec_b: pd.Series):
    """
    A function that calculates the p-value over
    the fake genomes.
    formula:
    # fake genomes of vec_a > # fake genomes of vec_b / # fake genomes
    """
    p_val = pd.Series(vec_a[vec_a>vec_b]).shape[0]/vec_a.shape[0]
    if p_val < 0.01: p,a = '0.01>p-value', '***'
    elif p_val < 0.05: p,a = '0.05>p-value', '**'
    elif p_val < 0.1: p,a = '0.1>p-value','*'
    else: p,a = 'not significant', ''
    return (p_val, p, a)

"""
Looking at overall astimetion over 'fake' genomes
"""
path_to_df_cat1 = f"{args['catpath1']}/EMmej_overall_estimations.tsv"
cat1_df = pd.read_csv(path_to_df_cat1, sep='\t')
cat1_df.dropna(inplace=True)
cat1_df['cat'] =args['cat1']

path_to_df_cat2 = f"{args['catpath2']}/EMmej_overall_estimations.tsv"
cat2_df = pd.read_csv(path_to_df_cat2, sep='\t')
cat2_df.dropna(inplace=True)
cat2_df['cat'] =args['cat2']

README = f"""
This folder contains all plots that has been generated using:
/home/labs/alevy/guyta/guy_master_project/scripts/EMmej/EMmej_util/EMmej_blockbootstrap/EMmej_results_plotting_general.py
using the following data:
{path_to_df_cat1}
{path_to_df_cat2}
"""
with open(f"{path_to_output}/README.txt", 'w') as f:
        f.write(README)

df = pd.concat([cat1_df, cat2_df])

plt_df = pd.melt(df, id_vars=['cat','boot_n'], value_vars=['MMEJ',
            'NHEJ', 'snap', 'loop', 'SDMMEJ','unclassified_ins'])
print(plt_df)
mechanisms = ['MMEJ',
    'NHEJ', 'snap', 'loop', 'SDMMEJ','unclassified_ins']

color_dict = {'MMEJ':'#fb8072', 'NHEJ':'#8dd3c7', 
    'snap':'#fdb462','loop':'#b3de69', 'unclassified_ins': '#bc80bd',
    args['cat1']: '#ef8a62', args['cat2']: '#67a9cf'}


"""
Deletions
"""
ax = plt.figure(figsize = (7,5), dpi=80) # , layout='tight'
ax = sns.boxplot(data=plt_df.loc[plt_df['variable'].isin(['MMEJ', 'NHEJ'])],
                 y='value', x='variable', 
                 linewidth=1, hue='cat', palette=color_dict)  # , order= mechanisms
ax.set_title('EM proportion estimation for MMEJ and NHEJ') ,ax.set_ylabel('Proportions'), ax.set_xlabel('')

MMEJ_theta_wilco = wilcoxon(df.loc[(df['cat'] == args['cat1']), 'MMEJ'], 
    df.loc[(df['cat'] == args['cat2']), 'MMEJ'], correction=False)

p_val,p,ann = calc_p_val(vec_a=df.loc[(df['cat'] == args['cat1']), 'MMEJ'],
     vec_b=df.loc[(df['cat'] == args['cat2']), 'MMEJ'])

ann_df = pd.DataFrame(index=['MMEJ', 'NHEJ'], 
                columns=['p_val','p','ann'])
# for m in ['MMEJ', 'NHEJ']:
for m in mechanisms:
    ann_df.loc[m,['p_val','p','ann']] = calc_p_val(
        vec_a=df.loc[(df['cat'] == args['cat1']), m],
        vec_b=df.loc[(df['cat'] == args['cat2']), m])

# if one wants to add * to the plot to indicate significance:
# for idx,i in enumerate(['MMEJ', 'NHEJ']):
#     ax.annotate(ann_df.loc[i,'ann'], xy=(idx, 0.65))

NHEJ_theta_wilco = wilcoxon(df.loc[(df['cat'] == args['cat1']), 'NHEJ'], 
    df.loc[(df['cat'] == args['cat2']), 'NHEJ'], correction=False)
text = f"""
Block size: 2Mb
Bootstrap N: {boot_n}
Hypothesis: 
Higher mechanism proportions in 
{args['cat1']} then {args['cat2']}
MMEJ:
    P-valus: {ann_df.loc['MMEJ', 'p_val']}
NHEJ:
    P-valus: {ann_df.loc['NHEJ', 'p_val']}
"""
# ax.text(1.55,.35,text)
sns.move_legend(
    ax, loc="center left", ncol=1,
    bbox_to_anchor=(1.04, 0.95), 
    title=None, frameon=False,
)
plt.ylim(bottom=0, top=min(plt_df.loc[plt_df['variable'].isin(['MMEJ', 'NHEJ']), 'value'].max()*1.1, 1.05))
ax.text(1.55,plt_df.loc[plt_df['variable'].isin(['MMEJ', 'NHEJ']), 'value'].mean()*0.65,text)

ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation=45)
fig = ax.get_figure()
fig.savefig(f"{path_to_output}/{args['cat1']}_vs_{args['cat2']}_box_deletions.svg", bbox_inches='tight')


"""
Deletions MMEJ only
"""
ax = plt.figure(figsize = (4,5), dpi=80)
ax = sns.boxplot(data=plt_df.loc[plt_df['variable'].isin(['MMEJ'])],
                 y='value', x='variable', 
                 linewidth=1, hue='cat', palette=color_dict)  # , order= mechanisms
ax.set_title('EM proportion estimation for MMEJ') ,ax.set_ylabel('Proportions'), ax.set_xlabel('')
sns.move_legend(
    ax, loc="center left", ncol=1,
    bbox_to_anchor=(1.04, 0.95), 
    title=None, frameon=False,
)

MMEJ_theta_wilco = wilcoxon(df.loc[(df['cat'] == args['cat1']), 'MMEJ'], 
    df.loc[(df['cat'] == args['cat2']), 'MMEJ'], correction=False)

p_val,p,ann = calc_p_val(vec_a=df.loc[(df['cat'] == args['cat1']), 'MMEJ'],
     vec_b=df.loc[(df['cat'] == args['cat2']), 'MMEJ'])

ann_df = pd.DataFrame(index=['MMEJ'], 
                columns=['p_val','p','ann'])
# for m in ['MMEJ', 'NHEJ']:
for m in mechanisms:
    ann_df.loc[m,['p_val','p','ann']] = calc_p_val(
        vec_a=df.loc[(df['cat'] == args['cat1']), m],
        vec_b=df.loc[(df['cat'] == args['cat2']), m])

# if one wants to add * to the plot to indicate significance:
# for idx,i in enumerate(['MMEJ', 'NHEJ']):
#     ax.annotate(ann_df.loc[i,'ann'], xy=(idx, 0.65))

text = f"""
Block size: 2Mb
Bootstrap N: {boot_n}
Hypothesis: 
Higher mechanism proportions in 
{args['cat1']} then {args['cat2']}
MMEJ:
    P-valus: {ann_df.loc['MMEJ', 'p_val']}

"""
plt.ylim(bottom=0, top=plt_df.loc[plt_df['variable'].isin(['MMEJ']), 'value'].mean()*1.5)
ax.text(0.55,plt_df.loc[plt_df['variable'].isin(['MMEJ']), 'value'].mean()*0.65,text)

ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation=45)
fig = ax.get_figure()
fig.savefig(f"{path_to_output}/{args['cat1']}_vs_{args['cat2']}_box_MMEJ_only.svg", bbox_inches='tight')


"""
Insertions
"""
ax = plt.figure(figsize = (7,5), dpi=80)
ax = sns.boxplot(data=plt_df.loc[plt_df['variable'].isin(['snap', 'loop', 'SDMMEJ','unclassified_ins'])],
                 y='value', x='variable', 
                 linewidth=1, hue='cat', palette=color_dict)  
ax.set_title('Proportion estimation of insertions') ,ax.set_ylabel(f"Proportions"), ax.set_xlabel('')
ax.set_xticks(range(len(['snap', 'loop', 'SDMMEJ','unclassified_ins']))) # <--- set the ticks first
ax.set_xticklabels(['SD-Snap-Back', 'SD-Loop-Out', 'Snap-Back and Loop-Out', 'Unclassified insertions'])
# ax.set(ylim=(0, 0.7))
# sns.move_legend(
#     ax, "upper center",
#     bbox_to_anchor=(1.2, 1), ncol=1, title=None, frameon=False,
# )

snap_wilco = wilcoxon(df.loc[(df['cat'] == args['cat1']), 'snap'], 
    df.loc[(df['cat'] == args['cat2']), 'snap'], correction=False)

loop_wilco = wilcoxon(df.loc[(df['cat'] == args['cat1']), 'loop'], 
    df.loc[(df['cat'] == args['cat2']), 'loop'], correction=False)
unclassified_ins_wilco = wilcoxon(df.loc[(df['cat'] == args['cat1']), 'unclassified_ins'], 
    df.loc[(df['cat'] == args['cat2']), 'unclassified_ins'], correction=False)
text = f"""
Block size: 2Mb
Bootstrap N: {boot_n}
Hypothesis: 
Higher mechanism proportions 
in {args['cat1']} then {args['cat2']}
SD-Snap back:
    P-valus: {ann_df.loc['snap', 'p_val']}
SD-Loop out:
    P-valus: {ann_df.loc['loop', 'p_val']}
Snap-Back and Loop-Out:
    P-valus: {ann_df.loc['SDMMEJ', 'p_val']}
Unclassified insertions:
    P-valus: {ann_df.loc['unclassified_ins', 'p_val']}

"""
# ax.text(3.6,.12,text)
sns.move_legend(
    ax, loc="center left", ncol=1,
    bbox_to_anchor=(1.04, 0.95), 
    title=None, frameon=False,
)
plt.ylim(bottom=0, top=min(plt_df.loc[plt_df['variable'].isin(['snap', 'loop', 'SDMMEJ', 'unclassified_ins']), 'value'].max()*1.1, 1.05))
ax.text(3.6,plt_df.loc[plt_df['variable'].isin(['snap', 'loop', 'SDMMEJ', 'unclassified_ins']), 'value'].mean()*0.5,text)
ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation=45)
fig = ax.get_figure()
fig.savefig(f"{path_to_output}/{args['cat1']}_vs_{args['cat2']}_box_insertions.svg", bbox_inches='tight')

df.to_csv(f"{path_to_output}/{args['cat1']}_vs_{args['cat2']}_data_for_overall_proportions_plots.tsv",
    sep='\t', index=False)




"""
Insertions _no_unclassified
"""
ax = plt.figure(figsize = (7,5), dpi=80)
ax = sns.boxplot(data=plt_df.loc[plt_df['variable'].isin(['snap', 'loop', 'SDMMEJ'])],
                 y='value', x='variable', 
                 linewidth=1, hue='cat', palette=color_dict)  
ax.set_title('Proportion estimation of insertions') ,ax.set_ylabel(f"Proportions"), ax.set_xlabel('')
ax.set_xticks(range(len(['snap', 'loop', 'SDMMEJ']))) # <--- set the ticks first
ax.set_xticklabels(['SD-Snap-Back', 'SD-Loop-Out', 'Snap-Back and Loop-Out'])
# ax.set(ylim=(0, 0.7))
sns.move_legend(
    ax, "upper center",
    bbox_to_anchor=(1.15, 1), ncol=1, title=None, frameon=False,
)

snap_wilco = wilcoxon(df.loc[(df['cat'] == args['cat1']), 'snap'], 
    df.loc[(df['cat'] == args['cat2']), 'snap'], correction=False)

loop_wilco = wilcoxon(df.loc[(df['cat'] == args['cat1']), 'loop'], 
    df.loc[(df['cat'] == args['cat2']), 'loop'], correction=False)

text = f"""
Block size: 2Mb
Bootstrap N: {boot_n}
Hypothesis: 
Higher mechanism proportions 
in {args['cat1']} then {args['cat2']}
SD-Snap back:
    P-valus: {ann_df.loc['snap', 'p_val']}
SD-Loop out:
    P-valus: {ann_df.loc['loop', 'p_val']}
Snap-Back and Loop-Out:
    P-valus: {ann_df.loc['SDMMEJ', 'p_val']}

"""
# ax.text(2.65,plt_df.loc[plt_df['variable'].isin(['snap', 'loop', 'SDMMEJ']), 'value'].mean(),text)
sns.move_legend(
    ax, loc="center left", ncol=1,
    bbox_to_anchor=(1.04, 0.95), 
    title=None, frameon=False,
)
plt.ylim(bottom=0, top=min(plt_df.loc[plt_df['variable'].isin(['snap', 'loop', 'SDMMEJ']), 'value'].max()*1.1, 1.05))
ax.text(2.65,plt_df.loc[plt_df['variable'].isin(['snap', 'loop', 'SDMMEJ']), 'value'].mean()*0.65,text)

ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation=45)
fig = ax.get_figure()
fig.savefig(f"{path_to_output}/{args['cat1']}_vs_{args['cat2']}_box_insertions_no_unclassified.svg", bbox_inches='tight')

df.to_csv(f"{path_to_output}/{args['cat1']}_vs_{args['cat2']}_data_for_overall_proportions_plots_no_unclassified.tsv",
    sep='\t', index=False)


"""
repair mechanism proportions stratified by indel length bins
"""
path_to_ins_cat1_by_indel_length = f"{args['catpath1']}/EMmej_estimations_by_indel_length_bins.tsv"

cat1_df = pd.read_csv(path_to_ins_cat1_by_indel_length, sep='\t')
cat1_df['cat'] =args['cat1']
cat1_df.sort_values(by=['indel_len_bin'], inplace=True)


path_to_ins_cat2_by_indel_length = f"{args['catpath2']}/EMmej_estimations_by_indel_length_bins.tsv"
cat2_df = pd.read_csv(path_to_ins_cat2_by_indel_length, sep='\t')
cat2_df.sort_values(by=['indel_len_bin'], inplace=True)
cat2_df['cat'] =args['cat2']


df = pd.concat([cat1_df, cat2_df])
df.reset_index(inplace=True, drop=True)

plt_df = pd.melt(df, id_vars=['cat','boot_n', 'indel_len_bin'], value_vars=['MMEJ',
            'NHEJ', 'snap', 'loop', 'SDMMEJ','unclassified_ins'])

insertions_lables = [str(i) for i in range(1,11)]# + ['[11:inf]']
insertions_lables = [str(i) for i in range(2,11)]# + ['[11:inf]']
if '[11:inf]' in plt_df.loc[(plt_df['variable'].isin(['snap', 'loop', 'SDMMEJ','unclassified_ins'])),'indel_len_bin'].to_list():
    insertions_lables = insertions_lables + ['[11:inf]']

deletions_lables = [str((-1*i)) for i in range(1,11)]# + ['[-50:-inf]']
deletions_lables = [str((-1*i)) for i in range(2,11)]# + ['[-50:-inf]']
print(plt_df.loc[(plt_df['variable'].isin(['MMEJ', 'NHEJ'])),'indel_len_bin'].value_counts())
if '[inf:-11]' in plt_df.loc[(plt_df['variable'].isin(['MMEJ', 'NHEJ'])),'indel_len_bin'].to_list():
    deletions_lables = deletions_lables + ['[inf:-11]']

plt_df.to_csv(f"{path_to_output}/plt_df_by_indel_len_bin.tsv", 
    sep='\t', index=False)

# ------------------------------------------------------
ax = plt.figure(figsize = (7,5), dpi=80)
ax = sns.boxplot(data=plt_df.loc[(plt_df['variable'].isin(['MMEJ']) 
                & (plt_df['indel_len_bin'].isin(deletions_lables))), :],
                 y='value', x='indel_len_bin', order=deletions_lables, 
                 linewidth=1, hue='cat', palette=color_dict)  # , order= mechanisms
ax.set_title('MMEJ proportion estimetion (\u03F4) for different indel lengths') ,ax.set_ylabel('\u03F4 MMEJ'), ax.set_xlabel('Indel length')
sns.move_legend(
    ax, "upper center",
    bbox_to_anchor=(1.15, 1), ncol=1, title=None, frameon=False,
)

# MMEJ_theta_wilco = wilcoxon(df.loc[(df['cat'] == 'Exon'), 'MMEJ'], 
#     df.loc[(df['cat'] == 'Non_Exon'), 'MMEJ'], correction=False)

# NHEJ_theta_wilco = wilcoxon(df.loc[(df['cat'] == 'Exon'), 'NHEJ'], 
#     df.loc[(df['cat'] == 'Non_Exon'), 'NHEJ'], correction=False)
# text = f"""
# Block size: 2Mb
# Bootstrap N: {boot_n}
# Wilcoxon test results:
# MMEJ:
#     P-valus: {MMEJ_theta_wilco.pvalue}
#     Statistic: {MMEJ_theta_wilco.statistic}
# NHEJ:
#     P-valus: {NHEJ_theta_wilco.pvalue}
#     Statistic: {MMEJ_theta_wilco.statistic}
# Filters applied:
#     - No repeated regions
#     - Non Exon is also non
#       any other genic feature
# """
# ax.text(1.55,.35,text)
ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation=45)
fig = ax.get_figure()
fig.savefig(f"{path_to_output}/{args['cat1']}_vs_{args['cat2']}_box_deletions_by_indel_length_bins.svg", bbox_inches='tight')




# -------- Barplot --------------
# look at guy_master_project/scripts/soybean/Liu.et.al_2020/data_exploration/0845_20220411_merged_chr_df_analysis.ipynb
def ci_calc(vec: pd.Series, p=0.05):
    avg = vec.mean()
    ci_l = vec.quantile(q=(p/2))
    ci_h = vec.quantile(q=1-(p/2))
    return (avg-ci_l, ci_h-avg)
    

"""
Deletions
"""
barWidth = 0.3
linewidth=0.8
elinewidth=0.8
error_kw=dict(lw=1, capsize=3, capthick=1)
mech = ['MMEJ', 'NHEJ']
for m in mech:
    # Creating a dataframe for the plot:
    del_cat1 = plt_df.loc[(plt_df['variable'].isin([m]) &
        (plt_df['indel_len_bin'].isin(deletions_lables)) &
        (plt_df['cat'] == args['cat1'])), ['value','indel_len_bin']]

    del_cat2 = plt_df.loc[(plt_df['variable'].isin([m]) &
        (plt_df['indel_len_bin'].isin(deletions_lables)) &
        (plt_df['cat'] == args['cat2'])), ['value','indel_len_bin']]
    del_cat1_ci = [ci_calc(del_cat1.loc[del_cat1['indel_len_bin'] == bin,'value']) for bin in deletions_lables]
    del_cat2_ci = [ci_calc(del_cat2.loc[del_cat2['indel_len_bin'] == bin,'value']) for bin in deletions_lables]
    
    p_df = pd.DataFrame(columns=['indel_len_bin', 'avg'])
    p_df['indel_len_bin'] = deletions_lables
    p_df['avg'] = [del_cat1.loc[del_cat1['indel_len_bin'] == bin, 'value'].mean()
        for bin in deletions_lables]
    p_df['CI'] = del_cat1_ci
    p_df['cat'] = args['cat1']

    p_df2 = pd.DataFrame(columns=['indel_len_bin', 'avg'])
    p_df2['indel_len_bin'] = deletions_lables
    p_df2['avg'] = [del_cat2.loc[del_cat2['indel_len_bin'] == bin, 'value'].mean()
        for bin in deletions_lables]
    p_df2['CI'] = del_cat2_ci 
    p_df2['cat'] = args['cat2']
    p_df = pd.concat([p_df, p_df2]).reset_index(drop=True)
    
    # setting the x axis positions
    r1 = np.arange(len(deletions_lables))
    r2 = [x + barWidth for x in r1]

    # Creating the bars for repeats
    # plt.style.use('seaborn-darkgrid')
    plt.style.use('seaborn-white')
    ax = plt.figure(figsize = (7,5), dpi=80)
    
    plt.bar(r1, p_df.loc[p_df['cat'] == args['cat1'], 'avg'], width=barWidth, color='#ef8a62', 
        yerr = [[i[0] for i in del_cat1_ci], [i[1] for i in del_cat1_ci]], 
            edgecolor='black',linewidth=linewidth ,
            error_kw=error_kw, label=args['cat1'])
    
    # Creating the bars for non repeats
    plt.bar(r2, p_df.loc[p_df['cat'] == args['cat2'], 'avg'], width=barWidth, color='#67a9cf', 
            yerr = [[i[0] for i in del_cat2_ci], [i[1] for i in del_cat2_ci]], 
            edgecolor='black',linewidth=linewidth ,
            error_kw=error_kw, label=args['cat2'])
    
    # Set plot parameters
    plt.xticks([r + barWidth for r in range(len(deletions_lables))], deletions_lables, rotation = 45)
    plt.legend(bbox_to_anchor=(1.01, 1.05))
    plt.ylabel(f'{m} Proportions')
    plt.xlabel('Deletion length (bp)')
    plt.title(f'{m} proportion estimetion for different indel lengths')
    fig = ax.get_figure()
    fig.savefig(f"{path_to_output}/{args['cat1']}_vs_{args['cat2']}_bar_{m}_by_indel_length_bins.svg", bbox_inches='tight')
    
    p_df.to_csv(f"{path_to_output}/{args['cat1']}_vs_{args['cat2']}_bar_{m}_by_indel_length_bins.tsv", 
        sep='\t', index=False)
"""
Insertions
"""
mech = ['snap', 'loop', 'unclassified_ins','SDMMEJ']
mech_name = ['SD-Snap-Back', 'SD-Loop-Out', 'Unclassified-insertion','Snap-Back and Loop-Out']
for m,n in zip(mech,mech_name):
    # Creating a dataframe for the plot:
    ins_cat1 = plt_df.loc[(plt_df['variable'].isin([m]) &
        (plt_df['indel_len_bin'].isin(insertions_lables)) &
        (plt_df['cat'] == args['cat1'])), ['value','indel_len_bin']]

    ins_cat2 = plt_df.loc[(plt_df['variable'].isin([m]) &
        (plt_df['indel_len_bin'].isin(insertions_lables)) &
        (plt_df['cat'] == args['cat2'])), ['value','indel_len_bin']]
    ins_cat1_ci = [ci_calc(ins_cat1.loc[ins_cat1['indel_len_bin'] == bin,'value']) for bin in insertions_lables] 
    ins_cat2_ci = [ci_calc(ins_cat2.loc[ins_cat2['indel_len_bin'] == bin,'value']) for bin in insertions_lables]

    p_df = pd.DataFrame(columns=['indel_len_bin', 'avg'])
    p_df['indel_len_bin'] = insertions_lables
    p_df['avg'] = [ins_cat1.loc[ins_cat1['indel_len_bin'] == bin, 'value'].mean()
             for bin in insertions_lables]

    p_df['CI'] = ins_cat1_ci
    p_df['cat'] = args['cat1']

    p_df2 = pd.DataFrame(columns=['indel_len_bin', 'avg'])
    p_df2['indel_len_bin'] = insertions_lables
    p_df2['avg'] = [ins_cat2.loc[ins_cat2['indel_len_bin'] == bin, 'value'].mean()
                        for bin in p_df['indel_len_bin'].unique()]

    p_df2['CI'] = ins_cat2_ci 
    p_df2['cat'] = args['cat2']
    p_df = pd.concat([p_df, p_df2]).reset_index(drop=True)
    
    # setting the x axis positions
    r1 = np.arange(len(insertions_lables))
    r2 = [x + barWidth for x in r1]
    # Creating the bars for repeats
    # plt.style.use('seaborn-darkgrid')
    plt.style.use('seaborn-white')
    
    ax = plt.figure(figsize = (7,5), dpi=80)
    plt.bar(r1, p_df.loc[p_df['cat'] == args['cat1'], 'avg'], width=barWidth, color='#ef8a62', 
        yerr = [[i[0] for i in ins_cat1_ci], [i[1] for i in ins_cat1_ci]], 
            edgecolor='black',linewidth=linewidth ,
            error_kw=error_kw, label=args['cat1'])
    
    # Creating the bars for non repeats
    plt.bar(r2, p_df.loc[p_df['cat'] == args['cat2'], 'avg'], width=barWidth, color='#67a9cf', 
            yerr = [[i[0] for i in ins_cat2_ci], [i[1] for i in ins_cat2_ci]], 
            edgecolor='black',linewidth=linewidth ,
            error_kw=error_kw, label=args['cat2'])
    
    # Set plot parameters
    plt.xticks([r + barWidth for r in range(len(insertions_lables))], insertions_lables, rotation = 45)
    plt.legend(bbox_to_anchor=(1.25, 1.01))
    plt.ylabel(f'{n} Proportions')
    plt.xlabel('Insertion length (bp)')
    plt.title(f'{n} proportion estimetion for different indel lengths')
    fig = ax.get_figure()
    fig.savefig(f"{path_to_output}/{args['cat1']}_vs_{args['cat2']}_bar_{n}_by_indel_length_bins.svg", bbox_inches='tight')
    # fig.savefig(f"{path_to_output}/{args['cat1']}_vs_{args['cat2']}_bar_{n}_by_indel_length_bins.png", bbox_inches='tight')
    p_df.to_csv(
        f"{path_to_output}/{args['cat1']}_vs_{args['cat2']}_bar_{n}_by_indel_length_bins.tsv", 
        sep='\t', index=False)

# print(plt.style.available)


# """
# Ploting the proportions of repair mechanisms along the chromosomes
# """
# # laoding the data
# # Creating a list of all folders that contains blockbootstraping runs
# exons_blocks = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/GDL_indel_EMmej_format_non_repeated_regions_exon_no_MT/20221021'
# dir_ls = glob.glob(f"{exons_blocks}/*")
# cols=['CHR', 'blockN',
#         'POS', 'original_pos','variant_id', 'direction','ancestral_indel', 'derived_indel', 
#         'del_mmej_lk', 'del_mmej_p_val', 
#         'snap','SD_snap_back_p_val',
#         'loop', 'SD_loop_out_p_val', 
#         'NHEJ_lk','NHEJ_p_val', 
#         'unclassified_ins', 'unclassified_ins_p_val',
#         'indel_len','del_mmej_cand','del_mmej_cand_len']
# EMmej_MM_outputs = pd.DataFrame(columns=cols)

# # Empty list to get names of runs that havn't finished yet
# undone_runs = []

# for i, fol in enumerate(dir_ls):
#       output_ls = glob.glob(f"{fol}/*")    
#       """
#       Markov model merging
#       """
#       MM_outputs = [f for f in os.listdir(f"{fol}/") if re.match(r'(.+?)_markov_output.tsv', f)]
#       if len(MM_outputs) > 0:
#             MM_output_df = pd.read_csv(f"{fol}/{MM_outputs[0]}", 
#                   sep='\t')
#             chr_name = re.search(r'chr_(.+?)_blockN', fol).group(1)
#             block_n = re.search(r'blockN(.+?)_', fol).group(1)
#             MM_output_df['blockN'] = block_n
#             MM_output_df['blockID'] = f"{chr_name}_{block_n}"
#             EMmej_MM_outputs = pd.concat([EMmej_MM_outputs, MM_output_df])
#       else: undone_runs.append(fol)

# boot_df = pd.DataFrame(columns=['CHR', 'blockN',
#         'snap', 'loop', 'unclassified_ins'])




# Command example:
# python3 EMmej_results_plotting_general.py -cat1 Repeat -cat1p /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/GDL_indel_EMmej_format_repeated_regions/20221026_blockbootstrap_output -cat2 Non_Repeat -cat2p /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/GDL_indel_EMmej_format_non_repeated_regions/20221026_blockbootstrap_output -o /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/plots/repeats_vs_non_repeats/20221026_blockbootstrap