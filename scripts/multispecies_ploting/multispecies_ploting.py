import os

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import seaborn as sns

out_path = '/home/labs/alevy/guyta/guy_master_project/results/multispecies_ploting'
pth = '/home/labs/alevy/guyta/guy_master_project/results/multispecies_data_to_plot/non_repeated_region_mechanisms_estimation'
datasets = ['Clark', 'maize', 'rice', 'A_thaliana','tomato_RGTR']
org = ['A.thaliana', 'D.melanogaster', 'Maize', 'Rice', 'Tomato']
org = ['Drosophila\n(D.melanogaster)', 'Maize\n(Zea mays)', 'Rice\n(Oryza sativa L.)','Arabidopsis\n(A.thaliana)', 'Tomato\n']
org = ['Drosophila', 'Maize', 'Rice','Arabidopsis', 'Tomato']
df = pd.DataFrame(columns=['NHEJ', 'MMEJ', 'boot_n', 'snap', 'loop', 'SDMMEJ', 'unclassified_ins','Organism'])
for d,o in zip(datasets,org):
    df_pth = f"{pth}/{d}_EMmej_overall_estimations.tsv"
    tmp_df = pd.read_csv(df_pth, sep='\t')
    tmp_df['Organism'] = o
    df = pd.concat([df, tmp_df])

print(df)

mechanisms = ['NHEJ', 'MMEJ', 'snap', 'loop', 'SDMMEJ', 'unclassified_ins']
sum_df = pd.DataFrame(columns=['Organism', 'mechanism','mean', 
                        'ci_001', 'ci_005','ci_snap_sdmmej_001','ci_snap_sdmmej_005',
                        'ci_snap_sdmmej_loop_001', 'ci_snap_sdmmej_loop_005'])
for o in org:
    tmp_df = df.loc[df['Organism'] == o, :]
    for m in mechanisms:
        i = sum_df.shape[0]
        sum_df.loc[i,'Organism'] = o
        sum_df.loc[i,'mechanism'] = m
        sum_df.loc[i,'mean'] = tmp_df[m].mean()
        sum_df.loc[i,'ci_005'] = [(tmp_df[m].mean() - tmp_df[m].quantile(q=0.025)), (tmp_df[m].quantile(q=0.975) - tmp_df[m].mean())]
        sum_df.loc[i,'ci_001'] = [(tmp_df[m].mean() - tmp_df[m].quantile(q=0.005)), (tmp_df[m].quantile(q=0.995) - tmp_df[m].mean())]
        
        snap_sdmmej = tmp_df['snap'] + tmp_df['SDMMEJ']
        sum_df.loc[i,'ci_snap_sdmmej_005'] = [(snap_sdmmej.mean() - snap_sdmmej.quantile(q=0.025)), 
                                            (snap_sdmmej.quantile(q=0.975) - snap_sdmmej.mean())]
        sum_df.loc[i,'ci_snap_sdmmej_001'] = [(snap_sdmmej.mean() - snap_sdmmej.quantile(q=0.005)), 
                                            (snap_sdmmej.quantile(q=0.995) - snap_sdmmej.mean())]
        
        snap_sdmmej_loop = tmp_df['snap'] + tmp_df['SDMMEJ'] + tmp_df['loop']
        sum_df.loc[i,'ci_snap_sdmmej_loop_005'] = [(snap_sdmmej_loop.mean() - snap_sdmmej_loop.quantile(q=0.025)), 
                                            (snap_sdmmej_loop.quantile(q=0.975) - snap_sdmmej_loop.mean())]
        sum_df.loc[i,'ci_snap_sdmmej_loop_001'] = [(snap_sdmmej_loop.mean() - snap_sdmmej_loop.quantile(q=0.005)), 
                                            (snap_sdmmej_loop.quantile(q=0.995) - snap_sdmmej_loop.mean())]
        
sum_df.loc[sum_df['mechanism'] == 'unclassified_ins', 'mechanism'] = "Non_SDMMEJ"
print(sum_df)
sum_df.to_csv(f"{out_path}/sum_df.tsv", sep='\t', index=False)

barWidth = 0.3
linewidth=0.8
elinewidth=0.8
top = 1.05
bot = 0
al=0.75
error_kw=dict(lw=1, capsize=3, capthick=1)
fs=16

r1 = np.arange(len(mechanisms)-1) 
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]
r4 = [x + barWidth for x in r3]
r5 = [x + barWidth for x in r4]
r6 = [x + barWidth for x in r5]

# ax = figure(figsize=(16, 6)) # num=None,
ax = figure(num=None, figsize=(16, 6))
plt.subplot(1, 3, 1)
ci = sum_df.loc[sum_df['mechanism'] == 'MMEJ', 'ci_005'].to_list()
plt.bar(r1, sum_df.loc[(sum_df['mechanism'] == 'MMEJ'), 'mean']
            , width=barWidth, color='#ca0020', alpha=al, #d7191c
            yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,#hatch="x",
                error_kw=error_kw, label='MMEJ')
ci = sum_df.loc[sum_df['mechanism'] == 'NHEJ', 'ci_005'].to_list()
plt.bar(r2, sum_df.loc[(sum_df['mechanism'] == 'NHEJ'), 'mean']
            , width=barWidth, color='#f4a582',alpha=al, #  #d6adad
            yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,#hatch="x",
                error_kw=error_kw, label='cNHEJ')
plt.legend(bbox_to_anchor=(1.01, 1))
plt.ylim(bot, top)
plt.xticks([r + (barWidth*0.5) for r in range(len(org))], 
    org, rotation = 45)

plt.subplot(1, 3, 2)
ci = sum_df.loc[sum_df['mechanism'] == 'loop', 'ci_005'].to_list()
plt.bar(r1, sum_df.loc[(sum_df['mechanism'] == 'loop'), 'mean']
            , width=barWidth, color='#80cdc1', alpha=al, #009966 #ffffbf '#f7f7f7'
            yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label='SD-Loop-Out')

ci = sum_df.loc[sum_df['mechanism'] == 'snap', 'ci_005'].to_list()
plt.bar(r2, sum_df.loc[(sum_df['mechanism'] == 'snap'), 'mean']
            , width=barWidth, color='#92c5de', alpha=al, # #662222 #cc80ff
            yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label='SD-Snap-Back')


ci = sum_df.loc[sum_df['mechanism'] == 'Non_SDMMEJ', 'ci_005'].to_list()
plt.bar(r3, sum_df.loc[(sum_df['mechanism'] == 'Non_SDMMEJ'), 'mean']
            , width=barWidth, color='#0571b0',alpha=al,  #006699
            yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,  
                error_kw=error_kw, label='Non SD-MMEJ')
plt.hlines(0.05, -0.08,4.75, linestyles='dotted', colors='red', alpha=0.7)
plt.xticks([r + barWidth for r in range(len(org))], 
    org, rotation = 45)
plt.ylim(bot, top)
plt.legend(bbox_to_anchor=(1.01, 1))
fig = ax.get_figure()
fig.savefig(f"{out_path}/tst1.svg", bbox_inches='tight')



# plt_df = pd.melt(df, id_vars=['Organism'], value_vars=['MMEJ',
#             'NHEJ', 'snap', 'loop', 'SDMMEJ','unclassified_ins'])
# print(plt_df)
# ax = plt.figure(figsize = (7,5), dpi=80) # , layout='tight'
# ax = sns.boxplot(data=plt_df, # df.loc[plt_df['variable'].isin(['MMEJ'])]
#                  y='value', x='variable', 
#                  linewidth=1, hue='Organism')  # , order= mechanisms , palette=color_dict

# fig = ax.get_figure()
# fig.savefig(f"{out_path}/tst1.svg", bbox_inches='tight')



# Stacked barplot with 0.05 ci
ax = figure(figsize=(16, 6)) # num=None,
ci = sum_df.loc[sum_df['mechanism'] == 'MMEJ', 'ci_005'].to_list()
print(ci)
mmej = sum_df.loc[(sum_df['mechanism'] == 'MMEJ'), 'mean']
plt.bar(r1, mmej
            , width=barWidth, color='#08519c', alpha=al, #ca0020  #a50f15
            yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,#hatch="x",
                error_kw=error_kw, label='MMEJ')


nhej=sum_df.loc[(sum_df['mechanism'] == 'NHEJ'), 'mean']
plt.bar(r1,nhej , bottom=mmej
            , width=barWidth, color='#fed976',alpha=al, #  #d6adad #f4a582 #6baed6
                edgecolor='black',linewidth=linewidth ,#hatch="x",
                error_kw=error_kw, label='cNHEJ')

# Snap
ci = sum_df.loc[sum_df['mechanism'] == 'snap', 'ci_005'].to_list()
snap = sum_df.loc[(sum_df['mechanism'] == 'snap'), 'mean']
plt.bar(r2, snap
            , width=barWidth, color='#016c59', alpha=al, # #662222 #cc80ff #92c5de #de2d26  #3182bd
            yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label='SD-Snap-Back')
# SD-MMEJ

SDMMEJ=sum_df.loc[(sum_df['mechanism'] == 'SDMMEJ'), 'mean']
# ci = sum_df.loc[sum_df['mechanism'] == 'Non_SDMMEJ', 'ci'].to_list()
ci = sum_df.loc[sum_df['mechanism'] == 'Non_SDMMEJ', 'ci_snap_sdmmej_005'].to_list()
plt.bar(r2, SDMMEJ, bottom=snap
            , width=barWidth, color='#1c9099',alpha=al,  #006699 #018571 #fb6a4a #6baed6
            yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,  
                error_kw=error_kw, label='SD-MMEJ')

# Loop
# ci = sum_df.loc[sum_df['mechanism'] == 'loop', 'ci'].to_list()
ci = sum_df.loc[sum_df['mechanism'] == 'loop', 'ci_snap_sdmmej_loop_005'].to_list()

loop = sum_df.loc[(sum_df['mechanism'] == 'loop'), 'mean']
plt.bar(r2, loop,bottom=np.add(SDMMEJ.to_list(),snap.to_list()).tolist()
            , width=barWidth, color='#9ecae1', alpha=al, #009966 #ffffbf '#f7f7f7' #80cdc1 #fc9272
            yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label='SD-Loop-Out')


# Non SD-MMEJ
Non_SDMMEJ=sum_df.loc[(sum_df['mechanism'] == 'Non_SDMMEJ'), 'mean']
# ci = sum_df.loc[sum_df['mechanism'] == 'Non_SDMMEJ', 'ci'].to_list()
sdmmej_bars = np.add(loop.to_list(),snap.to_list()).tolist()
sdmmej_bars = np.add(sdmmej_bars,SDMMEJ.to_list()).tolist()

plt.bar(r2, Non_SDMMEJ, bottom=sdmmej_bars
            , width=barWidth, color='#ffffb2',alpha=al,  #006699 #0571b0 #9ecae1
                edgecolor='black',linewidth=linewidth ,  
                error_kw=error_kw, label='Non SD-MMEJ')

plt.legend(bbox_to_anchor=(1.01, 1))
plt.ylim(bot, top)
plt.xticks([r + (barWidth*0.5) for r in range(len(org))], 
    org, rotation = 45)

plt.hlines(0.05, -0.15,4.45, linestyles='dotted', colors='black', alpha=0.7)

fig = ax.get_figure()
fig.savefig(f"{out_path}/muli_species_non_repeated_regions_mechanisms_prop_ci_005_with_dashed.svg", bbox_inches='tight')
# fig.savefig(f"{out_path}/muli_species_non_repeated_regions_mechanisms_prop_ci_005.svg", bbox_inches='tight')




# Stacked barplot with 0.01 ci
ax = figure(figsize=(16, 6)) # num=None,
ci = sum_df.loc[sum_df['mechanism'] == 'MMEJ', 'ci_001'].to_list()
print(ci)
mmej = sum_df.loc[(sum_df['mechanism'] == 'MMEJ'), 'mean']
plt.bar(r1, mmej
            , width=barWidth, color='#08519c', alpha=al, #ca0020  #a50f15
            yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,#hatch="x",
                error_kw=error_kw, label='MMEJ')


nhej=sum_df.loc[(sum_df['mechanism'] == 'NHEJ'), 'mean']
plt.bar(r1,nhej , bottom=mmej
            , width=barWidth, color='#fed976',alpha=al, #  #d6adad #f4a582 #6baed6
            # yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,#hatch="x",
                error_kw=error_kw, label='cNHEJ')

# Snap
ci = sum_df.loc[sum_df['mechanism'] == 'snap', 'ci_001'].to_list()
snap = sum_df.loc[(sum_df['mechanism'] == 'snap'), 'mean']
plt.bar(r2, snap
            , width=barWidth, color='#016c59', alpha=al, # #662222 #cc80ff #92c5de #de2d26  #3182bd
            yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label='SD-Snap-Back')
# SD-MMEJ

SDMMEJ=sum_df.loc[(sum_df['mechanism'] == 'SDMMEJ'), 'mean']
# ci = sum_df.loc[sum_df['mechanism'] == 'Non_SDMMEJ', 'ci'].to_list()
ci = sum_df.loc[sum_df['mechanism'] == 'Non_SDMMEJ', 'ci_snap_sdmmej_001'].to_list()
plt.bar(r2, SDMMEJ, bottom=snap
            , width=barWidth, color='#1c9099',alpha=al,  #006699 #018571 #fb6a4a #6baed6
            yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,  
                error_kw=error_kw, label='SD-MMEJ')

# Loop
# ci = sum_df.loc[sum_df['mechanism'] == 'loop', 'ci'].to_list()
ci = sum_df.loc[sum_df['mechanism'] == 'loop', 'ci_snap_sdmmej_loop_001'].to_list()

loop = sum_df.loc[(sum_df['mechanism'] == 'loop'), 'mean']
plt.bar(r2, loop,bottom=np.add(SDMMEJ.to_list(),snap.to_list()).tolist()
            , width=barWidth, color='#9ecae1', alpha=al, #009966 #ffffbf '#f7f7f7' #80cdc1 #fc9272
            yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label='SD-Loop-Out')


# Non SD-MMEJ
Non_SDMMEJ=sum_df.loc[(sum_df['mechanism'] == 'Non_SDMMEJ'), 'mean']
# ci = sum_df.loc[sum_df['mechanism'] == 'Non_SDMMEJ', 'ci'].to_list()
sdmmej_bars = np.add(loop.to_list(),snap.to_list()).tolist()
sdmmej_bars = np.add(sdmmej_bars,SDMMEJ.to_list()).tolist()

plt.bar(r2, Non_SDMMEJ, bottom=sdmmej_bars
            , width=barWidth, color='#ffffb2',alpha=al,  #006699 #0571b0 #9ecae1
            # yerr = [[i[0] for i in ci], [i[1] for i in ci]], 
                edgecolor='black',linewidth=linewidth ,  
                error_kw=error_kw, label='Non SD-MMEJ')

plt.legend(bbox_to_anchor=(1.01, 1))
plt.ylim(bot, top)
plt.xticks([r + (barWidth*0.5) for r in range(len(org))], 
    org, rotation = 45)

plt.hlines(0.05, -0.15,4.45, linestyles='dotted', colors='black', alpha=0.7)

fig = ax.get_figure()
fig.savefig(f"{out_path}/muli_species_non_repeated_regions_mechanisms_prop_ci_001_with_dashed.svg", bbox_inches='tight')
# fig.savefig(f"{out_path}/muli_species_non_repeated_regions_mechanisms_prop_ci_001.svg", bbox_inches='tight')
