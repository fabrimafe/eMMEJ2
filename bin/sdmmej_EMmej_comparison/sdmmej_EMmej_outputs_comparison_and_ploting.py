import os

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import seaborn as sns
# sns.set_style("whitegrid")
# sns.set_style(style="ticks")

pd.options.display.max_colwidth = 3500
# pd.set_option('display.max_rows', 1000)
pd.set_option("display.max_columns", None)


def get_likelier_mechanism(df:pd.DataFrame):
    """
    This function generates a dataset of most
    likely repair scenario out of the two directions
    and the three repair pathways
    that has been looked for (0 and 1) by taking the
    max of the two
    """
    local_df = df.copy()
    local_df['likeliest_mechanism'] = np.nan
    likelier_idx_ls = []
    # df['tmp_variant_id'] = df['variant_id'].str.slice(start=0, stop=-2)

    for variant in df['sliced_variant_id'].unique():
        variant_df = df.loc[(df['sliced_variant_id'] == variant), :]
        """
        In cases were snap or loop pattern are found only on one direction
        but not the other, then the other direction max lk is aoutomaticly
        will be unclasified_ins. Since we dont know if its the ground truth,
        we give priority to the direction at which a pattern was found. 
        """
        snap, loop = False, False
        if variant_df.loc[:, 'SD_snap_back'].any():
            snap=True
            variant_df_snap = variant_df.loc[variant_df['SD_snap_back'], :]
        if variant_df.loc[:, 'SD_loop_out'].any():
            loop=True
            variant_df_loop = variant_df.loc[variant_df['SD_loop_out'], :]
        if snap == loop == True:
            c_df = pd.concat([variant_df_snap, variant_df_loop], axis=0)
            c_df.drop_duplicates(inplace=True)
        elif (snap == True) & (loop == False): c_df = variant_df.loc[(variant_df['SD_snap_back'] == True), :]
        elif (snap == False) & (loop == True): c_df = variant_df.loc[(variant_df['SD_loop_out'] == True), :]
        else: c_df = variant_df
        variant_df = c_df

        # finding the mechanism with the highest lk
        likelyhood_dict = {}
        likelier_snap = variant_df.loc[
            (variant_df['SD_snap_back_lk'] == variant_df['SD_snap_back_lk'].max()), :]
        likelier_loop = variant_df.loc[
            (variant_df['SD_loop_out_lk'] == variant_df['SD_loop_out_lk'].max()), :]
        likelier_unclassified_ins_lk = variant_df.loc[
            (variant_df['unclassified_ins_lk'] == variant_df['unclassified_ins_lk'].max()), :]

        likelyhood_dict['SD_snap_back_lk'] = likelier_snap['SD_snap_back_lk'].to_list()[0]
        likelyhood_dict['SD_loop_out_lk'] = likelier_loop['SD_loop_out_lk'].to_list()[0]
        likelyhood_dict['unclassified_ins_lk'] = likelier_unclassified_ins_lk['unclassified_ins_lk'].to_list()[0]

        likeliest_mechanism = max(likelyhood_dict, key=likelyhood_dict.get)
        likeliest_mechanism_lk = max(likelyhood_dict.values())
        likelier_direction_idx = min(variant_df.loc[(variant_df[likeliest_mechanism] == likeliest_mechanism_lk), 
                             ['SD_snap_back_lk', 'SD_loop_out_lk', 'unclassified_ins_lk']].index.to_list())
        likelier_idx_ls.append(likelier_direction_idx)
        local_df.loc[likelier_direction_idx,'likeliest_mechanism'] = likeliest_mechanism          
                
    return local_df.loc[likelier_idx_ls, :]

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3
def parallel(lst1, lst2):
    lst3 = [value for value in lst1 if value not in lst2]
    return lst3

def norm_approx_interval(p:float,n:int):
    """
    This function calculates the confidence 
    interval using normal approximation interval
    (Wlad interval) using the following formula:
    p^ = z*((p*(1-p)/n)**0.5)
    were z=1.96 for 0.95 confidence interval.
    """
    ci_095 = 1.96*((p*(1-p)/n)**0.5)
    # return (p-ci_095, p+ci_095)
    return np.array([ci_095, ci_095], dtype=object)


# date = '20221028'
# date = '20221030'
date = '20221117'
date = '20221129'
date = '20221130'
run_path = '/home/labs/alevy/guyta/guy_master_project/results/sdmmej_EMmej_comparison/'
EMmej_output_path = f'{run_path}/{date}'
# path_to_save_plots = f'{run_path}plots/{date}'
# EMmej_output_path = f'{run_path}/20221201'
path_to_save_plots = f'{run_path}plots/20221201'
path_to_save_plots = f'{run_path}plots/20221222'
path_to_save_plots = f'{run_path}plots/20221224'


dir = os.path.join(f"{path_to_save_plots}")
if not os.path.exists(dir):
    os.mkdir(dir)

mutants_res_df = pd.DataFrame(columns=['mechanism','indel_type', 'Algorithm', 'Proportions', 'ci','Mutant'])
mutants = ['R0_WT', 'R0_POLQ', 'R0_Lig4'] 
for mut in mutants: 
    print(f'-------------- {mut} ---------------')
    EMmej_input_tsv_path = f'/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/{mut}/20221030/20221030_{mut}_sdmmej_output_EMmej_input.tsv'
    file_prefix = '20221030'
    # file_prefix = date
    RMdetector_output = f'{EMmej_output_path}/{file_prefix}_{mut}_sdmmej_output_EMmej_input_EMmej_output/{file_prefix}_{mut}_sdmmej_output_EMmej_input_RMdetector_output.tsv' 
    del_repair_proportions = f'{EMmej_output_path}/{file_prefix}_{mut}_sdmmej_output_EMmej_input_EMmej_output/{file_prefix}_{mut}_sdmmej_output_EMmej_input_markov.tsv'
    path_to_MM_output = f'{EMmej_output_path}/{file_prefix}_{mut}_sdmmej_output_EMmej_input_EMmej_output/{file_prefix}_{mut}_sdmmej_output_EMmej_input_markov_output.tsv' 

    EMmej_input_tsv = pd.read_csv(EMmej_input_tsv_path, sep='\t')
    EMmej_input_tsv['sliced_variant_id'] = EMmej_input_tsv.apply(lambda row: f"{row['#CHR']}_{row['POS']}", axis=1)

    emmej_res = pd.read_csv(RMdetector_output, sep='\t')
    MM_output = pd.read_csv(path_to_MM_output, sep='\t')
    EMmej_RM_MM = emmej_res.join(MM_output, rsuffix='_MM')
    EMmej_RM_MM.drop(columns = ['CHR_MM', 'POS_MM',
           'original_pos_MM', 'variant_id_MM', 'direction_MM',
           'ancestral_indel_MM', 'derived_indel_MM'], inplace=True)
    
    EMmej_RM_MM['sliced_variant_id'] = EMmej_RM_MM.apply(lambda row: f"{row['CHR']}_{row['POS']}", axis=1)
    
    # Merging the EMmej and sdmmej results
    # emmej_sdmmej_merged = EMmej_RM_MM.merge(EMmej_input_tsv, left_on='CHR', right_on='#CHR')
    emmej_sdmmej_merged = EMmej_RM_MM.merge(EMmej_input_tsv, left_on='sliced_variant_id', right_on='sliced_variant_id')
    emmej_sdmmej_merged.drop(columns=['REF', 'ALT', 'ID'], inplace=True)
    
    emmej_sdmmej_merged.to_csv(
            f"{path_to_save_plots}/{mut}_emmej_sdmmej_merged.tsv", index=False, sep='\t')
    emmej_sdmmej_merged.loc[emmej_sdmmej_merged['indel_type'] == 'DEL',
        ['del_mmej_cand', 'MICROHOMOLOGY']].to_csv(
            f"{path_to_save_plots}/{mut}_MMEJ_patterns.tsv", index=False, sep='\t')
    
    emmej_sdmmej_merged['sliced_variant_id'] = emmej_sdmmej_merged['variant_id'].str.slice(start=0, stop=-2)
    print(emmej_sdmmej_merged['sliced_variant_id'].unique().shape)
    # emmej_sdmmej_merged['loop_dist_between_reps'].fillna(value=1, inplace=True)
    # emmej_sdmmej_merged['snap_dist_between_reps'].fillna(value=1, inplace=True)
    # emmej_sdmmej_merged = emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['snap_dist_between_reps']<65) &
    #                      (emmej_sdmmej_merged['loop_dist_between_reps']<65)), :]
    # print(emmej_sdmmej_merged['sliced_variant_id'].unique().shape)


    # if mut == 'R0_WT':
    #     # emmej_sdmmej_merged.to_csv(f"{path_to_save_plots}/emmej_sdmmej_merged.tsv", sep='\t', index=False)
    #     # print(emmej_sdmmej_merged.loc[emmej_sdmmej_merged['indel_type'] == 'DEL', 'del_mmej_lk'].describe())
    #     a = figure(num=None, figsize=(16, 6))
    #     sns.histplot(
    #         emmej_sdmmej_merged.loc[emmej_sdmmej_merged['indel_type'] == 'DEL', 'del_mmej_lk']
    #         ,bins=100, stat='proportion')
    #     plt.title(f'Deletion MMEJ Probability distribution for R0-WT')
    #     fig = a.get_figure()
    #     fig.savefig(f"{path_to_save_plots}/R0_WT_del_mmej_lk_PDF.svg", bbox_inches='tight')
    #     c = ['CHR', 'POS_x', 'original_pos', 'variant_id', 'direction',
    #    'ancestral_indel', 'derived_indel', 'indel_type', 'indel_len',
    #    'del_mmej', 'del_mmej_cand', 'del_mmej_marked',
    #    'del_mmej_marked_on_ref','del_mmej_lk','NHEJ_lk','REPAIR_TYPE', 'consistency', 'CLASS',
    #    'CLASS_final', 'READS', 'MICROHOMOLOGY', 'MH_Length']
        
        # emmej_sdmmej_merged.loc[
        #     emmej_sdmmej_merged['indel_type'] == 'DEL', c].to_csv(f"{path_to_save_plots}/R0_WT_deletions.tsv",
        #         sep='\t', index= False)
    del_n = emmej_sdmmej_merged.loc[emmej_sdmmej_merged['indel_type'] == 'DEL','sliced_variant_id'].unique().shape[0]
    ins_n = emmej_sdmmej_merged.loc[(emmej_sdmmej_merged['indel_type'] == 'INS'), 'sliced_variant_id'].unique().shape[0]
    
    if mut == 'R0_WT': WT_del_n, WT_ins_n = del_n, ins_n
    if mut == 'R0_POLQ': POLQ_del_n, POLQ_ins_n = del_n, ins_n
    if mut == 'R0_Lig4': Lig_4_del_n, Lig_4_ins_n = del_n, ins_n
    
    """
    EMmej results: discrete meassuring (counting mechanisms using the boolian indications)
    """
    # Deletions:
    # making a list of cases of which del MMEJ was detected on at list one direction
    mmej = emmej_sdmmej_merged.loc[
        ((emmej_sdmmej_merged['indel_type'] == 'DEL') &
         (emmej_sdmmej_merged['del_mmej'] == True)), 'CHR'].unique().tolist()
    EMmej_mmej = emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['indel_type'] == 'DEL') & 
                             (emmej_sdmmej_merged['del_mmej'] == True) & 
                            (emmej_sdmmej_merged['CHR'].isin(mmej))), :]
    EMmej_mmej_count = EMmej_mmej['CHR'].unique().shape[0]
    EMmej_mmej_count_prop = EMmej_mmej_count/del_n

    EMmej_nhej = emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['indel_type'] == 'DEL') & 
                             (emmej_sdmmej_merged['del_mmej'] == False) & 
                            (~emmej_sdmmej_merged['CHR'].isin(mmej))), :]
    EMmej_nhej_count = EMmej_nhej['CHR'].unique().shape[0]
    EMmej_nhej_count_prop = EMmej_nhej_count/del_n
    
    # Insertions
    likeliest_ins = get_likelier_mechanism(df=emmej_sdmmej_merged.loc[(emmej_sdmmej_merged['indel_type'] == 'INS'),:])
    
    likeliest_ins.to_csv(f"{path_to_save_plots}/{mut}_likeliest_ins.tsv", index=False, sep='\t')
    likeliest_ins.loc[likeliest_ins['indel_type'] == 'INS',
        ['direction','snap_mmej_marked', 'snap_repeat_pat', 'snap_inv_comp_repeat_pat','Snap-back']].to_csv(
            f"{path_to_save_plots}/{mut}_snap_patterns.tsv", index=False, sep='\t')

    likeliest_ins.loc[likeliest_ins['indel_type'] == 'INS',
        ['direction','loop_mmej_marked', 'loop_repeat_pat', 'Loop-out']].to_csv(
            f"{path_to_save_plots}/{mut}_loop_patterns.tsv", index=False, sep='\t')



    # Snap back
    snap = likeliest_ins.loc[
        ((likeliest_ins['indel_type'] == 'INS') &
            (likeliest_ins['SD_snap_back'] == True)),
            'variant_id'].unique().tolist()

    # loop out
    loop = likeliest_ins.loc[
        ((likeliest_ins['indel_type'] == 'INS') &
            (likeliest_ins['SD_loop_out'] == True)),
            'variant_id'].unique().tolist()

    snap_only = parallel(snap, loop)
    loop_only = parallel(loop, snap)
    snap_loop = intersection(snap, loop)
    
    EMmej_snap = likeliest_ins.loc[likeliest_ins['variant_id'].isin(snap_only), :]
    EMmej_snap_est = (EMmej_snap['SD_snap_back_lk'].sum()/ins_n)
    EMmej_snap_count = len(parallel(snap,loop))
    EMmej_snap_count_prop = EMmej_snap_count/ins_n

    # Loop out
    EMmej_loop = likeliest_ins.loc[likeliest_ins['variant_id'].isin(loop_only), :]
    EMmej_loop_est = (EMmej_loop['SD_loop_out_lk'].sum()/ins_n)
    EMmej_loop_count = len(parallel(loop, snap))
    EMmej_loop_count_prop = EMmej_loop_count/ins_n

    # SD-MMEJ (Snap + Loop)
    # EMmej_snap_loop = pd.concat([EMmej_snap, EMmej_loop], axis=0)
    EMmej_snap_loop = likeliest_ins.loc[likeliest_ins['variant_id'].isin(snap_loop), :]
    EMmej_snap_loop_est = ((EMmej_snap_loop['SD_snap_back_lk'].sum()/ins_n) + (EMmej_snap_loop['SD_loop_out_lk'].sum()/ins_n))
    EMmej_snap_loop_count = len(intersection(snap, loop))
    EMmej_snap_loop_count_prop = EMmej_snap_loop_count/ins_n

    EMmej_SDMMEJ_prop = EMmej_snap_count_prop + EMmej_loop_count_prop + EMmej_snap_loop_count_prop
    EMmej_SDMMEJ_est = EMmej_snap_est+EMmej_loop_est+EMmej_snap_loop_est

    # Unclassified ins
    # Unclassifed insertions
    EMmej_unclassified_ins = likeliest_ins.loc[~likeliest_ins['variant_id'].isin(snap_only+loop_only+snap_loop),:]
    EMmej_unclassified_ins_est = (EMmej_unclassified_ins['unclassified_ins_lk'].sum()/ins_n)
    EMmej_unclassified_ins_count = EMmej_unclassified_ins['variant_id'].unique().shape[0]
    EMmej_unclassified_ins_count_prop = EMmej_unclassified_ins_count/ins_n

    """
    EMmej results: Estimating quantities of each mechanism.
        For deletions: quantities are given by EM estimations * number of deletions
        For insertions: quantities are given by summing the whole lk columns, normalizing them
            to get proportions and then * number of insertions
    """
    # Deletions using EM
    if mut != 'R0_POLQ':
        del_proportions = pd.read_csv(del_repair_proportions, sep='\t')
        EMmej_mmej_proportion = del_proportions.loc[del_proportions.index.max(), 'MMEJ'] 
        EMmej_nhej_proportion = del_proportions.loc[del_proportions.index.max(), 'NHEJ'] 
    else: EMmej_mmej_proportion, EMmej_nhej_proportion = np.nan, np.nan
    # Deletions using sum of lk
    mmej_lk_sum = emmej_sdmmej_merged.loc[emmej_sdmmej_merged['indel_type'] == 'DEL', 'del_mmej_lk'].sum()
    nhej_lk_sum = emmej_sdmmej_merged.loc[emmej_sdmmej_merged['indel_type'] == 'DEL', 'NHEJ_lk'].sum()
    norm_mmej_lk_sum = (mmej_lk_sum/(mmej_lk_sum+nhej_lk_sum))*del_n
    norm_nhej_lk_sum = (nhej_lk_sum/(mmej_lk_sum+nhej_lk_sum))*del_n

    """
    sdmmej results:
    """
    # Deletions
    sdmmej_Consistent_MHJ_reads = emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['indel_type'] == 'DEL') & 
                                                 (emmej_sdmmej_merged['REPAIR_TYPE'] == 'MHJ') &
                                                 (emmej_sdmmej_merged['consistency'] == True)), 'sliced_variant_id'].unique().shape[0]
    sdmmej_MMEJ_reads = emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['indel_type'] == 'DEL') & 
                                                 (emmej_sdmmej_merged['REPAIR_TYPE'] == 'MHJ') &
                                                 (emmej_sdmmej_merged['consistency'] == False)), 'sliced_variant_id'].unique().shape[0]

    sdmmej_mmej_reads = sdmmej_Consistent_MHJ_reads + sdmmej_MMEJ_reads
    sdmmej_mmej_reads_prop = sdmmej_mmej_reads/del_n
    sdmmej_Consistent_ABJ_reads = emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['indel_type'] == 'DEL') & 
                                                 (emmej_sdmmej_merged['REPAIR_TYPE'] == 'ABJ') &
                                                 (emmej_sdmmej_merged['consistency'] == True)), 'sliced_variant_id'].unique().shape[0]
    sdmmej_NHEJ_reads = emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['indel_type'] == 'DEL') & 
                                                 (emmej_sdmmej_merged['REPAIR_TYPE'] == 'ABJ') &
                                                 (emmej_sdmmej_merged['consistency'] == False)), 'sliced_variant_id'].unique().shape[0]

    sdmmej_nhej_abj_reads = sdmmej_Consistent_ABJ_reads + sdmmej_NHEJ_reads
    sdmmej_nhej_abj_reads_prop = sdmmej_nhej_abj_reads/del_n

    # Insertions
    emmej_sdmmej_merged['SDMMEJ_snap'] = False
    emmej_sdmmej_merged['SDMMEJ_loop'] = False
    emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['Snap-back'] !='0') & (emmej_sdmmej_merged['Snap-back'].isna() == False)), 'SDMMEJ_snap'] = True
    emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['Loop-out'] !='0') & (emmej_sdmmej_merged['Loop-out'].isna() == False)), 'SDMMEJ_loop'] = True

    SDMMEJ_snap_reads = emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['SDMMEJ_snap'] == True) & 
                                                 (emmej_sdmmej_merged['indel_type'] == 'INS')), 'sliced_variant_id'].unique().shape[0]
    SDMMEJ_loop_reads = emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['SDMMEJ_loop'] == True) & 
                                                 (emmej_sdmmej_merged['indel_type'] == 'INS')), 'sliced_variant_id'].unique().shape[0]


    SDMMEJ_snap = emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['indel_type'] == 'INS') & (emmej_sdmmej_merged['SDMMEJ_snap'] == True)), 'sliced_variant_id'].unique().tolist()
    SDMMEJ_loop = emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['indel_type'] == 'INS') & (emmej_sdmmej_merged['SDMMEJ_loop'] == True)), 'sliced_variant_id'].unique().tolist()
    SDMMEJ_snap_count_prop = len(parallel(SDMMEJ_snap, SDMMEJ_loop))/ins_n
    SDMMEJ_loop_count_prop = len(parallel(SDMMEJ_loop, SDMMEJ_snap))/ins_n
    SDMMEJ_snap_and_loop_count_prop = len(intersection(SDMMEJ_loop, SDMMEJ_snap))/ins_n
    SDMMEJ_SDMMEJ_prop = SDMMEJ_snap_and_loop_count_prop+SDMMEJ_snap_count_prop+SDMMEJ_loop_count_prop
    Unclassified_ins_Insertion_sdmmej = (emmej_sdmmej_merged.loc[((emmej_sdmmej_merged['indel_type'] == 'INS') &
                (~emmej_sdmmej_merged['sliced_variant_id'].isin(
                    intersection(SDMMEJ_loop, SDMMEJ_snap)+SDMMEJ_snap+SDMMEJ_loop
                    ))), 'sliced_variant_id'].unique().shape[0]/ins_n)

    # ---------------------------------------------------------------------------------------------
    # Insertions using Bonferroni
    emmej_sdmmej_merged['variant_id'] = emmej_sdmmej_merged['variant_id'].str.slice(0,-2)
    emmej_sdmmej_merged['variant_id_N'] = 2

    def parallel(lst1, lst2):
        lst3 = [value for value in lst1 if value not in lst2]
        return lst3

    snap_ins_cutoff = 0.005
    loop_ins_cutoff = 0.003
    # Isolate insertons:
    ins_df = emmej_sdmmej_merged.loc[
        (emmej_sdmmej_merged['ancestral_indel'].str.len() < emmej_sdmmej_merged['derived_indel'].str.len()),
         :].copy()
    insertions_proportions = pd.DataFrame(columns=['snap', 'loop', 'SD-MMEJ','unclassified_ins'])
    if ins_df.shape[0] == 0:
        snap_markov_proportion = np.nan
        loop_markov_proportion = np.nan
        SDMMEJ_markov_proportion = np.nan
    else:
        # Snap back
        snap_df = ins_df.loc[(~ins_df['snap_repeat_pat'].isna()), :].copy()
        snap_df['bonferroni_cutoff'] = snap_ins_cutoff / snap_df['variant_id_N']
        snap_variant_ids = snap_df.loc[(snap_df['SD_snap_back_p_val'] <= snap_df['bonferroni_cutoff']), 
                            'variant_id'].unique()
        
        # Loop out
        loop_df = ins_df.loc[(~ins_df['loop_repeat_pat'].isna()), :].copy()
        loop_df['bonferroni_cutoff'] = loop_ins_cutoff / loop_df['variant_id_N']
        loop_variant_ids = loop_df.loc[(loop_df['SD_loop_out_p_val'] <= loop_df['bonferroni_cutoff']), 
                            'variant_id'].unique()
        
        # Crearing mutually exclusive variant_id lists
        only_snap = parallel(lst1=snap_variant_ids, lst2=loop_variant_ids)
        only_loop = parallel(lst1=loop_variant_ids, lst2=snap_variant_ids)
        SDMMEJ = (set(pd.Series(snap_variant_ids)) & set(pd.Series(loop_variant_ids)))
                
        
        snap_markov_proportion = len(only_snap)/ins_df['variant_id'].unique().shape[0]
        loop_markov_proportion = len(only_loop)/ins_df['variant_id'].unique().shape[0]
        SDMMEJ_markov_proportion = len(SDMMEJ)/ins_df['variant_id'].unique().shape[0]
        

    insertions_proportions.loc[0,'snap'] = snap_markov_proportion
    insertions_proportions.loc[0,'loop'] = loop_markov_proportion
    insertions_proportions.loc[0,'SD-MMEJ'] = SDMMEJ_markov_proportion
    insertions_proportions.loc[0,'unclassified_ins'] = (1 - snap_markov_proportion - loop_markov_proportion - SDMMEJ_markov_proportion)

    # ---------------------------------------------------------------------------------------------

    mechanism_proportion = pd.DataFrame(columns=['mechanism','indel_type', 'Algorithm', 'Proportions', 'ci','Mutant'])

    mechanism_proportion.loc[0,:] = 'NHEJ', 'Deletion', 'sdmmej' , sdmmej_nhej_abj_reads_prop, np.nan,mut
    mechanism_proportion.loc[0,'ci'] = norm_approx_interval(p=sdmmej_nhej_abj_reads_prop, n=del_n)
    mechanism_proportion.loc[1,:] = 'NHEJ', 'Deletion', 'EMmej' , EMmej_nhej_count_prop, np.nan, mut
    mechanism_proportion.loc[1,'ci'] =norm_approx_interval(p=EMmej_nhej_count_prop, n=del_n)

    mechanism_proportion.loc[2,:] = 'NHEJ', 'Deletion', 'EMmej proportions estimation' , EMmej_nhej_proportion, np.nan, mut
    
    mechanism_proportion.loc[3,:] = 'MMEJ', 'Deletion', 'sdmmej' , sdmmej_mmej_reads_prop, np.nan,mut
    mechanism_proportion.loc[3,'ci'] = norm_approx_interval( p=sdmmej_mmej_reads_prop, n=del_n)
    mechanism_proportion.loc[4,:] = 'MMEJ', 'Deletion', 'EMmej' , EMmej_mmej_count_prop, np.nan,mut
    mechanism_proportion.loc[4,'ci'] =norm_approx_interval(p=EMmej_mmej_count_prop, n=del_n)
    mechanism_proportion.loc[5,:] = 'MMEJ', 'Deletion', 'EMmej proportions estimation' , EMmej_mmej_proportion, np.nan, mut
    
    mechanism_proportion.loc[6,:] = 'Snap back', 'Insertion', 'sdmmej' , SDMMEJ_snap_count_prop/SDMMEJ_SDMMEJ_prop, np.nan,mut
    mechanism_proportion.loc[6,'ci'] =norm_approx_interval(p=SDMMEJ_snap_count_prop/SDMMEJ_SDMMEJ_prop, n=ins_n)

    mechanism_proportion.loc[7,:] = 'Snap back', 'Insertion', 'EMmej', EMmej_snap_count_prop/EMmej_SDMMEJ_prop, np.nan, mut
    mechanism_proportion.loc[7,'ci'] =norm_approx_interval(p=EMmej_snap_count_prop/EMmej_SDMMEJ_prop, n=ins_n)
    # mechanism_proportion.loc[8,:] = 'Snap back', 'Insertion', 'EMmej proportions estimation' , EMmej_snap_est/EMmej_SDMMEJ_est, np.nan ,mut
    mechanism_proportion.loc[8,:] = 'Snap back', 'Insertion', 'EMmej proportions estimation' , snap_markov_proportion, np.nan ,mut
    mechanism_proportion.loc[8,'ci'] =norm_approx_interval(p=snap_markov_proportion/EMmej_SDMMEJ_prop, n=ins_n)

    mechanism_proportion.loc[9,:] = 'Loop out', 'Insertion', 'sdmmej' , SDMMEJ_loop_count_prop/SDMMEJ_SDMMEJ_prop, np.nan, mut
    mechanism_proportion.loc[9,'ci'] =norm_approx_interval(p=SDMMEJ_loop_count_prop/SDMMEJ_SDMMEJ_prop, n=ins_n)
    mechanism_proportion.loc[10,:] = 'Loop out', 'Insertion', 'EMmej', EMmej_loop_count_prop/EMmej_SDMMEJ_prop, np.nan, mut
    mechanism_proportion.loc[10,'ci'] =norm_approx_interval(p=EMmej_loop_count_prop/EMmej_SDMMEJ_prop, n=ins_n)
    # mechanism_proportion.loc[11,:] = 'Loop out', 'Insertion', 'EMmej proportions estimation' , EMmej_loop_est/EMmej_SDMMEJ_est, np.nan, mut
    mechanism_proportion.loc[11,:] = 'Loop out', 'Insertion', 'EMmej proportions estimation' , loop_markov_proportion, np.nan, mut
    mechanism_proportion.loc[11,'ci'] =norm_approx_interval(p=loop_markov_proportion/EMmej_SDMMEJ_prop, n=ins_n)
    
    mechanism_proportion.loc[12,:] = 'Unclassified ins', 'Insertion', 'sdmmej' ,Unclassified_ins_Insertion_sdmmej , np.nan, mut
    mechanism_proportion.loc[12,'ci'] = norm_approx_interval(p=Unclassified_ins_Insertion_sdmmej, n=ins_n)

    
    mechanism_proportion.loc[13,:] = 'Unclassified ins', 'Insertion', 'EMmej' , EMmej_unclassified_ins_count_prop, np.nan, mut
    mechanism_proportion.loc[13,'ci'] =norm_approx_interval(p=EMmej_unclassified_ins_count_prop, n=ins_n)

    # mechanism_proportion.loc[14,:] = 'Unclassified ins', 'Insertion', 'EMmej proportions estimation',  EMmej_unclassified_ins_est, np.nan,mut
    mechanism_proportion.loc[14,:] = 'Unclassified ins', 'Insertion', 'EMmej proportions estimation',  (1 - snap_markov_proportion - loop_markov_proportion - SDMMEJ_markov_proportion), np.nan,mut

    mechanism_proportion.loc[15,:] = 'SD-MMEJ', 'Insertion', 'sdmmej' , SDMMEJ_SDMMEJ_prop, np.nan, mut
    mechanism_proportion.loc[15,'ci'] =norm_approx_interval(p=SDMMEJ_SDMMEJ_prop, n=ins_n)

    mechanism_proportion.loc[16,:] = 'SD-MMEJ', 'Insertion', 'EMmej' , EMmej_SDMMEJ_prop, np.nan, mut
    mechanism_proportion.loc[16,'ci'] =norm_approx_interval(p=EMmej_SDMMEJ_prop, n=ins_n)
    # mechanism_proportion.loc[17,:] = 'SD-MMEJ', 'Insertion', 'EMmej proportions estimation',  EMmej_SDMMEJ_est, np.nan,mut
    mechanism_proportion.loc[17,:] = 'SD-MMEJ', 'Insertion', 'EMmej proportions estimation',  (SDMMEJ_markov_proportion+loop_markov_proportion+snap_markov_proportion), np.nan,mut
    mechanism_proportion.loc[17,'ci'] =norm_approx_interval(p=(SDMMEJ_markov_proportion+loop_markov_proportion+snap_markov_proportion), n=ins_n) ######

    mechanism_proportion.loc[18,:] = 'Snap back and Loop out', 'Insertion', 'sdmmej' , SDMMEJ_snap_and_loop_count_prop/SDMMEJ_SDMMEJ_prop, np.nan, mut
    mechanism_proportion.loc[18,'ci'] = norm_approx_interval(p=SDMMEJ_snap_and_loop_count_prop/SDMMEJ_SDMMEJ_prop, n=ins_n)
    mechanism_proportion.loc[19,:] = 'Snap back and Loop out', 'Insertion', 'EMmej' , EMmej_snap_loop_count_prop/EMmej_SDMMEJ_prop, np.nan, mut
    mechanism_proportion.loc[19,'ci'] = norm_approx_interval(p=EMmej_snap_loop_count_prop/EMmej_SDMMEJ_prop, n=ins_n)
    # mechanism_proportion.loc[20,:] = 'Snap back and Loop out', 'Insertion', 'EMmej proportions estimation',  EMmej_snap_loop_est/EMmej_SDMMEJ_est, np.nan,mut
    mechanism_proportion.loc[20,:] = 'Snap back and Loop out', 'Insertion', 'EMmej proportions estimation',  (SDMMEJ_markov_proportion), np.nan,mut
    mechanism_proportion.loc[20,'ci'] = norm_approx_interval(p=SDMMEJ_markov_proportion, n=ins_n)
    
    print(mechanism_proportion)
    mutants_res_df = pd.concat([mutants_res_df, mechanism_proportion])
    
    if mut=='R0_WT':
        """
        First plot: 
        """
        barWidth = 0.3
        linewidth=0.8
        elinewidth=0.8
        top = 1.05
        bot = 0
        error_kw=dict(lw=1, capsize=3, capthick=1)
        fs=16
        # setting the x axis positions
        plt_df = mechanism_proportion.loc[mechanism_proportion['mechanism'].isin([ 'MMEJ', 'SD-MMEJ']), :] # 'NHEJ', 'Unclassified ins',
        first_plt_mech = ['MMEJ'] # 'NHEJ',
        r1 = np.arange(len(first_plt_mech))
        r2 = [x + barWidth for x in r1]
        r3 = [x + barWidth for x in r2]

        # Creating the bars for sdmmej
        # plt.style.use('seaborn-whitegrid')
        plt.style.use('seaborn-ticks')
        
        ax = figure(num=None, figsize=(16, 6))
        ax.suptitle(f"EMmej and sdmmej comparison, {mut}", fontsize=16)
        plt.subplot(1, 3, 1)
        sdmmej_ci = plt_df.loc[((plt_df['Algorithm'] == 'sdmmej') & (plt_df['indel_type']=='Deletion')), 'ci'].to_list()
        plt.bar(r1, plt_df.loc[((plt_df['Algorithm'] == 'sdmmej') & 
            (plt_df['indel_type']=='Deletion')), 'Proportions'], width=barWidth, color='#ef8a62', 
            yerr = [[i[0] for i in sdmmej_ci], [i[1] for i in sdmmej_ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label='sdmmej')
        
        # Creating the bars for non EMmej
        EMmej_ci = plt_df.loc[((plt_df['Algorithm'] == 'EMmej') & (plt_df['indel_type']=='Deletion')), 'ci'].to_list()
        plt.bar(r2, plt_df.loc[((plt_df['Algorithm'] == 'EMmej') & 
            (plt_df['indel_type']=='Deletion')), 'Proportions'], width=barWidth, color='#67a9cf', 
                yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label='Patterns') # label='EMmej'
        
        # Creating the bars for non EMmej proportions estimation
        plt.bar(r3, plt_df.loc[((plt_df['Algorithm'] == 'EMmej proportions estimation') & 
            (plt_df['indel_type']=='Deletion')), 'Proportions'], width=barWidth, color='#abdda4', 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label=f"EMmej") # label=f"EMmej $\hat{'P'}$"
        # Set plot parameters
        plt.xticks([r + barWidth for r in range(len(first_plt_mech))], first_plt_mech, rotation = 45)
        plt.ylabel(f'Proportion of MMEJ', fontsize=fs)
        plt.ylim(bot, top)
        plt.xlabel('')
        plt.title(f'Deletions')

        plt.subplot(1, 3, 2)
        first_plt_mech = ['SD-MMEJ'] # 'Unclassified ins',
        r1 = np.arange(len(first_plt_mech))
        r2 = [x + barWidth for x in r1]
        r3 = [x + barWidth for x in r2]
        sdmmej_ci = plt_df.loc[((plt_df['Algorithm'] == 'sdmmej') & (plt_df['indel_type']=='Insertion')), 'ci'].to_list()
        plt.bar(r1, plt_df.loc[((plt_df['Algorithm'] == 'sdmmej') & 
            (plt_df['indel_type']=='Insertion')), 'Proportions'], width=barWidth, color='#ef8a62', 
            yerr = [[i[0] for i in sdmmej_ci], [i[1] for i in sdmmej_ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label='sdmmej')
        
        # Creating the bars for non EMmej
        EMmej_ci = plt_df.loc[((plt_df['Algorithm'] == 'EMmej') & (plt_df['indel_type']=='Insertion')), 'ci'].to_list()
        plt.bar(r2, plt_df.loc[((plt_df['Algorithm'] == 'EMmej') & 
            (plt_df['indel_type']=='Insertion')), 'Proportions'], width=barWidth, color='#67a9cf', 
                yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label='Patterns') # label='EMmej'
        
        # Creating the bars for non EMmej proportions estimation
        EMmej_ci = plt_df.loc[((plt_df['Algorithm'] == 'EMmej proportions estimation') & (plt_df['indel_type']=='Insertion')), 'ci'].to_list()
        plt.bar(r3, plt_df.loc[((plt_df['Algorithm'] == 'EMmej proportions estimation') & 
            (plt_df['indel_type']=='Insertion')), 'Proportions'], width=barWidth, color='#abdda4',
            yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label=f"EMmej") # label=f"EMmej $\hat{'P'}$"

        
        # Set plot parameters
        plt.xticks([r + barWidth for r in range(len(first_plt_mech))], first_plt_mech, rotation = 45)
        plt.ylabel(f'Proportion of SD-MMEJ', fontsize=fs)
        plt.ylim(bot, top)
        plt.xlabel('')
        plt.title(f'Insertions')

        plt.subplot(1, 3, 3)
        SDMMEJ_zoomin_mech = ['Snap back','Loop out', 'Snap back and Loop out'] # 'SD-MMEJ',
        plt_df = mechanism_proportion.loc[mechanism_proportion['mechanism'].isin(SDMMEJ_zoomin_mech), :]
        # set x axis positions
        r1 = np.arange(len(SDMMEJ_zoomin_mech))
        r2 = [x + barWidth for x in r1]
        r3 = [x + barWidth for x in r2]
        r4 = [x + barWidth for x in r3]

        # Creating the bars for sdmmej
        plt.style.use('seaborn-ticks')
        
        sdmmej_ci = plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'ci'].to_list()
        plt.bar(r1, plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'Proportions'], 
            width=barWidth, color='#ef8a62', 
            yerr = [[i[0] for i in sdmmej_ci], [i[1] for i in sdmmej_ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label='sdmmej')
        
        # Creating the bars for non EMmej
        EMmej_ci = plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'ci'].to_list()
        plt.bar(r2, plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'Proportions'], 
                width=barWidth, color='#67a9cf', 
                yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label='Patterns') # label='EMmej'
        
        # Creating the bars for non EMmej proportions estimation
        EMmej_ci = plt_df.loc[plt_df['Algorithm'] == 'EMmej proportions estimation', 'ci'].to_list()
        print(EMmej_ci)
        plt.bar(r3, plt_df.loc[plt_df['Algorithm'] == 'EMmej proportions estimation', 'Proportions'], 
                width=barWidth, color='#abdda4', 
                yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
                edgecolor='black',linewidth=linewidth ,
                error_kw=error_kw, label=f"EMmej") # label=f"EMmej $\hat{'P'}$"

        # Set plot parameters
        plt.xticks([r + barWidth for r in range(len(SDMMEJ_zoomin_mech))], SDMMEJ_zoomin_mech, rotation = 45)
        plt.legend(bbox_to_anchor=(1.35, 1))
        plt.ylabel(f'Proportions (normalized to SD-MMEJ)', fontsize=fs)
        plt.ylim(bot, top)
        plt.xlabel('')
        plt.title(f'SD-MMEJ stratified to sub-mechanisms')
        fig = ax.get_figure()
        fig.savefig(f"{path_to_save_plots}/R0_WT_comparison_subplots.svg", bbox_inches='tight')

mutants_res_df.reset_index(drop=True, inplace=True)
mutants_res_df.to_csv(f"{path_to_save_plots}/mutants_plt_df.tsv", index=False, sep='\t')


"""
Mutant comparison:
"""
fs=16
third_plot_mech = [ 'SD-MMEJ'] # 'Unclassified ins',
plt_df = mutants_res_df.loc[((mutants_res_df['mechanism'].isin(third_plot_mech)) & 
            (mutants_res_df['Mutant'] == 'R0_WT')), :]
# set x axis positions
r1 = np.arange(len(third_plot_mech))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]

# Creating the bars for sdmmej
# plt.style.use('seaborn-darkgrid')
plt.style.use('seaborn-ticks')
ax = figure(num=None, figsize=(9, 6))
ax.suptitle(f"EMmej and sdmmej mutant comparison", fontsize=16)
plt.subplot(1, 3, 1)
sdmmej_ci = plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'ci'].to_list()
plt.bar(r1, plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'Proportions'], width=barWidth, color='#ef8a62', 
    yerr = [[i[0] for i in sdmmej_ci], [i[1] for i in sdmmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label='sdmmej')

# Creating the bars for non EMmej
EMmej_ci = plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'ci'].to_list()
plt.bar(r2, plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'Proportions'], width=barWidth, color='#67a9cf', 
        yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label='Patterns') # label='EMmej'

# Creating the bars for non EMmej proportions estimation
EMmej_ci = plt_df.loc[((plt_df['Algorithm'] == 'EMmej proportions estimation') & (plt_df['indel_type']=='Insertion')), 'ci'].to_list()
plt.bar(r3, plt_df.loc[plt_df['Algorithm'] == 'EMmej proportions estimation', 'Proportions'], 
        width=barWidth, color='#abdda4', 
        yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label=f"EMmej") # label=f"EMmej $\hat{'P'}$"

# Set plot parameters
plt.xticks([r + barWidth for r in range(len(third_plot_mech))], third_plot_mech, rotation = 45)
plt.ylabel(f'SD-MMEJ Proportion',fontsize=fs)
plt.ylim(bot, top)
plt.xlabel('WT',fontsize=fs)
plt.title(f'R0-WT (Deletions N: {WT_del_n}, Insertions N: {WT_ins_n})',fontsize=7)


plt_df = mutants_res_df.loc[((mutants_res_df['mechanism'].isin(third_plot_mech)) & 
            (mutants_res_df['Mutant'] == 'R0_POLQ')), :]
plt.subplot(1, 3, 2)
sdmmej_ci = plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'ci'].to_list()
plt.bar(r1, plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'Proportions'], width=barWidth, color='#ef8a62', 
    yerr = [[i[0] for i in sdmmej_ci], [i[1] for i in sdmmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label='sdmmej')

# Creating the bars for non EMmej
EMmej_ci = plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'ci'].to_list()
plt.bar(r2, plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'Proportions'], width=barWidth, color='#67a9cf', 
        yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label='Patterns') # label='EMmej'

# Creating the bars for non EMmej proportions estimation
EMmej_ci = plt_df.loc[((plt_df['Algorithm'] == 'EMmej proportions estimation') & (plt_df['indel_type']=='Insertion')), 'ci'].to_list()
plt.bar(r3, plt_df.loc[plt_df['Algorithm'] == 'EMmej proportions estimation', 'Proportions'], 
        width=barWidth, color='#abdda4', 
        yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label=f"EMmej") # label=f"EMmej $\hat{'P'}$"

# Set plot parameters
plt.xticks([r + barWidth for r in range(len(third_plot_mech))], third_plot_mech, rotation = 45)
plt.ylabel(f'')
plt.ylim(bot, top)
plt.xlabel('PolQ',fontsize=fs)
plt.title(f'R0-POLQ (Deletions N: {POLQ_del_n}, Insertions N: {POLQ_ins_n})',fontsize=7)



plt_df = mutants_res_df.loc[((mutants_res_df['mechanism'].isin(third_plot_mech)) & 
            (mutants_res_df['Mutant'] == 'R0_Lig4')), :]
print(plt_df)
# print(plt_df.loc[plt_df['mechanism'] == 'SD-MMEJ', :])
# print('---------------------------------')
plt.subplot(1, 3, 3)
sdmmej_ci = plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'ci'].to_list()
plt.bar(r1, plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'Proportions'],
    width=barWidth, color='#ef8a62', 
    yerr = [[i[0] for i in sdmmej_ci], [i[1] for i in sdmmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label='sdmmej')

# Creating the bars for non EMmej
EMmej_ci = plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'ci'].to_list()
plt.bar(r2, plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'Proportions'],
        width=barWidth, color='#67a9cf', 
        yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label='Patterns') # label='EMmej'

# Creating the bars for non EMmej proportions estimation
EMmej_ci = plt_df.loc[((plt_df['Algorithm'] == 'EMmej proportions estimation') & (plt_df['indel_type']=='Insertion')), 'ci'].to_list()
plt.bar(r3, plt_df.loc[plt_df['Algorithm'] == 'EMmej proportions estimation', 'Proportions'], 
        width=barWidth, color='#abdda4', 
        yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label=f"EMmej") # label=f"EMmej $\hat{'P'}$"

# Set plot parameters
plt.xticks([r + barWidth for r in range(len(third_plot_mech))], third_plot_mech, rotation = 45)
plt.legend(bbox_to_anchor=(1.85, 1.05))
plt.ylabel(f'')
plt.ylim(bot, top)
plt.xlabel('Lig4',fontsize=fs)
plt.title(f'R0-Lig4 (Deletions N: {Lig_4_del_n}, Insertions N: {Lig_4_ins_n})',fontsize=7)

fig = ax.get_figure()
fig.savefig(f"{path_to_save_plots}/mutant_comparispn.svg", bbox_inches='tight')

"""
Fourth plot:
"""

plt_df = mutants_res_df.loc[((mutants_res_df['mechanism'].isin(SDMMEJ_zoomin_mech)) & 
            (mutants_res_df['Mutant'] == 'R0_WT')), :]
# set x axis positions
r1 = np.arange(len(SDMMEJ_zoomin_mech))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]

# Creating the bars for sdmmej
# plt.style.use('seaborn-darkgrid')
plt.style.use('seaborn-ticks')
ax = figure(num=None, figsize=(16, 6))
ax.suptitle(f"EMmej and sdmmej mutant comparison, SDMMEJ stratified to sub-mechanisms", fontsize=16)
plt.subplot(1, 3, 1)
sdmmej_ci = plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'ci'].to_list()
plt.bar(r1, plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'Proportions'], width=barWidth, color='#ef8a62', 
        yerr = [[i[0] for i in sdmmej_ci], [i[1] for i in sdmmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label='sdmmej')

# Creating the bars for non EMmej
EMmej_ci = plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'ci'].to_list()
plt.bar(r2, plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'Proportions'], width=barWidth, color='#67a9cf', 
        yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label='Patterns') # label='EMmej'

# Creating the bars for non EMmej proportions estimation
EMmej_ci = plt_df.loc[((plt_df['Algorithm'] == 'EMmej proportions estimation') & (plt_df['indel_type']=='Insertion')), 'ci'].to_list()
plt.bar(r3, plt_df.loc[plt_df['Algorithm'] == 'EMmej proportions estimation', 'Proportions'], width=barWidth, color='#abdda4', 
        edgecolor='black',linewidth=linewidth ,
        yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
        error_kw=error_kw, label=f"EMmej") # label=f"EMmej $\hat{'P'}$"

# Set plot parameters
plt.xticks([r + barWidth for r in range(len(SDMMEJ_zoomin_mech))], SDMMEJ_zoomin_mech, rotation = 45)
plt.ylabel(f'Proportions (normalized to WT SD-MMEJ)',fontsize=12)
plt.ylim(bot, top)
plt.xlabel('WT',fontsize=fs)
plt.title(f'R0-WT (Deletions N: {WT_del_n}, Insertions N: {WT_ins_n})')


plt_df = mutants_res_df.loc[((mutants_res_df['mechanism'].isin(SDMMEJ_zoomin_mech)) & 
            (mutants_res_df['Mutant'] == 'R0_POLQ')), :]
plt.subplot(1, 3, 2)
sdmmej_ci = plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'ci'].to_list()
plt.bar(r1, plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'Proportions'], width=barWidth, color='#ef8a62', 
    yerr = [[i[0] for i in sdmmej_ci], [i[1] for i in sdmmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label='sdmmej')

# Creating the bars for non EMmej
EMmej_ci = plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'ci'].to_list()
plt.bar(r2, plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'Proportions'], width=barWidth, color='#67a9cf', 
        yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label='Patterns') # label='EMmej'

# Creating the bars for non EMmej proportions estimation
EMmej_ci = plt_df.loc[((plt_df['Algorithm'] == 'EMmej proportions estimation') & (plt_df['indel_type']=='Insertion')), 'ci'].to_list()
plt.bar(r3, plt_df.loc[plt_df['Algorithm'] == 'EMmej proportions estimation', 'Proportions'], width=barWidth, color='#abdda4', 
        edgecolor='black',linewidth=linewidth ,
        yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
        error_kw=error_kw, label=f"EMmej") # label=f"EMmej $\hat{'P'}$"

# Set plot parameters
plt.xticks([r + barWidth for r in range(len(SDMMEJ_zoomin_mech))], SDMMEJ_zoomin_mech, rotation = 45)
plt.ylabel(f'Proportions (normalized to POLQ SD-MMEJ)',fontsize=12)
plt.ylim(bot, top)
plt.xlabel('PolQ',fontsize=fs)
plt.title(f'R0-POLQ (Deletions N: {POLQ_del_n}, Insertions N: {POLQ_ins_n})')

plt_df = mutants_res_df.loc[((mutants_res_df['mechanism'].isin(SDMMEJ_zoomin_mech)) & 
            (mutants_res_df['Mutant'] == 'R0_Lig4')), :]

plt.subplot(1, 3, 3)
sdmmej_ci = plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'ci'].to_list()
plt.bar(r1, plt_df.loc[plt_df['Algorithm'] == 'sdmmej', 'Proportions'], width=barWidth, color='#ef8a62', 
    yerr = [[i[0] for i in sdmmej_ci], [i[1] for i in sdmmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label='sdmmej')

# Creating the bars for non EMmej
EMmej_ci = plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'ci'].to_list()
plt.bar(r2, plt_df.loc[plt_df['Algorithm'] == 'EMmej', 'Proportions'], width=barWidth, color='#67a9cf', 
        yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
        edgecolor='black',linewidth=linewidth ,
        error_kw=error_kw, label='Patterns') # label='EMmej'

# Creating the bars for non EMmej proportions estimation
EMmej_ci = plt_df.loc[((plt_df['Algorithm'] == 'EMmej proportions estimation') & (plt_df['indel_type']=='Insertion')), 'ci'].to_list()
plt.bar(r3, plt_df.loc[plt_df['Algorithm'] == 'EMmej proportions estimation', 'Proportions'], width=barWidth, color='#abdda4', 
        edgecolor='black',linewidth=linewidth ,
        yerr = [[i[0] for i in EMmej_ci], [i[1] for i in EMmej_ci]], 
        error_kw=error_kw, label=f"EMmej") # label=f"EMmej $\hat{'P'}$"

# Set plot parameters
plt.xticks([r + barWidth for r in range(len(SDMMEJ_zoomin_mech))], SDMMEJ_zoomin_mech, rotation = 45)
plt.legend(bbox_to_anchor=(1.35, 1.05))
plt.ylabel(f'Proportions (normalized to Lig4 SD-MMEJ)',fontsize=12)
plt.ylim(bot, top)
plt.xlabel('Lig4', fontsize=fs)
plt.title(f'R0-Lig4 (Deletions N: {Lig_4_del_n}, Insertions N: {Lig_4_ins_n})')

fig = ax.get_figure()
fig.savefig(f"{path_to_save_plots}/mutant_comparison_SDMMEJ_stratification_to_single_mechanisms.svg", bbox_inches='tight')

