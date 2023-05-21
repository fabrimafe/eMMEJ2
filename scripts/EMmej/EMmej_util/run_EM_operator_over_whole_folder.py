#



# importing relevant libreries
import argparse
from asyncio import subprocess
import datetime
import os
import glob
from subprocess import Popen, PIPE, STDOUT

import re

import pandas as pd
import numpy as np

# Construct an argument parser
all_args = argparse.ArgumentParser()
all_args.add_argument("-v", "--bootstrap_pth", required=True,
   help="path to input vcf (string)")
all_args.add_argument("-r", "--ref", required=False,
      help="path to indexed reference genome .fa file (string)")
all_args.add_argument("-b", "--bootnumber", 
      required=False, type=int, default=100,
      help="path to output file (string)")
all_args.add_argument("-m", "--mem", 
      required=False, type=int, default=1500,
      help="memory required")


args = vars(all_args.parse_args())

# def intersection(lst1, lst2):
#     lst3 = [value for value in lst1 if value in lst2]
#     return lst3
def parallel(lst1, lst2):
    lst3 = [value for value in lst1 if value not in lst2]
    return lst3

# date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))

# Creating output folders
# dir = os.path.join(f"{args['outputpath']}/{date}_blockbootstrap_output")
# if not os.path.exists(dir):
#     os.mkdir(dir)
#     os.mkdir(f"{dir}/bootstrap_dfs")
#     os.mkdir(f"{dir}/EM_res")
#     os.mkdir(f"{dir}/insertions_proportions")




"""
Running the EM over the bootstrap dfs
"""
dir = args['bootstrap_pth']
boot_dfs = glob.glob(f"{dir}/bootstrap_dfs/*")
EMmej_submition_to_cluster = '/EMmej_submition_to_cluster.sh'
script_path ='/home/labs/alevy/guyta/guy_master_project/scripts/EMmej'

# boot_dfs = boot_dfs[:2]




for idx,b_df in enumerate(boot_dfs):
    filename = b_df.split(sep="/")[-1]
    # Estimating deletions proportions
    bsub = ['bsub', '-q', 'new-short', '-m', 'public_hosts', '-R', f"rusage[mem={args['mem']}]",
     '-e', '/home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt' ,
     '-o' ,'/home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt']
    Em_command = ['python3', f"{script_path}/EM_operator.py", '-v', b_df, 
      '-o', f"{dir}/EM_res/{filename[:-4]}.tsv", '-r',  args['ref'], '-e', '12']
    EM_submition_to_cluster_command = bsub + Em_command
    print(f'----------------------- {idx} --------------------')
    print(EM_submition_to_cluster_command)
    Popen(EM_submition_to_cluster_command, 
            stdout=PIPE, stderr=PIPE, text=True)






for idx,b_df in enumerate(boot_dfs):
    filename = b_df.split(sep="/")[-1]
    # Estimating deletions proportions
    # bsub = ['bsub', '-q', 'new-short', '-m', 'public_hosts', '-R', f"rusage[mem={args['mem']}]",
    #  '-e', '/home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt' ,
    #  '-o' ,'/home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt']
    # Em_command = ['python3', f"{script_path}/EM_operator.py", '-v', b_df, 
    #   '-o', f"{dir}/EM_res/{filename[:-4]}.tsv", '-r',  args['ref'], '-e', '12']
    # EM_submition_to_cluster_command = bsub + Em_command
    # print(f'----------------------- {idx} --------------------')
    # print(EM_submition_to_cluster_command)
    # Popen(EM_submition_to_cluster_command, 
    #         stdout=PIPE, stderr=PIPE, text=True)
    
    # Estimating Insertions proportions
    # Setting cutoffs based on simulations results
    snap_ins_cutoff = 0.005
    loop_ins_cutoff = 0.003
    # Isolate insertons:
    df = pd.read_csv(b_df, sep='\t')
    ins_df = df.loc[(df['ancestral_indel'].str.len() < df['derived_indel'].str.len()), :].copy()
    insertions_proportions = pd.DataFrame(columns=['snap', 'loop', 'SDMMEJ','unclassified_ins'])
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
    insertions_proportions.loc[0,'SDMMEJ'] = SDMMEJ_markov_proportion
    insertions_proportions.loc[0,'unclassified_ins'] = (1 - snap_markov_proportion - loop_markov_proportion - SDMMEJ_markov_proportion)
    insertions_proportions.to_csv(
        f"{dir}/insertions_proportions/{filename[:-4]}_insertions_proportions.tsv", 
                                                sep='\t', index=False)


# python3 run_EM_operator_over_whole_folder.py -v /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/20221126/Clark_GDL_EMmej_formated_non_repeated_regions_genic/20221128_blockbootstrap_output -r /home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/ref_genome/ISO1_GCF_000001215.4.fa -m 2000
# bsub -q new-short -e /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt -o /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt -m public_hosts -R "rusage[mem=4000]" "python3 run_EM_operator_over_whole_folder.py -v /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/20221126/Clark_GDL_EMmej_formated_non_repeated_regions_genic/20221128_blockbootstrap_output -r /home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/ref_genome/ISO1_GCF_000001215.4.fa -m 3500"