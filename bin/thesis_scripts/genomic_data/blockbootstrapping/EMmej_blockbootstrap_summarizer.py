"""
This script summerize the results of multiple EMmej
runs into one dataframe. Usful to summerize bootstrap/batching
runs
"""

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
all_args.add_argument("-i", "--bootstraps_res", required=True,
   help="path to a folder that contains bootstrap inputs (string)")
all_args.add_argument("-sil", "--stratify_by_indel_length", 
            required=False, action='store_true',default=False, 
      help="a flag to turn on stratification by indel length")

args = vars(all_args.parse_args())
date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))


# def get_lk_proportions(df: pd.DataFrame):
#       """
#       A function that calculates the proportions of
#       mehcnaisms relative to their indel type as well
#       as relative to the whole dataset by summing
#       and normalizing lk
#       """
#       # deletions
#       del_mmej_lk_sum = df['del_mmej_lk'].sum()
#       nhej_lk_sum = df['NHEJ_lk'].sum()
#       # deletions proportions out of deletions only
#       del_mmej_del_prop = del_mmej_lk_sum/(del_mmej_lk_sum+nhej_lk_sum)
#       nhej_del_prop = nhej_lk_sum/(del_mmej_lk_sum+nhej_lk_sum)

#       # insertions proportions
#       SD_snap_back_lk_sum = df['SD_snap_back_lk'].sum()
#       SD_loop_out_lk_sum = df['SD_loop_out_lk'].sum()
#       unclassified_ins_lk_sum = df['unclassified_ins_lk'].sum()
#       # insertions proportions out of insertions only
#       SD_snap_back_lk_ins_prop = SD_snap_back_lk_sum/(SD_snap_back_lk_sum +
#                                           SD_loop_out_lk_sum + unclassified_ins_lk_sum)
#       SD_loop_out_lk_ins_prop = SD_loop_out_lk_sum/(SD_snap_back_lk_sum +
#                                           SD_loop_out_lk_sum + unclassified_ins_lk_sum)
#       unclassified_ins_lk_ins_prop = unclassified_ins_lk_sum/(SD_snap_back_lk_sum +
#                                           SD_loop_out_lk_sum + unclassified_ins_lk_sum)
      
#       tot_lk = (del_mmej_lk_sum + nhej_lk_sum +
#                   SD_snap_back_lk_sum + SD_loop_out_lk_sum +
#                   unclassified_ins_lk_sum)
#       # deletions proportions out of deletions all indels
#       del_mmej_lk_tot_prop = del_mmej_lk_sum/tot_lk
#       nhej_lk_tot_prop = nhej_lk_sum/tot_lk
#       # insertions proportions out of deletions all indels
#       SD_snap_back_lk_tot_prop = SD_snap_back_lk_sum/tot_lk
#       SD_loop_out_lk_tot_pop = SD_loop_out_lk_sum/tot_lk
#       unclassified_ins_lk_tot_pop = unclassified_ins_lk_sum/tot_lk
#       return (del_mmej_del_prop, nhej_del_prop, 
#             SD_snap_back_lk_ins_prop, SD_loop_out_lk_ins_prop,
#             unclassified_ins_lk_ins_prop, del_mmej_lk_tot_prop,
#             nhej_lk_tot_prop, SD_snap_back_lk_tot_prop, 
#             SD_loop_out_lk_tot_pop, unclassified_ins_lk_tot_pop)



dir_ls = glob.glob(f"{args['bootstraps_res']}/EM_res/*")
# print(dir_ls)
"""
Summerizing the EM bootstraping results from:
EMmej_blockbootstrap.py
For deletions
"""
EM_estimates = pd.Series([f for f in dir_ls if re.match(r'(.+?).tsv', f)])
EM_estimates = EM_estimates[~EM_estimates.str.contains('bin')]
EM_bootstrap_estimates = pd.DataFrame(columns=['NHEJ', 'MMEJ', 'boot_n'])
for r in EM_estimates:
      boot_res = pd.read_csv(r, sep='\t')
      boot_res['boot_n'] = re.search(r'bootnumber_(.[0-9]*)', r).group(1)
      EM_bootstrap_estimates = pd.concat([EM_bootstrap_estimates, boot_res])

EM_bootstrap_estimates.reset_index(drop=True, inplace=True)
# print(EM_bootstrap_estimates['boot_n'].to_list())
# for i in range(101):
#       if i not in EM_bootstrap_estimates['boot_n'].to_list(): print(i)
"""
Summerizing the EM bootstraping results from:
EMmej_blockbootstrap.py
For insertions
"""
ins_dir_ls = glob.glob(f"{args['bootstraps_res']}/insertions_proportions/*")
ins_estimates = pd.Series([f for f in ins_dir_ls if re.match(r'(.+?).tsv', f)])
ins_estimates = ins_estimates[~ins_estimates.str.contains('bin')]
ins_bootstrap_estimates = pd.DataFrame(columns=['snap', 'loop', 'SDMMEJ'])
for r in ins_estimates:
      boot_res = pd.read_csv(r, sep='\t')
      boot_res['boot_n'] = re.search(r'bootnumber_(.[0-9]*)', r).group(1)
      ins_bootstrap_estimates = pd.concat([ins_bootstrap_estimates, boot_res])
print(ins_bootstrap_estimates)
ins_bootstrap_estimates.reset_index(drop=True, inplace=True)
EMmej_bootstraps = EM_bootstrap_estimates.merge(ins_bootstrap_estimates, left_on='boot_n', right_on='boot_n')
EMmej_bootstraps.to_csv(f"{args['bootstraps_res']}/EMmej_overall_estimations.tsv",sep='\t', index=False)





if args['stratify_by_indel_length']:
      print('====================')
      """
      Summerizing the bootstraping results from:
      EMmej_blockbootstrap.py
      For different indel lengths
      Deletions
      """
      EM_bootstrap_estimates = pd.DataFrame(
            columns=['NHEJ', 'MMEJ', 'snap', 'loop', 'SDMMEJ','unclassified_ins', 'indel_len_bin','boot_n'])
      if args['stratify_by_indel_length']:
            EM_estimates = [f for f in dir_ls if re.match(r'(.+?)_bin_(.+?).tsv', f)]
            for idx,r in enumerate(EM_estimates):
                  boot_n = re.search(r'bootnumber_(.[0-9]*)', r).group(1)
                  indel_len_bin = re.search(rf'bootnumber_{boot_n}_bin_(.+?).tsv', r).group(1)
                  boot_res = pd.read_csv(r, sep='\t')
                  
                  if ':' in indel_len_bin: indel_len_bin_lable = f"[{indel_len_bin}]"
                  else: indel_len_bin_lable = indel_len_bin
                  
                  EM_bootstrap_estimates.loc[idx,'NHEJ'] = boot_res.loc[0,'NHEJ']
                  EM_bootstrap_estimates.loc[idx,'MMEJ'] = boot_res.loc[0,'MMEJ']
                  EM_bootstrap_estimates.loc[idx,'indel_len_bin'] = indel_len_bin_lable
                  EM_bootstrap_estimates.loc[idx,'boot_n'] = boot_n
            
            EM_estimates = [f for f in ins_dir_ls if re.match(r'(.+?)_bin_(.+?).tsv', f)]
            for idx,r in enumerate(EM_estimates):
                  boot_n = re.search(r'bootnumber_(.[0-9]*)_bin', r).group(1)
                  indel_len_bin = re.search(rf'bootnumber_{boot_n}_bin_(.+?)_insertions_proportions.tsv', r).group(1)
                  boot_res = pd.read_csv(r, sep='\t')
                  
                  if ':' in indel_len_bin: indel_len_bin_lable = f"[{indel_len_bin}]"
                  else: indel_len_bin_lable = indel_len_bin
                  
                  EM_bootstrap_estimates.loc[idx,'snap'] = boot_res.loc[0,'snap']
                  EM_bootstrap_estimates.loc[idx,'loop'] = boot_res.loc[0,'loop']
                  EM_bootstrap_estimates.loc[idx,'SDMMEJ'] = boot_res.loc[0,'SDMMEJ']
                  EM_bootstrap_estimates.loc[idx,'unclassified_ins'] = boot_res.loc[0,'unclassified_ins']
                  EM_bootstrap_estimates.loc[idx,'indel_len_bin'] = indel_len_bin_lable
                  EM_bootstrap_estimates.loc[idx,'boot_n'] = boot_n

            EM_bootstrap_estimates.to_csv(f"{args['bootstraps_res']}/EMmej_estimations_by_indel_length_bins.tsv", 
                                    sep='\t', index=False)



# """
# Summerizing the bootstraping results from:
# EMmej_blockbootstrap.py
# For insertions
# """
# bootstraps = pd.Series(
#       glob.glob(f"{args['bootstraps_res']}/bootstrap_dfs/*"))
# bootstraps = bootstraps[~bootstraps.str.contains('bin')]
# # print(bootstraps)
# EMmej_bootstraps = pd.DataFrame(
#                   columns=['MMEJ_EM_theta', 'NHEJ_EM_theta','SD_snap_back_prop', 
#                   'SD_loop_out_prop', 'unclassified_ins_prop',  'boot_n'], # 'SDMMEJ_prop',
#                   index=[i for i in range(len(bootstraps))])
# for idx, r in enumerate(bootstraps):
#       boot_res = pd.read_csv(r, sep='\t')
#       # boot_res.drop(columns=['index'], inplace=True)
#       ins_n = boot_res.loc[(boot_res['ancestral_indel'].str.len() < boot_res['derived_indel'].str.len()), :].shape[0]
#       del_n = boot_res.loc[(boot_res['ancestral_indel'].str.len() > boot_res['derived_indel'].str.len()), :].shape[0]
#       boot_n = re.search(r'bootnumber_(.[0-9]*)', r).group(1)
#       boot_res['boot_n'] = boot_n
#       boot_res['snap_cutoff'] = 0.001 / boot_res['variant_id_N']
#       boot_res['loop_cutoff'] = 0.0001 / boot_res['variant_id_N']
#       boot_res['unclassified_cutoff'] = 0.001 / boot_res['variant_id_N']
#       boot_res['SDMMEJ_cutoff'] = 0.001 / boot_res['variant_id_N']
      
#       EMmej_bootstraps.loc[idx, 'SD_snap_back_prop'] = (boot_res.loc[
#             (boot_res['SD_snap_back_p_val'] <= boot_res['snap_cutoff']), 
#             'variant_id'].unique().shape[0] / boot_res['variant_id'].unique().shape[0])
#       EMmej_bootstraps.loc[idx, 'SD_loop_out_prop'] = (boot_res.loc[
#             (boot_res['SD_loop_out_p_val'] <= boot_res['loop_cutoff']), 
#             'variant_id'].unique().shape[0] / boot_res['variant_id'].unique().shape[0])
#       # EMmej_bootstraps.loc[idx, 'SDMMEJ_prop'] = (boot_res.loc[
#       #       (boot_res['SDMMEJ_lk'] <= boot_res['SDMMEJ_cutoff']), 
#       #       'variant_id'].unique().shape[0] / boot_res['variant_id'].unique().shape[0])
#       EMmej_bootstraps.loc[idx, 'unclassified_ins_prop'] = (boot_res.loc[
#             (boot_res['unclassified_ins_p_val'] <= boot_res['unclassified_cutoff']), 
#             'variant_id'].unique().shape[0] / boot_res['variant_id'].unique().shape[0])
      
      
      
#       EMmej_bootstraps.loc[idx, 'boot_n'] = re.search(r'bootnumber_(.[0-9]*)', r).group(1)
#       # EMmej_bootstraps.loc[idx, 'SD_snap_back_lk'] = boot_res['SD_snap_back_lk'].sum()/ins_n
#       # EMmej_bootstraps.loc[idx, 'SD_loop_out_lk'] = boot_res['SD_loop_out_lk'].sum()/ins_n
#       # EMmej_bootstraps.loc[idx, 'unclassified_ins_lk'] = boot_res['unclassified_ins_lk'].sum()/ins_n
#       # EMmej_bootstraps.loc[idx, 'SDMMEJ_lk'] = boot_res['SDMMEJ_lk'].sum()/ins_n
      
#       EMmej_bootstraps.loc[idx, 'MMEJ_EM_theta'] = float(EM_bootstrap_estimates.loc[
#                               EM_bootstrap_estimates['boot_n'] == boot_n, 'MMEJ'])
#       EMmej_bootstraps.loc[idx, 'NHEJ_EM_theta'] = float(EM_bootstrap_estimates.loc[
#                               EM_bootstrap_estimates['boot_n'] == boot_n, 'NHEJ'])
            
# EMmej_bootstraps.reset_index(drop=True, inplace=True)
# EMmej_bootstraps.to_csv(f"{args['bootstraps_res']}/EMmej_overall_estimations.tsv",sep='\t', index=False)
# print(EMmej_bootstraps)

# if args['stratify_by_indel_length']:
#       """
#       Summerizing the bootstraping results from:
#       EMmej_blockbootstrap.py
#       For different indel lengths
#       Deletions
#       """
#       EM_bootstrap_estimates = pd.DataFrame(columns=['NHEJ', 'MMEJ', 'indel_len_bin','boot_n'])
#       if args['stratify_by_indel_length']:
#             EM_estimates = [f for f in dir_ls if re.match(r'(.+?)_bin_(.+?).tsv', f)]
#             for idx,r in enumerate(EM_estimates):
#                   boot_n = re.search(r'bootnumber_(.[0-9]*)', r).group(1)
#                   indel_len_bin = re.search(rf'bootnumber_{boot_n}_bin_(.+?).tsv', r).group(1)
#                   boot_res = pd.read_csv(r, sep='\t')
                  
#                   if ':' in indel_len_bin: indel_len_bin_lable = f"[{indel_len_bin}]"
#                   else: indel_len_bin_lable = indel_len_bin
                  
#                   EM_bootstrap_estimates.loc[idx,'NHEJ'] = boot_res.loc[0,'NHEJ']
#                   EM_bootstrap_estimates.loc[idx,'MMEJ'] = boot_res.loc[0,'MMEJ']
#                   EM_bootstrap_estimates.loc[idx,'indel_len_bin'] = indel_len_bin_lable
#                   EM_bootstrap_estimates.loc[idx,'boot_n'] = boot_n

#       """
#       Summerizing the bootstraping results from:
#       EMmej_blockbootstrap.py
#       For different indel lengths
#       Insertions
#       """
#       bootstraps = pd.Series(glob.glob(f"{args['bootstraps_res']}/bootstrap_dfs/*"))
#       bootstraps = bootstraps[bootstraps.str.contains('bin')]
#       EMmej_bootstraps = pd.DataFrame(columns=['MMEJ_EM_theta', 'NHEJ_EM_theta','SD_snap_back_lk', 
#                         'SD_loop_out_lk', 'unclassified_ins_lk', 'indel_len_bin','boot_n'],
#                         index=[i for i in range(len(bootstraps))])

#       for idx, r in enumerate(bootstraps):
#             boot_n = re.search(r'bootnumber_(.[0-9]*)', r).group(1)
#             # try:
#             indel_len_bin = re.search(rf'bootnumber_{boot_n}_bin_(.+?).tsv', r).group(1)
#             if ':' in indel_len_bin: indel_len_bin_lable = f"[{indel_len_bin}]"
#             else: indel_len_bin_lable = indel_len_bin
#             # except: 
#             #       continue
#             boot_res = pd.read_csv(r, sep='\t')
#             ins_n = boot_res.loc[
#                   (boot_res['ancestral_indel'].str.len() < 
#                   boot_res['derived_indel'].str.len()), :].shape[0]
#             del_n = boot_res.loc[
#                   (boot_res['ancestral_indel'].str.len() >
#                   boot_res['derived_indel'].str.len()), :].shape[0]
#             EMmej_bootstraps.loc[idx, 'boot_n'] = boot_n
#             EMmej_bootstraps.loc[idx, 'indel_len_bin'] = indel_len_bin_lable

#             if '-' in indel_len_bin: # Set EM estimetion for deletions
#                   EMmej_bootstraps.loc[idx, 'MMEJ_EM_theta'] = float(EM_bootstrap_estimates.loc[
#                         ((EM_bootstrap_estimates['boot_n'] == boot_n) & 
#                         (EM_bootstrap_estimates['indel_len_bin'] == indel_len_bin_lable)), 'MMEJ'])
#                   EMmej_bootstraps.loc[idx, 'NHEJ_EM_theta'] = float(EM_bootstrap_estimates.loc[
#                         ((EM_bootstrap_estimates['boot_n'] == boot_n) & 
#                         (EM_bootstrap_estimates['indel_len_bin'] == indel_len_bin_lable)), 'NHEJ'])
            
#             else: # Set likelihood based estimetion for insertions
#                   EMmej_bootstraps.loc[idx, ['MMEJ_EM_theta','NHEJ_EM_theta']] = np.nan, np.nan
#                   boot_res['boot_n'] = boot_n
#                   EMmej_bootstraps.loc[idx, 'SD_snap_back_lk'] = boot_res['SD_snap_back_lk'].sum()/ins_n
#                   EMmej_bootstraps.loc[idx, 'SD_loop_out_lk'] = boot_res['SD_loop_out_lk'].sum()/ins_n
#                   EMmej_bootstraps.loc[idx, 'unclassified_ins_lk'] = boot_res['unclassified_ins_lk'].sum()/ins_n
#                   EMmej_bootstraps.loc[idx, 'SDMMEJ_lk'] = boot_res['SDMMEJ_lk'].sum()/ins_n
#       EMmej_bootstraps.reset_index(drop=True, inplace=True)
#       EMmej_bootstraps.to_csv(f"{args['bootstraps_res']}/EMmej_estimations_by_indel_length_bins.tsv", 
#                                     sep='\t', index=False)


# Command example:
# python3 EMmej_blockbootstrap_summarizer.py -i /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/blockbootstraping/GDL_indel_EMmej_format_non_repeated_regions_exon_no_MT/20221023_blockbootstrap_output -sil
