

import argparse
import datetime
import os
import glob
import re

import pandas as pd
import numpy as np


# Construct an argument parser
all_args = argparse.ArgumentParser()
all_args.add_argument("-v", "--vcf", required=True,
   help="path to input vcf (string)")
all_args.add_argument("-o", "--outputpath", required=True,
      help="path to output file (string)")

args = vars(all_args.parse_args())
date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))


# Creating a list of all folders that contains blockbootstraping runs
dir_ls = glob.glob(f"{args['vcf']}/*")
# # Creating a dataframe for the summary
# cols = ['#CHR', 'blockN', 'bootN', 'MMEJ_theta', 'NHEJ_theta', 
#                         'del_mmej_del_prop', 'nhej_del_prop', 
#             'SD_snap_back_lk_ins_prop', 'SD_loop_out_lk_ins_prop',
#             'unclassified_ins_lk_ins_prop', 'del_mmej_lk_tot_prop',
#             'nhej_lk_tot_prop', 'SD_snap_back_lk_tot_prop', 
#             'SD_loop_out_lk_tot_pop', 'unclassified_ins_lk_tot_pop']
EMmej_logs = pd.DataFrame(index=range(len(dir_ls))) # columns=cols,
# Empty list to get names of runs that havn't finished yet
undone_runs = []

for i, fol in enumerate(dir_ls):
    output_ls = glob.glob(f"{fol}/output/*")
    EM_log = [f for f in os.listdir(f"{fol}/output_files/") if re.match(r'[0-9]*_EM_log.tsv', f)]
    if len(EM_log) > 0:
        MM_output = [f for f in os.listdir(f"{fol}/output/") if re.match(r'[0-9]*MM_output.tsv', f)]



