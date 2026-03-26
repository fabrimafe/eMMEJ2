# 15.05.2022
"""
This script takes the realigned indels (made by parsevcf_v2.py, 
then assigned to repair mechanisms by 20220504_MMEJ_detection_pipline.py)
and aggregate them into one indel per original location.
"""
# importing libreries
import re
import sys
sys.path.append('/home/labs/alevy/guyta/guy_master_project/src')
import io
import datetime
from time import gmtime, strftime, localtime
import argparse

import pandas as pd
import numpy as np

from data_exploration_util import *

# setting pandas display options
pd.options.display.max_colwidth = 2200
pd.set_option("display.max_columns", None)

mem_reporter()
date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
local_h = strftime("%H%M", localtime())


# laod data
# Construct an argument parser
all_args = argparse.ArgumentParser()

# Add arguments to the parser
all_args.add_argument("-p", "--path_to_data", required=True,
   help="path to input data")
all_args.add_argument("-f", "--filename", required=True,
   help="path to input data")
all_args.add_argument("-o", "--outputfile", required=True,
      help="output file in vcf-like format")
all_args.add_argument("-a", "--aggby", required=True,
      help="aggregation method (max,mean)")

args = vars(all_args.parse_args())

def get_max_prob(data: pd.DataFrame):
   """
   This function takes a dataframe that contains only one original_pos
   and extract the indel/indels that have thiere 
   max(prob among mechanisms) == max(max(prob among mechanisms) among indels)
   Args:
      data: pd.Dataframe -> A dataframe that contains only one original_pos
   Return:
      max_prob_df: pd.DataFrame -> A vector containing the desiered indel.
   """
   data['max_prob'] = data.loc[:, ['del_mmej_prob', 
                  'trans_mmej_prob', 'SD_snap_back_prob', 'SD_loop_out_prob',
                   'NHEJ_prob', 'unclassified_ins_prob']].max(axis=1)
   max_prob_df = data.loc[data['max_prob'] == data['max_prob'].max(), :]
   
   ### ASK FABRIZIO: in some cases there are more then one indels that has the same probability
   ### That == max probability among mechanisms. Which one of them should I take? All of them?
   # if max_prob_df.shape[0] > 1:
   #    max_prob_df = max_prob_df.loc[max_prob_df['POS'] == max_prob_df['original_pos'], :]
   return max_prob_df


def get_avg_prob(data: pd.DataFrame):
   """
   This function takes a dataframe that contains only one original_pos
   and averege thier probabilities, then assign the avg probabilities to
   all indels. 
   Args:
      data: pd.Dataframe -> A dataframe that contains only one original_pos
   Return:
      data: pd.DataFrame -> A dataframe with all indels from the input
         data but with avg probabilities instead of the original ones.
   """
   mechanisms_cols = ['del_mmej_prob', 
                  'trans_mmej_prob', 'SD_snap_back_prob', 'SD_loop_out_prob',
                   'NHEJ_prob', 'unclassified_ins_prob']
   avg_prob = data.loc[:,mechanisms_cols].mean().to_list()
   data.loc[:,mechanisms_cols] = avg_prob
   return data


def get_nalignments_cutoff(data: pd.DataFrame, cutoff: float) -> float:
   """
   This set a probability cutoff for a group of realignments
   that correspond to the same indel based on a given cutoff
   and the number of realigmnents.
   Args:
      data (pd.DataFrame):
      cutoff (float): 
   returns:
      new_cutoff (float): 
   """
   pass

def aggindels(data: pd.DataFrame,pos: int):
   """
   A function that find the indel with the highest certeainty to be
   one repair mechanism using get_max_prob().
   Args:
      data: pd.DataFrame -> data to operate on
      pos: int -> a unique original pos
   Return:
      max_prob_df: pd.DataFrame -> A vector containing the desiered indel.
   """
   pos_df = data.loc[(data['original_pos'] == pos),: ].copy()
   if len(pos_df['CHR'].unique())<2:
      if args["aggby"] == 'max':
         return get_max_prob(data=pos_df)
      elif args["aggby"] == 'mean':
         return get_avg_prob(data=pos_df)
   else:
      if args["aggby"] == 'max':
         return pd.concat(
         [get_max_prob(data=pos_df.loc[(pos_df['CHR']==Chr), :].copy()) for Chr in pos_df['CHR'].unique()])
      elif args["aggby"] == 'mean':
         return pd.concat(
         [get_avg_prob(data=pos_df.loc[(pos_df['CHR']==Chr), :].copy()) for Chr in pos_df['CHR'].unique()])

Dtypes = {'CHR': str,'POS': int, 'original_pos': float, 'REF': str, 'ALT': str,
 'del_mmej_prob': float ,'trans_mmej_prob': float, 'SD_snap_back_prob': float, 
 'SD_loop_out_prob': float, 'NHEJ_prob': float, 'unclassified_ins_prob': float}

realigned_df = pd.read_csv(f'{args["path_to_data"]}/{args["filename"]}', sep='\t', dtype=Dtypes, header=0)
realigned_df['original_pos'] = realigned_df['original_pos'].astype('str').str.extract('(^\d*)').astype('int')

agg_df = pd.DataFrame(columns=realigned_df.columns)
for p in realigned_df['original_pos'].unique():
   agg_df = pd.concat([agg_df, aggindels(data=realigned_df, pos=p)])

print(agg_df.head(30))
agg_df.to_csv(f'{args["outputfile"]}/{args["aggby"]}_aggregated_{args["filename"]}', sep='\t', index=False)
mem_reporter()

# Example: python3 20220515_realigned_indels_aggregation.py -p /home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/20220512/2312_20220512_main_pipline_output_test9/test9_realinged/ -f 2312_20220512_test9_realined_classified_repair_mechanism_prob.csv -o /home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/20220517/ -a mean