# python3 EM_operator.py -v 20220718_MM_output.csv -o /home/labs/alevy/guyta/guy_master_project/scripts/general_scripts/RM_wraped_software -w 15 -r /home/labs/alevy/guyta/guy_master_project/data/indels_simulations/from_fabrizio/indel_sumlations/simulations/GWHAAEV00000000.1.genome.fasta -b 10 -bs 0.7

import sys
# sys.path.append('/home/labs/alevy/guyta/guy_master_project/src')
sys.path.append('src')
from random import uniform
from io import StringIO
import datetime
from time import gmtime, strftime, localtime
import argparse
import logging
import regex as re
import pandas as pd
import numpy as np
from pysam import FastaFile
import math

from ExpectationMaximization_q import EMq
from pysam_getfasta import *

pd.options.display.max_colwidth = 3100
pd.set_option('display.max_columns', 500)

# # setting pandas display options
# pd.options.display.max_colwidth = 3500
date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
local_h = strftime("%H%M", localtime())

# Construct an argument parser
all_args = argparse.ArgumentParser()
# Add arguments to the parser
all_args.add_argument("-v", "--vcf", required=True, 
   help="path to input vcf (string)")
all_args.add_argument("-u", "--vcf2", required=False,default=False, 
   help="path to input vcf2 (string)")
all_args.add_argument("-o", "--outputfile", required=False, default="EM_output.tsv",
      help="path to output file (string)")
all_args.add_argument("-r", "--ref", required=True,
   help="path to reference genome in fasta format (string)")
all_args.add_argument("-vr", "--verbose", required=False, 
        action='store_true',default=False, 
      help="a flag to turn on verbose mode")
all_args.add_argument("-w", "--windowsize", required=False,
      default=300, type=int,
      help="size (in bp) of the context window (int)")
all_args.add_argument("-b", "--bootstrap", required=False, 
        default=0, type=int,
      help="number of bootstraps to perform (int)")
all_args.add_argument("-e", "--EMalgorithm", required=False,
    default=3, type=int,
      help="type of Expectation Maximization used. choose 1 for q+markov mixed implementation, 2 for raw p-value, 3 for conditional probabilities estimated using p-values, 5 for presence/absence of MH, others are experimental.")
all_args.add_argument("-bs", "--bootsize", required=False,
    default=100, type=int,
      help="number of variants for each bootstrap")
all_args.add_argument("-bs2", "--bootsize2", required=False,
    default=100, type=int,
      help="number of variants each bootstrap from vcf2")
all_args.add_argument("-i", "--initial_proportion_NHEJ", required=False,
    default=0.9, type=float,
      help="number of bootstraps to perform (float)")
all_args.add_argument("-c", "--convergence", required=False,
        default=0.0000001,  type=float,
            help="convergence threshold (float)")
all_args.add_argument("-l", "--logs", 
    required=False, default="EM_logs.tsv",
      help="logs file (string)")
all_args.add_argument("-pd", "--posteriordecoding", 
    required=False, default="EM_posterior_decoding.tsv",
      help="logs file (string)")

args = vars(all_args.parse_args())

# EMalgorithm=int(args['EMalgorithm'])
vcf2 = args['vcf2']
prob_NHEJ = float(args['initial_proportion_NHEJ'])
convergence_threshold=args['convergence']
outfile=args['outputfile']
logs_outputfile=args['logs']
posteriordecodingfile=args['posteriordecoding']

Dtypes = {'CHR': str, 'POS': int, 'original_pos': str, 'variant_id': str,'variant_id_N': str,
             'direction':int,'ANC': str, 
            'DER': str, 'indel_len':int, 'del_mmej_lk':float, 'del_mmej_p_val':float, 
        'SDMMEJ_lk': float, 'SD_snap_back_lk':float, 'SD_snap_back_lk_bc':float,'SD_snap_back_p_val':float, 
         'SD_loop_out_lk':float, 'SD_loop_out_lk_bc':float,'SD_loop_out_p_val':float,  
         'NHEJ_lk':float, 'NHEJ_p_val':float, 'unclassified_ins_lk':float, 
        'unclassified_ins_p_val':float,  'del_mmej_cand':str,'del_mmej_cand_len':float,'del_mmej_lk_MH2':float, 
        'del_mmej_p_val_MH2':float,
        'NHEJ_lk_MH2':float,'NHEJ_p_val_MH2':float, 'pChance_MMEJ':float,
        'pChance_MMEJ_MH2':float, 'pChance_loop':float, 'pChance_snap':float, 'motif_pos':str,
        'loop_repeat_pat':str, 'loop_repeat_pat_len':float,
        'snap_repeat_pat':str, 'snap_repeat_pat_len':float,
        'loop_motif_pos':str, 'loop_freq_small_window':float, 'loop_freq_large_window':float,
        'snap_motif_pos':str, 'snap_freq_small_window':float, 'snap_freq_large_window':float,
        'del_mmej_freq_motif_eq':float,'loop_freq_motif_eq':float ,'snap_freq_motif_eq':float}

col_names = ['CHR', 
        'POS', 'original_pos','variant_id', 'variant_id_N','direction','ANC', 'DER', 
        'del_mmej_lk', 'del_mmej_p_val', 

        'SDMMEJ_lk','SD_snap_back_lk','SD_snap_back_lk_bc','SD_snap_back_p_val', 'pChance_snap',
        'snap_repeat_pat', 'snap_repeat_pat_len','snap_freq_motif_eq',
        'SD_loop_out_lk', 'SD_loop_out_lk_bc', 'SD_loop_out_p_val','pChance_loop', 'loop_freq_motif_eq',
        'loop_repeat_pat', 'loop_repeat_pat_len',
        'NHEJ_lk','NHEJ_p_val', 

        'unclassified_ins_lk', 'unclassified_ins_p_val',
        'indel_len','del_mmej_cand','del_mmej_cand_len','del_mmej_freq_motif_eq',
        'del_mmej_lk_MH2', 'del_mmej_p_val_MH2','pChance_MMEJ',
        'pChance_MMEJ_MH2', 'NHEJ_lk_MH2', 'NHEJ_p_val_MH2', 
        'del_mmej_motif_pos', 'del_mmej_freq_small_window', 'del_mmej_freq_large_window',
        'loop_motif_pos', 'loop_freq_small_window', 'loop_freq_large_window',
        'snap_motif_pos', 'snap_freq_small_window', 'snap_freq_large_window']

path_to_data = args["vcf"]
#df = pd.read_csv(path_to_data, sep="\t", header=0, dtype=Dtypes, names=col_names) # skiprows=1,
df = pd.read_csv(path_to_data, sep="\t") #, header=0, dtype=Dtypes, names=col_names) # skiprows=1,

# take deletions only
df = df.loc[(df['ANC'].str.len() > df['DER'].str.len()), :]
fastafile=args['ref']
refFA=FastaFile(fastafile)

##########################################
### EM using indel length distribution ###
##########################################
if args['EMalgorithm']==1:
    #print(df.loc[2,])
    #print(df.loc[2,'del_mmej_motif_pos'])
    # Getting the reference genome around the indel
    print("Running EM based on indel lengths")
    df.loc[(df['direction'] == 1), 'POS'] = df.loc[(df['direction'] == 1), :].apply(lambda row: (row['POS']-row['indel_len']+1),axis=1)
    df.loc[:, f"ref_context_seq_{args['windowsize']}bp"] = df.apply(lambda x: get_ref_context(refFA=refFA, chrom=x['CHR'], indel_pos=x['POS'], indel_seq=x['DER'],context_window_size=args['windowsize']), axis=1)

    df.loc[:, f"ref_context_seq_{args['windowsize']}bp"] = df.loc[:, f"ref_context_seq_{args['windowsize']}bp"].str.upper()
    df.loc[(df['direction'] == 1), f"ref_context_seq_{args['windowsize']}bp"] = df.loc[(df['direction'] == 1), f"ref_context_seq_{args['windowsize']}bp"].apply(lambda row: row[::-1].upper())
    df['indel_len'] = df['indel_len'].abs()
    #print(df)
    EMq_obj = EMq(data=df,
                initial_theta=0.1,
                convergence_threshold=args['convergence'],
                window_size=args['windowsize']) #,
#                MM_lk='del_mmej_lk')

    # print(f'EMq final theta MMEJ = {EMq_obj.theta_a}')
    EMq_obj.log.to_csv(outfile,sep='\t', index=False)
    EMq_obj.indel_length_dist_log.to_csv(logs_outputfile,sep='\t', index=False)
    df['del_MMEJ_EMq_post_decoding'] = EMq_obj.del_MMEJ_post_decoding
    df['del_NHEJ_EMq_post_decoding'] = EMq_obj.del_NHEJ_post_decoding
    df.to_csv(posteriordecodingfile, sep='\t', index=False)
    

    if args['bootstrap']>0:
        #boot_s = int(round(boot_size*(len(df)),0))
        boot_s = args['bootsize']
        boot_df = pd.DataFrame(columns=['boot_iter', 'MMEJ_theta', 'NHEJ_theta'])
        for b in range((args['bootstrap']+1)):
            sampled_variants = df['variant_id'].sample(n=boot_s,replace=True)
            tmp_df = df.loc[df['variant_id'].isin(sampled_variants), :].copy()
            EMq_obj = EMq(data=tmp_df,
                    initial_theta=0.1,
                    convergence_threshold=args['convergence'],
                    window_size=args['windowsize'],
                    MM_lk='del_mmej_lk')
            boot_df.loc[b, 'boot_iter'] = b
            boot_df.loc[b, 'MMEJ_theta'] = EMq_obj.theta_a
            boot_df.loc[b, 'NHEJ_theta'] = (1-EMq_obj.theta_a)

        boot_df.to_csv(posteriordecodingfile,sep='\t', index=False)


########################################
##### EM using motif probabilities #####
########################################

xvcf_0=df

if (args['EMalgorithm']!=1) & (xvcf_0.shape[0] > 3):
    print("Running EM type 2")
    proportions_df = pd.DataFrame(columns=['NHEJ', 'MMEJ'], 
            index=[i for i in range((args['bootstrap']+1))])

    unique_variants_1=np.unique(np.array(xvcf_0['variant_id']))
    if vcf2:
        xvcf2_0=pd.read_csv(vcf2, sep='\t') 
        unique_variants_2=np.unique(np.array(xvcf2_0['variant_id']))

    for iboostrap in range(0,args['bootstrap']+1):
        if args['bootstrap']>0:
            print("running bootstrap")
            
            # Hanaling boostraps when the data has less datapoints then the defult bootsize
            if args['bootsize'] > len(unique_variants_1): 
                args['bootsize'] = len(unique_variants_1)
                print(f"""\n#### Data is too small for defult bootsize, 
            setting bootsize to uniq variant_id number: {len(unique_variants_1)} ####\n""")
            
            sampled_vars=np.random.choice(unique_variants_1, size=args['bootsize'],replace=False)
            xvcf=xvcf_0[xvcf_0['variant_id'].isin(sampled_vars)].copy()
            print("total number variants in file1 is: " + str(xvcf.shape[0]))
            if vcf2:
                sampled_vars=np.random.choice(unique_variants_2, size=args['bootsize2'],replace=False)
                xvcf2_b=xvcf2_0[xvcf2_0['variant_id'].isin(sampled_vars)]
                xvcf = pd.concat([xvcf, xvcf2_b], axis=0)
                print("total number variants in file1+file2 is: " + str(xvcf.shape[0]))
        else:
            xvcf=xvcf_0
        xvcf.sort_values(['POS', 'CHR'], inplace=True)
        n_variants=len(np.unique(np.array(xvcf['variant_id'])))
        print("total number variants is: " + str(n_variants) + "in vcf with:" + str(xvcf.shape))
        print("Start EM")
        prob_NHEJ = float(args['initial_proportion_NHEJ'])
        # prob_MMEJ = round(uniform(0.05, 0.95), 1)
        print(f"{prob_NHEJ}")
        if (args['EMalgorithm']==2):
            xvcf.loc[:,'prob_NHEJ']=xvcf['del_mmej_p_val'].astype('float64')
            xvcf.loc[xvcf['prob_NHEJ']<0.49,'prob_NHEJ']=xvcf.loc[xvcf['prob_NHEJ']<0.49,'prob_NHEJ']+10**(-10)
            xvcf.loc[xvcf['prob_NHEJ']>0.51,'prob_NHEJ']=xvcf.loc[xvcf['prob_NHEJ']>0.51,'prob_NHEJ']-10**(-10)
        elif (args['EMalgorithm']==3):
            xvcf.loc[:,'prob_NHEJ']=xvcf['del_mmej_p_val'].astype('float64')
            xvcf['hasMH']=-(xvcf['del_mmej_cand_len'].isnull().astype('int')-1)
            xvcf['prob_NHEJ']=xvcf['prob_NHEJ']/(xvcf['prob_NHEJ']+xvcf['hasMH']+10**(-16))
        elif (args['EMalgorithm']==4):
            xvcf.loc[:,'prob_NHEJ']=1-xvcf['del_mmej_lk'].astype('float64')
            xvcf['hasMH']=-(xvcf['del_mmej_cand_len'].isnull().astype('int')-1)
            xvcf['prob_NHEJ']=xvcf['prob_NHEJ']/(xvcf['prob_NHEJ']+xvcf['hasMH']+10**(-16))
        elif (args['EMalgorithm']==5):
            xvcf['hasMH']=1-xvcf['del_mmej_cand_len'].isnull().astype('int')
            xvcf['prob_NHEJ']=1-xvcf['hasMH']+10**(-6)
        elif (args['EMalgorithm']==6):
            #xvcf.loc[:,'prob_NHEJ']=xvcf['del_mmej_lk_MH2'].astype('float64')
            xvcf.loc[:,'prob_NHEJ']=xvcf['pChance_MMEJ'].astype('float64')
            xvcf.loc[xvcf['prob_NHEJ'].isnull(),'prob_NHEJ']=1
            xvcf['hasMH']=-(xvcf['del_mmej_cand_len'].isnull().astype('int')-1)
            xvcf['prob_NHEJ']=xvcf['prob_NHEJ']/(xvcf['prob_NHEJ']+xvcf['hasMH']+10**(-16))
        elif (args['EMalgorithm']==7):
            #xvcf.loc[:,'prob_NHEJ']=xvcf['del_mmej_lk_MH2'].astype('float64')
            xvcf.loc[:,'prob_NHEJ']=1-xvcf['pChance_MMEJ'].astype('float64')
            xvcf.loc[xvcf['prob_NHEJ'].isnull(),'prob_NHEJ']=1
        elif (args['EMalgorithm']==8):
            #xvcf.loc[:,'prob_NHEJ']=xvcf['del_mmej_lk_MH2'].astype('float64')
            xvcf.loc[:,'prob_NHEJ']=xvcf['pChance_MMEJ_MH2'].astype('float64')
            xvcf.loc[xvcf['prob_NHEJ'].isnull(),'prob_NHEJ']=1
            xvcf['hasMH']=1-xvcf['pChance_MMEJ_MH2'].isnull().astype('int')
            xvcf['prob_NHEJ']=xvcf['prob_NHEJ']/(xvcf['prob_NHEJ']+xvcf['hasMH']+10**(-16))
        elif (args['EMalgorithm']==9):
            xvcf.loc[:,'prob_NHEJ']=xvcf['del_mmej_lk'].astype('float64')+10**(-16)
            xvcf.loc[:,'prob_MMEJ']=xvcf['prob_NHEJ']/(xvcf['del_mmej_p_val'].astype('float64')+10**(-16))
            #xvcf['hasMH']=-(xvcf['del_mmej_cand_len'].isnull().astype('int')-1)
            xvcf['prob_NHEJ']=xvcf['prob_NHEJ']/(xvcf['prob_NHEJ']+xvcf['prob_MMEJ']) #+10**(-16))
        #xvcf.loc[xvcf['prob_NHEJ']<0.49,'prob_NHEJ']=xvcf.loc[xvcf['prob_NHEJ']<0.49,'prob_NHEJ']+10**(-10)
        #xvcf.loc[xvcf['prob_NHEJ']>0.51,'prob_NHEJ']=xvcf.loc[xvcf['prob_NHEJ']>0.51,'prob_NHEJ']-10**(-10)
        elif (args['EMalgorithm']==10):
            xvcf.loc[:,'prob_NHEJ']=1-xvcf['del_mmej_lk'].astype('float64')
            xvcf.loc[xvcf['prob_NHEJ']<0.5,'prob_NHEJ']=xvcf.loc[xvcf['prob_NHEJ']<0.5,'prob_NHEJ']**2
            xvcf.loc[xvcf['prob_NHEJ']>0.5,'prob_NHEJ']=1-((1-xvcf.loc[xvcf['prob_NHEJ']>0.5,'prob_NHEJ'])**2)
            xvcf['hasMH']=-(xvcf['del_mmej_cand_len'].isnull().astype('int')-1)
            xvcf['prob_NHEJ']=xvcf['prob_NHEJ']/(xvcf['prob_NHEJ']+xvcf['hasMH']+10**(-16))
        elif (args['EMalgorithm']==11):
            xvcf.loc[xvcf['del_mmej_motif_pos'].isnull(),'del_mmej_motif_pos']="-"
            xvcf['freq_motif']=xvcf['del_mmej_motif_pos'].apply(lambda row: row.count(".0"))+1.0
            xvcf['freq_motif']=xvcf['freq_motif']/(args['windowsize']-xvcf['del_mmej_cand_len']+1)
            xvcf.loc[xvcf['freq_motif'].isnull(),'freq_motif']=1
            # print(xvcf['freq_motif'])
            xvcf.loc[:,'prob_NHEJ']=(1-xvcf['del_mmej_lk'].astype('float64'))*xvcf['freq_motif']
            #xvcf.loc[xvcf['prob_NHEJ']<0.5,'prob_NHEJ']=xvcf.loc[xvcf['prob_NHEJ']<0.5,'prob_NHEJ']**2
            #xvcf.loc[xvcf['prob_NHEJ']>0.5,'prob_NHEJ']=1-((1-xvcf.loc[xvcf['prob_NHEJ']>0.5,'prob_NHEJ'])**2)
            xvcf['hasMH']=-(xvcf['del_mmej_cand_len'].isnull().astype('int')-1)
        elif (args['EMalgorithm']==12):
            xvcf['freq_motif']=xvcf['del_mmej_freq_motif_eq']
            xvcf.loc[xvcf['freq_motif'].isnull(),'freq_motif']=1
            # print(xvcf['freq_motif'])
            xvcf.loc[:,'prob_NHEJ']=(1-xvcf['del_mmej_lk'].astype('float64'))*xvcf['freq_motif']
            #xvcf.loc[xvcf['prob_NHEJ']<0.5,'prob_NHEJ']=xvcf.loc[xvcf['prob_NHEJ']<0.5,'prob_NHEJ']**2
            #xvcf.loc[xvcf['prob_NHEJ']>0.5,'prob_NHEJ']=1-((1-xvcf.loc[xvcf['prob_NHEJ']>0.5,'prob_NHEJ'])**2)
            xvcf['hasMH']=-(xvcf['del_mmej_cand_len'].isnull().astype('int')-1)
            xvcf['prob_NHEJ']=xvcf['prob_NHEJ']/(xvcf['prob_NHEJ']+xvcf['hasMH']+10**(-16))
            xvcf['prob_NHEJ']=xvcf['prob_NHEJ']/(xvcf['prob_NHEJ']+xvcf['hasMH']+10**(-16))
        xvcf['prob_NHEJ']=xvcf['prob_NHEJ'].round(16)
        
        xvcf=xvcf.reset_index()
        not_converged=True
        while not_converged:
            #Expectation step
            #Given the current estimated proportion of MMEJ, calculate the probability that each alignment is the real one
            xvcf['Estep']=xvcf['prob_NHEJ']*prob_NHEJ+(1-xvcf['prob_NHEJ'])*(1-prob_NHEJ)     
            Estep=xvcf.groupby('variant_id',sort=False).apply(lambda x: x['Estep']/math.fsum(x['Estep']))
    
            xvcf['Estep']=Estep.reset_index()['Estep']
            #Maximization step
            prob_NHEJ_t=sum(xvcf['Estep']*xvcf['prob_NHEJ'])/n_variants
            distance_from_convergence=(prob_NHEJ_t-prob_NHEJ) ** 2
            prob_NHEJ=prob_NHEJ_t
            print(f"{prob_NHEJ}")
            file1 = open(logs_outputfile, 'a')
            file1.write(str(prob_NHEJ)+"\n")
            file1.close()    
            if distance_from_convergence < convergence_threshold:
                not_converged=False
                
                proportions_df.loc[iboostrap, :] = prob_NHEJ, (1-prob_NHEJ)
                xvcf.drop(columns=['index'], inplace=True)
                xvcf.rename(columns={'Estep': 'prob_alignment'}, inplace=True)
                
                xvcf.to_csv(posteriordecodingfile,
                             sep="\t", index=False)

    proportions_df.to_csv(outfile, 
                    sep='\t', index=False)








