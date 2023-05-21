"""
This scrip performs blockbootstraping using
data that is already devided into blocks 
(or batchs) and analyzed with EMmej until
at list the likelihood calculation done,
EM results is'nt necessery for this script
"""

# importing relevant libreries
import argparse
from asyncio import subprocess
import datetime
import os
import os.path
import glob
from subprocess import Popen, PIPE, STDOUT

import re

import pandas as pd
import numpy as np

# Construct an argument parser
all_args = argparse.ArgumentParser()
all_args.add_argument("-v", "--vcf", required=True,
   help="path to input vcf (string)")
all_args.add_argument("-r", "--ref", required=False,
      help="path to indexed reference genome .fa file (string)")
all_args.add_argument("-o", "--outputpath", required=False,
      help="path to output file (string)")
all_args.add_argument("-b", "--bootnumber", 
      required=False, type=int, default=100,
      help="path to output file (string)")
all_args.add_argument("-m", "--mem", 
      required=False, type=int, default=1500,
      help="memory required")
all_args.add_argument("-d", "--date", 
      required=False, type=int, default=10,
      help="date")
all_args.add_argument("-sil", "--stratify_by_indel_length", 
            required=False, action='store_true',default=False, 
      help="a flag to turn on stratification by indel length")
all_args.add_argument("-filter_1bp", "--filter_1bp_indels", 
            required=False, action='store_true',default=False, 
      help="a flag to turn on stratification by indel length")


args = vars(all_args.parse_args())

# def intersection(lst1, lst2):
#     lst3 = [value for value in lst1 if value in lst2]
#     return lst3
def parallel(lst1, lst2):
    lst3 = [value for value in lst1 if value not in lst2]
    return lst3
if args['date'] == 10:
      date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
else: date = args['date']

# Creating output folders
dir = os.path.join(f"{args['outputpath']}/{date}_blockbootstrap_output")
if not os.path.exists(dir):
    os.mkdir(dir)
    os.mkdir(f"{dir}/bootstrap_dfs")
    os.mkdir(f"{dir}/EM_res")
    os.mkdir(f"{dir}/insertions_proportions")

# Creating a list of all folders that contains blockbootstraping runs
dir_ls = glob.glob(f"{args['vcf']}/*")
cols = ['CHR', 
        'POS', 'original_pos','variant_id', 'variant_id_N','direction','ancestral_indel', 'derived_indel', 
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
EMmej_MM_outputs = pd.DataFrame(columns=cols)

# Empty list to get names of runs that havn't finished yet
undone_runs = []

for i, fol in enumerate(dir_ls):
      output_ls = glob.glob(f"{fol}/*")    
      """
      Markov model merging
      """
      MM_outputs = [f for f in os.listdir(f"{fol}/") if re.match(r'(.+?)_markov_output.tsv', f)]
      if len(MM_outputs) > 0:
            MM_output_df = pd.read_csv(f"{fol}/{MM_outputs[0]}", 
                  sep='\t')
            chr_name = re.search(r'chr_(.+?)_blockN', fol).group(1)
            block_n = re.search(r'blockN(.+?)_', fol).group(1)
            MM_output_df['blockID'] = f"{chr_name}_{block_n}"
            EMmej_MM_outputs = pd.concat([EMmej_MM_outputs, MM_output_df])
      else: undone_runs.append(fol)

print(len(undone_runs))
print(undone_runs)
      
EMmej_MM_outputs.reset_index(drop=True, inplace=True)
if args['filter_1bp_indels']:
      EMmej_MM_outputs = EMmej_MM_outputs.loc[
            (~EMmej_MM_outputs['indel_len'].isin([1,-1])), :]
"""
If one wants to stratify by indel length bins,
then uncomment those lines
"""
if args['stratify_by_indel_length']: 
      conditions = [
      EMmej_MM_outputs['indel_len'].isin(range(-10,11,1)),
      EMmej_MM_outputs['indel_len'] > 10,
      EMmej_MM_outputs['indel_len'] < (-10)
      ]

      choices = [
      EMmej_MM_outputs['indel_len'],
      '[11:inf]',
      '[inf:-11]',
      ]
      EMmej_MM_outputs.loc[:,'indel_len_bin'] = np.select(conditions, choices)
      EMmej_MM_outputs.loc[:,'indel_len_bin'] = EMmej_MM_outputs.loc[
            :,'indel_len_bin'].astype('category')


"""
Creating the bootstraps (i.e. the fake genomes)
"""
file_list = []
for i in range(args['bootnumber']):
    fname = f"{dir}/bootstrap_dfs/bootnumber_{i}.tsv"
    file_list.append(fname)
    if not os.path.isfile(fname):
        print(fname)
        samp_blocks = np.random.choice(a=EMmej_MM_outputs['blockID'].unique(),
                                size=EMmej_MM_outputs['blockID'].unique().shape[0],
                                replace=True)
        blockbootstrap_df = pd.DataFrame(
                columns=EMmej_MM_outputs.columns.to_list().remove('blockID'))
        for idx,b in enumerate(samp_blocks):
            samp_df = EMmej_MM_outputs.loc[
                EMmej_MM_outputs['blockID'] == b, :].copy()
            bn = str(str(b).split('_')[-1])
            idx_s = str(idx)
            # samp_df['variant_id'] = samp_df.apply(
            #       lambda row: f"{row['variant_id']}_{str(b).split('_')[-1]}_{idx}", axis=1)
            samp_df['variant_id'] = samp_df['variant_id'] + '_' + bn + "_" + idx_s
            samp_df.drop(columns=['blockID'], inplace=True)
            blockbootstrap_df = pd.concat([blockbootstrap_df, samp_df])
        
        """
        If one wants to stratify by indel length bins,
        then uncomment those lines
        """
        if args['stratify_by_indel_length']:
            for indel_l_bin in blockbootstrap_df['indel_len_bin'].unique(): 
                if '[' in str(indel_l_bin): indel_l_bin_lable = indel_l_bin[1:-1]
                else: indel_l_bin_lable = indel_l_bin
                fname_bin = f"{dir}/bootstrap_dfs/bootnumber_{i}_bin_{indel_l_bin_lable}.tsv"
                file_list.append(fname_bin)
                if not os.path.isfile(fname_bin):
                    blockbootstrap_df_bin = blockbootstrap_df.loc[
                            (blockbootstrap_df['indel_len_bin'] == indel_l_bin), :].copy()
                    blockbootstrap_df_bin.drop(columns=['indel_len_bin'], 
                                            inplace=True)
                    blockbootstrap_df_bin.to_csv(
                            f"{dir}/bootstrap_dfs/bootnumber_{i}_bin_{indel_l_bin_lable}.tsv",
                            sep='\t', index=False)
        
        blockbootstrap_df.drop(columns=['indel_len_bin'], 
                                            inplace=True)
        blockbootstrap_df.to_csv(fname, 
                            sep='\t', index=False )
            

"""
Running the EM over the bootstrap dfs
"""
boot_dfs = glob.glob(f"{dir}/bootstrap_dfs/*")
# EMmej_submition_to_cluster = '/EMmej_submition_to_cluster.sh'
script_path ='/home/labs/alevy/guyta/guy_master_project/scripts/EMmej'

unfinished_files = []
unfinished_EM = []
for b_df in file_list:
    filename = b_df.split(sep="/")[-1]
    if not os.path.isfile(b_df):
        unfinished_files.append(filename)
        continue
    else:
        if not os.path.isfile(f"{dir}/insertions_proportions/{filename[:-4]}_insertions_proportions.tsv"):
            print(b_df)
            # Estimating Insertions proportions
            # Setting cutoffs based on simulations results
            
            if args['filter_1bp_indels']:
                snap_ins_cutoff = 0.01
                loop_ins_cutoff = 0.01
            else:
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
            print('====================')
            print(f"{dir}/insertions_proportions/{filename[:-4]}_insertions_proportions.tsv")


        if not os.path.isfile(f"{dir}/EM_res/{filename[:-4]}.tsv"):
            # unfinished_EM.append(filename)
            bsub = ['bsub', '-q', 'new-medium', '-m', 'public_hosts', '-R', f"rusage[mem={args['mem']}]",
            '-e', '/home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt' ,
            '-o' ,'/home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt']
            Em_command = ['python3', f"{script_path}/EM_operator.py", '-v', b_df, 
            '-o', f"{dir}/EM_res/{filename[:-4]}.tsv", '-r',  args['ref'], '-e', '12']
            EM_submition_to_cluster_command = bsub + Em_command
            Popen(EM_submition_to_cluster_command, 
                    stdout=PIPE, stderr=PIPE, text=True)



"""
Adding a README
"""
# if len(unfinished_files) != 0:
README = f"""This folder contains the output of: 
/home/labs/alevy/guyta/guy_master_project/scripts/EMmej/EMmej_util/EMmej_blockbootstrap.py
On the {date}, {{datetime.datetime.now().strftime('%H:%M:%S')}}
The arguments are:
{args}

"""
# EM commands are in the following format, but for each bootstrap seperatly:
# {EM_submition_to_cluster_command}

with open(f"{dir}/EMmej_blockbootstrap_README.txt", 'w') as f:
    f.write(README)        


if len(unfinished_files) != 0:
    undone = f"""
    List of unfinished files in /bootstrap_dfs:
    {unfinished_files}
    -------------------------------------------------------------------
    """

    with open(f"{dir}/EMmej_blockbootstrap_undone_list.txt", 'a') as f:
        f.write(undone)  

# if len(unfinished_EM) != 0:
#     undone = f"""
#     List of unfinished files in /EM_res:
#     {unfinished_EM}
#     """

#     with open(f"{dir}/EMmej_blockbootstrap_undone_EM.txt", 'w') as f:
#         f.write(undone)  

# Command example:
# python3 EMmej_blockbootstrap.py -v /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/GDL_indel_EMmej_format_non_repeated_regions_exon_no_MT/20221021 -o /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/blockbootstraping/GDL_indel_EMmej_format_non_repeated_regions_exon_no_MT/ -r /home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/ref_genome/ISO1_GCF_000001215.4.fa -b 100 -sil
# bsub -q new-short -e /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt -o /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt -m public_hosts -R "rusage[mem=5000]" "python3 EMmej_blockbootstrap.py -v /home/labs/alevy/guyta/guy_master_project/results/arabidopsis/1001_genomes/EMmej_runs/1001genomes_minor_allele_count_2_no1bp_del/1001genomes_non_repeats_exons_EMmej_formated_data/20221121_2000000bp_blocks_EMmej_res -o /home/labs/alevy/guyta/guy_master_project/results/arabidopsis/1001_genomes/EMmej_runs/1001genomes_minor_allele_count_2_no1bp_del/1001genomes_non_repeats_exons_EMmej_formated_data/ -r /home/labs/alevy/guyta/guy_master_project/data/arabidopsis/TAIR10_ref_genome/TAIR10_fasta/GCF_000001735.4_TAIR10.1_genomic.fna -b 100 -sil"
# bsub -q new-short -e /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt -o /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt -m public_hosts -R "rusage[mem=5000]" "python3 EMmej_blockbootstrap.py -v /home/labs/alevy/guyta/guy_master_project/results/arabidopsis/1001_genomes/EMmej_runs/20221126_all_intersections/1001genomes_non_repeats_exons_EMmej_formated_data/20221126_2000000bp_blocks_EMmej_res -o /home/labs/alevy/guyta/guy_master_project/results/arabidopsis/1001_genomes/EMmej_runs/20221126_all_intersections/1001genomes_non_repeats_exons_EMmej_formated_data/ -r /home/labs/alevy/guyta/guy_master_project/data/arabidopsis/TAIR10_ref_genome/TAIR10_fasta/GCF_000001735.4_TAIR10.1_genomic.fna -b 5 -sil"
# bsub -q new-short -e /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt -o /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt -m public_hosts -R "rusage[mem=5000]" "python3 EMmej_blockbootstrap.py -v /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/20230103_EMmej_with_50bp_markov_window_size/Clark_GDL_EMmej_formated_non_repeated_regions_non_genic/20221126_2000000bp_blocks_EMmej_res/ -r /home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/ref_genome/ISO1_GCF_000001215.4.fa -b 100 -sil"
# bsub -q new-short -e /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt -o /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt -m public_hosts -R "rusage[mem=5000]" "python3 EMmej_blockbootstrap.py -v /home/labs/alevy/guyta/guy_master_project/results/tomato/RTGR/EMmej_runs/RTGR_EMmej_formated_non_repeated_regions_exons/20230106_2000000bp_blocks_EMmej_res -r /home/labs/alevy/guyta/guy_master_project/data/tomato/SL3_ref_genome/GCF_000188115.4_SL3.0_genomic.fna -o /home/labs/alevy/guyta/guy_master_project/results/tomato/RTGR/EMmej_runs/RTGR_EMmej_formated_non_repeated_regions_exons -b 100 -sil"