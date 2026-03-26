"""
A script that divide the data into batchs
and then run EMmej
"""

# importing relevant libreries
import logging
import os
import psutil
import glob
from subprocess import Popen, PIPE, STDOUT
import argparse
import datetime

import re

import pandas as pd
import numpy as np


# Construct an argument parser
all_args = argparse.ArgumentParser()
all_args.add_argument("-vpth", "--vcf_pth", required=False,
    default=0,
   help="path to input vcfs (string)")
all_args.add_argument("-v", "--vcf", required=False,
   help="path to input vcf (string)")
all_args.add_argument("-r", "--ref", required=False,
      help="path to indexed reference genome .fa file (string)")
all_args.add_argument("-o", "--outputpath", required=False,
      help="path to output file (string)")
all_args.add_argument("-bo", "--batch_outputpath", required=False,
      help="path to output file (string)")
all_args.add_argument("-m", "--mem", 
      required=False, type=int, default=1500,
      help="path to output file (string)")
all_args.add_argument("-bs", "--Blocksize", 
    required=False, default=6,
   help="number of blocks to devide per chromosome (int)")
all_args.add_argument("-anc", "--ancestral", 
    required=False, action='store_true',default=False, 
   help="indicate whether ancestral state column is available (bool)")

args = vars(all_args.parse_args())
date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
block_size = 2*10**int(args['Blocksize'])

path_to_EMmej = '/home/labs/alevy/guyta/guy_master_project/scripts/EMmej'

if args['vcf_pth'] != 0:
    files_pths = glob.glob(f"{args['vcf_pth']}/*.vcf*")
    files = []
    for f in files_pths:
        f = str(f).split(sep="/")[-1]
        if f[-3:] == '.gz' : f = f[:-3]
        if f[-4:] == '.vcf' : f = f[:-4]
        files.append(f)
    
else:
    filename = str(args['vcf']).split(sep="/")[-1]
    if filename[-3:] == '.gz' : filename = filename[:-3]
    if filename[-4:] == '.vcf' : filename = filename[:-4]
    files = [filename]


def log_subprocess_output(pipe, mode, proc):
    """
    A function the wrigts stdout and stderr from 
    subprocesses to log file
    Args:
        pipe (subprocess object): the subprocess itslf
        mode (str): err/ not err
        proc (str): subprocess name
    """
    if mode=="err": logging.info(f'\n### Traceback from {proc}: ')
    for line in iter(pipe.readline, b''): 
        logging.info(str(line)[1:-3])


print('files')
print(files)
for filename in files:
    print('========================================================')
    print(filename)
    """
    1st step: devide the dataset into batches
    """
    dir = os.path.join(f"{args['batch_outputpath']}{filename}")
    print(dir)
    if not os.path.exists(dir):
        os.mkdir(dir)
    LOG_FILENAME = f"{dir}/{date}_{block_size}bp_blocks_EMmej_batch_err_and_stdout.txt"
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG, format='%(message)s')

    whole_genome_batching_by_blocks = ['python3', f"{path_to_EMmej}/EMmej_util/EMmej_blockbootstrap/whole_genome_batching_by_blocks.py", 
            '-v', f"{args['vcf_pth']}/{filename}.vcf.gz",
            '-o', dir,
            '-bs',  f"{args['Blocksize']}"]
    if args['ancestral']: 
        whole_genome_batching_by_blocks =  whole_genome_batching_by_blocks + ['-anc']
    # print(whole_genome_batching_by_blocks)
    logging.info(f'\n### whole_genome_batching_by_blocks Stdout: ')
    process = Popen(whole_genome_batching_by_blocks, stdout=PIPE,
                            stderr=PIPE)
    process.wait()
    log_subprocess_output(process.stdout,mode='stdout', proc='whole_genome_batching_by_blocks')
    log_subprocess_output(process.stderr, mode='err', proc='whole_genome_batching_by_blocks')

    """
    2nd step: running (on wexac) EMmej over batches
    """
    if not os.path.exists(f"{dir}/EMmej_res"):
        os.mkdir(f"{dir}/{date}_{block_size}bp_blocks_EMmej_res")
    os.chdir(path_to_EMmej)

    boot_dfs = glob.glob(f"{dir}/{date}_{block_size}bp_blocks/*.vcf")

    # Submiting to Wexac
    for b_df in boot_dfs:
        filename = b_df.split(sep="/")[-1]
        bsub = ['bsub', '-q', 'new-medium', '-m', 'public_hosts',
                '-R', f"rusage[mem={args['mem']}]", '-e',
                '/home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt',
                '-o', ' /home/labs/bs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt']
        EMmej = ['python3', f"{path_to_EMmej}/EMmej.py",
                '-v', b_df, 
                '-o', f"{dir}/{date}_{block_size}bp_blocks_EMmej_res",
                '-r', f"{args['ref']}",
                '-p', 'markov', '-db', '-vr']
        if args['ancestral']: 
            EMmej =  EMmej + ['-anc']

        EMmej_submition = bsub + EMmej
        print(' '.join(EMmej_submition))
        Popen(EMmej_submition, stdout=PIPE,
                                stderr=PIPE)
        # os.chdir(f"{path_to_EMmej}/EMmej_util/EMmej_blockbootstrap")

log = f"""
This folder contains the output of: 
/home/labs/alevy/guyta/guy_master_project/scripts/EMmej/EMmej_util/EMmej_blockbootstrap/batch_data_and_submit_EMmej_to_cluster.py
On the {date}, {{datetime.datetime.now().strftime('%H:%M:%S')}}
The arguments are:
{args}
EMmej.py commands are in the following format, but for each bootstrap seperatly:
{EMmej_submition}
"""
logging.info(log)


# command example:
# python3 batch_data_and_submit_EMmej_to_cluster.py -v /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/stratification_by_annotation/GDL_indel_EMmej_format_repeated_regions.vcf.gz -bo /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/ -anc -r /home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/ref_genome/ISO1_GCF_000001215.4.fa
# python3 batch_data_and_submit_EMmej_to_cluster.py -vpth /home/labs/alevy/guyta/guy_master_project/results/arabidopsis/1001_genomes/annotation_intersection -bo /home/labs/alevy/guyta/guy_master_project/results/arabidopsis/1001_genomes/EMmej_runs/20221126_all_intersections/ -anc -r /home/labs/alevy/guyta/guy_master_project/data/arabidopsis/TAIR10_ref_genome/TAIR10_fasta/GCF_000001735.4_TAIR10.1_genomic.fna