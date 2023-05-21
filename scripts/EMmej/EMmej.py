"""
This script is a python version of the EMmej.sh
script that operate the whole pipeline
"""

import argparse
import os
import shutil
import glob
import datetime
import logging
from subprocess import Popen, PIPE, STDOUT
import pandas as pd
# import numpy as np

# Construct an argument parser
all_args = argparse.ArgumentParser()

# Arguments that you must provide
all_args.add_argument("-v", "--vcf", required=True,
   help="path to input vcf (string)")
all_args.add_argument("-r", "--ref", required=True,
   help="path to reference genome in fasta format (string)")
all_args.add_argument("-o", "--outputpath", required=True,
      help="path to output file (string)")

# Argumrnts that are essential for EMmej functioning but has defult values
all_args.add_argument("-wp", "--patterns_windowsize", required=False, default=150,
      help="size (in bp) of the context window for pattern detection (int)")
all_args.add_argument("-wm", "--MM_windowsize", required=False, default=1000,
      help="size (in bp) of the context window for Markov model (int)")
all_args.add_argument("-p", "--probability", required=True,
      help="whether to speed up calculations by using a simplified 'Geometric' model rather than a full Markovian chain to calculate likelihoods 'Markov'")
all_args.add_argument("-s", "--steps", required=False, default="all",
      help="a flage that indicates what steps to perform {'all', 'PT' for pattern detection only, 'MM' for Markov model only, 'EM' for Expectation maximization only}, defult is all steps")

# optional arguments for RMdetector
all_args.add_argument("-mhl", "--MH_length_early_stop", required=False ,default=None, 
      help="a flag that specify the specific MH length that one wants to analyze, defult is all MH lengths")

# optional arguments, turned on by using the flag
all_args.add_argument("-vr", "--verbose", required=False, 
        action='store_true',default=False, 
      help="a flag to turn on verbose mode")
all_args.add_argument("-cr", "--CRISPR", required=False, action='store_true',default='', 
      help="a flag to indicate if the algorithm performes on CRISPR data")
all_args.add_argument("-anc", "--ancestral", required=False, action='store_true',default='', 
      help="a flag to indicate whether an ancestral state is available or not")
all_args.add_argument("-ic", "--include_context", required=False,  default='', 
      action='store_true',
      help="a flag to indicate whether to include context in output or not")
all_args.add_argument("-db", "--debugmode", required=False, action='store_true',default=False, 
      help="a flag to turn on debug mode")

# optional arguments for EM
all_args.add_argument("-b", "--bootstrap", required=False, 
        default=0, type=int,
      help="number of bootstraps to perform (int)")
all_args.add_argument("-e", "--EMalgorithm", required=False,
    default=12, type=int,
      help="type of Expectation Maximization used. choose 1 for q+markov mixed implementation (guy's version), 2 for fabrizio's version (originally developed for geometric)")
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
      help="path to posterior decoding file (string)")

args = vars(all_args.parse_args())

# handel optional args
args['CRISPR'] = '-cr' if args['CRISPR'] else None
args['ancestral'] = '-anc' if args['ancestral'] else None
args['include_context'] = '-ic' if args['include_context'] else None
args['verbose'] = '-vr' if args['verbose'] else None
flags = pd.Series([args['CRISPR'], args['ancestral'], 
        args['include_context'], args['verbose']])
flags = flags[~flags.isna()].to_list()

path_to_df = args['vcf']
path_to_output = args['outputpath']
filename = str(path_to_df.split(sep='/')[-1])[:-4]

# Creating an output folders
date = ''.join(str(datetime.datetime.today()).split()[0].split('-'))
dir = os.path.join(f"{args['outputpath']}/{filename}_EMmej_output")
if not os.path.exists(dir):
    os.mkdir(dir)
    # os.mkdir(f"{dir}/output_files")

# initiate log
# LOG_FILENAME = f"{dir}/output_files/EMmej_err_and_stdout.txt"
LOG_FILENAME = f"{dir}/EMmej_err_and_stdout.txt"
logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG, format='%(message)s')

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

def mem_reporter():
    """
    Reports the RSS and VMS of a process
    No args.
    No return, just reporting VMS and RSS at a current time point
    """
    import os, psutil
    mem = f"""Memory consumption:
    # VMS(Mb):{psutil.Process(os.getpid()).memory_info().vms / 1024 ** 2}
    # RSS(Mb):{psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2}
    """
    return mem
    # print(f'# VMS(Mb):{psutil.Process(os.getpid()).memory_info().vms / 1024 ** 2}')
    # print(f'# RSS(Mb):{psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2}')

with open(f"{dir}/EMmej_err_and_stdout.txt", 'a') as f:
        f.write(f"EMmej stderr and stdout\n")
        f.write(f"EMmej arguments:\n{args}\n")


if args['verbose']: 
    with open(f"{dir}/EMmej_err_and_stdout.txt", 'a') as f:
        f.write(f"----Verbose mode activated----\n")

if args['debugmode']:
    if not os.path.exists(f"{dir}/scripts"):
        os.mkdir(f"{dir}/scripts")
    # Copying relevant scripts:
    modules_loc = './src' ####  THIS'LL BE REMOVED FROM THE LAST VERSION ####
    scripts_loc = './scripts/EMmej' ####  THIS'LL BE REMOVED FROM THE LAST VERSION ####
    scripts = [f"{scripts_loc}/{i}" for i in ['RMdetector.py', 'Markov_model_operator.py', 
                'EM_operator.py', 'EMmej.py']]
    modules = [f"{modules_loc}/{i}" for i in ['ExpectationMaximization_q.py', 'MicroHomology_module_v3.py',
                                            'MMEJ_2nd_order_MM_v2.py', 'pysam_getfasta.py']]

    for s in scripts: shutil.copy2(s, f"{dir}/scripts")
    for m in modules: shutil.copy2(m, f"{dir}/scripts")

# Initiating general README file
README = f"""This folder contains the output of: 
/home/labs/alevy/guyta/guy_master_project/scripts/EMmej/EMmej.py
The inputs are:
input data: {path_to_df}
"""
with open(f"{dir}/EMmej_README.txt", 'w') as f:
   f.write(README)

"""
Pattern detection step
"""
if args['steps'] in ['all', 'PT']: 
    if args['steps'] == 'PT': print('Running pattern detection step only')
    RMdetector = ["python3", f'./scripts/EMmej/RMdetector.py', 
        '-v',f"{args['vcf']}", 
        '-o',f"{dir}/{filename}_RMdetector_output.tsv", 
        '-r',f"{args['ref']}", 
        '-w',f"{args['patterns_windowsize']}",
        '-mhl', f"{args['MH_length_early_stop']}"] + flags
    
    # print(print(' '.join(RMdetector)))
    
    # log start
    with open(f"{dir}/EMmej_err_and_stdout.txt", 'a') as f:
        f.write(f"\nRMdetector STARTS ({datetime.datetime.now().strftime('%H:%M:%S')})")
        f.write(f"\nArgs:\n{RMdetector}")
    
    if args['verbose']: logging.info(f'\n### Stdout from RMdetector: ')
    process = Popen(RMdetector, stdout=PIPE,
                        stderr=PIPE)
    process.wait()
    if args['verbose']: 
        with process.stdout:
            log_subprocess_output(process.stdout,mode='stdout', proc='RMdetector')
    if process.returncode != 0:
        with process.stderr:
            log_subprocess_output(process.stderr, mode='err', proc='RMdetector')
    # log end
    
    with open(f"{dir}/EMmej_err_and_stdout.txt", 'a') as f:
        f.write(f"\nRMdetector DONE ({datetime.datetime.now().strftime('%H:%M:%S')}), exit code: {process.returncode}\n")
    print(f"\nRMdetector DONE ({datetime.datetime.now().strftime('%H:%M:%S')}), exit code: {process.returncode}\n")
    
    


"""
Likelihood calculation step
"""
if args['steps'] in ['all', 'MM']:
    if args['steps'] == 'MM': print('Running Markov model step only')
    Markov_model_operator = [
        "python3", f'{scripts_loc}/Markov_model_operator.py', 
        '-v',f"{dir}/{filename}_RMdetector_output.tsv", 
        '-o',f"{dir}/{filename}_{args['probability']}_output.tsv" , 
        '-r',f"{args['ref']}", 
        '-w',f"{args['MM_windowsize']}", '-p', f"{args['probability']}"] 

    if args['include_context'] == '-ic': Markov_model_operator = Markov_model_operator + ['-ic']
    # log start
    with open(f"{dir}/EMmej_err_and_stdout.txt", 'a') as f:
        f.write(f"\nMarkov_model_operator STARTS ({datetime.datetime.now().strftime('%H:%M:%S')})")
        f.write(f"\nArgs:\n{Markov_model_operator}")
        f.write(f"\n{mem_reporter()}")
    
    if args['verbose']: logging.info(f'\n### Stdout from Markov_model_operator: ')
    process = Popen(Markov_model_operator, stdout=PIPE,
                        stderr=PIPE)
    process.wait()
    if args['verbose']: 
        with process.stdout:
            log_subprocess_output(process.stdout,mode='stdout', proc='Markov_model_operator')
    if process.returncode != 0:
        with process.stderr:
            log_subprocess_output(process.stderr, mode='err', proc='Markov_model_operator')

    # log end
    with open(f"{dir}/EMmej_err_and_stdout.txt", 'a') as f:
        f.write(f"\nMarkov_model_operator DONE ({datetime.datetime.now().strftime('%H:%M:%S')}), exit code: {process.returncode}\n")
        f.write(f"\n{mem_reporter()}")
    print(f"\nMarkov_model_operator DONE ({datetime.datetime.now().strftime('%H:%M:%S')}), exit code: {process.returncode}\n")


"""
Parameter estimation step
"""
if args['steps'] in ['all', 'EM']:
    if args['steps'] == 'EM': print('Running Expectation maximization step only')
    EM_flags = [
        "-e", f"{args['EMalgorithm']}", "-b" , f"{args['bootstrap']}", 
        "-bs", f"{args['bootsize']}", "-bs2",f"{args['bootsize2']}", 
        "-i", f"{args['initial_proportion_NHEJ']}", "-c", f"{args['convergence']}", 
        "-l", f"{dir}/{filename}_{args['probability']}_EM_output_log.tsv",
        "-pd", f"{dir}/{filename}_{args['probability']}_EM_output_posterior_decoding.tsv"]
    
    if args['verbose']: EM_flags = EM_flags + ['-vr']
    
    EM_operator = [
        "python3",f'{scripts_loc}/EM_operator.py', 
        '-v',f"{dir}/{filename}_{args['probability']}_output.tsv", 
        '-o',f"{dir}/{filename}_{args['probability']}.tsv", 
        '-r',f"{args['ref']}", 
        '-w', f"{args['patterns_windowsize']}" ] + EM_flags
    # log start
    with open(f"{dir}/EMmej_err_and_stdout.txt", 'a') as f:
        f.write(f"\nEM_operator STARTS ({datetime.datetime.now().strftime('%H:%M:%S')})")
        f.write(f"\nArgs:\n{EM_operator}")
    if args['verbose']: logging.info(f'\n### Stdout from EM_operator: ')
    process = Popen(EM_operator, stdout=PIPE,
                        stderr=PIPE)
    process.wait()
    if args['verbose']: 
        with process.stdout:
            log_subprocess_output(process.stdout,mode='stdout', proc='EM_operator')
    if process.returncode != 0:
        with process.stderr:
            log_subprocess_output(process.stderr, mode='err', proc='EM_operator')

    # log end
    with open(f"{dir}/EMmej_err_and_stdout.txt", 'a') as f:
        f.write(f"\nEM_operator DONE ({datetime.datetime.now().strftime('%H:%M:%S')}), exit code: {process.returncode}\n")
        f.write(f"\n{mem_reporter()}")
    print(f"\nEM_operator DONE ({datetime.datetime.now().strftime('%H:%M:%S')}), exit code: {process.returncode}\n")


# execution example: python3 EMmej.py -v /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/test_data/EMmej_test_data_indels.vcf -o /home/labs/alevy/guyta/guy_master_project/scripts/EMmej_merging_fabrizio_and_guy/EMmej_runs_outputs/EMmej_py_output -r /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/test_data_reference_genome/EMmej_test_data_ref_genome.fa -p markov -db -cr
# execution example: python3 EMmej.py -v /home/labs/alevy/fabrizio/workspace/guy/simulations/test12_subsampled.vcf -o /home/labs/alevy/guyta/guy_master_project/scripts/EMmej_merging_fabrizio_and_guy/EMmej_runs_outputs/EMmej_py_output -r /home/labs/alevy/fabrizio/workspace/guy/simulations/GWHAAEV00000000.1.genome.fasta -p markov -db -cr