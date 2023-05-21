import pandas as pd
from pysam import FastaFile
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import argparse
import random
import sys
#from DSBsimulate import *
import re as re
import numpy as np

import argparse

# Construct an argument parser
print("Parse arguments")
all_args = argparse.ArgumentParser()

# Add arguments to the parser
all_args.add_argument("-v", "--vcf", required=True,
   help="input vcf")
all_args.add_argument("-w", "--vcf2", required=False,default=False,
   help="input vcf 2 for simulated mixtures")
all_args.add_argument("-b", "--bootstrap", required=False,type=int,default=0,
   help="number of bootstraps")
all_args.add_argument("-x", "--nvariants_vcf1", required=False,type=int,default=100,
      help="nvariants from bootstrapped vcf1")
all_args.add_argument("-y", "--nvariants_vcf2", required=False,type=int,default=100,
      help="nvariants from bootstrapped vcf2")
all_args.add_argument("-r", "--ref", required=True,
   help="reference genome in fasta format")
all_args.add_argument("-o", "--outputfile", required=True,
      help="output file")
all_args.add_argument("-c", "--convergence", required=False,type=float,default=0.0000001,
      help="convergence threshold")
all_args.add_argument("-l", "--logs", required=False,default=False,
      help="print logs in stout")

args = vars(all_args.parse_args())

fastafile=args['ref']
vcffile=args['vcf']
vcffile2=args['vcf2']
nvariants_vcf1=args['nvariants_vcf1']
nvariants_vcf2=args['nvariants_vcf2']
outputfile=args['outputfile']
nbootstraps=args['bootstrap']
convergence_threshold=args['convergence']
printlogs=args['logs']

def count_overlapping(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count


#i=0
#chrom=xvcf.iloc[i,0]
#pos=xvcf.iloc[i,1]
#indel=xvcf.iloc[i,3]
#REF=xvcf.iloc[i,2]
#original_pos=xvcf.iloc[i,4]

def find_longest_MH(refFA,chrom,pos,indel):
    """Find the longest possible MH motif around indel
    """
    l_indel=len(indel)-1
    l_MH=0
    isMH=True
    while isMH:
        l_MH=l_MH+1
        MH_before=refFA.fetch(chrom,pos-l_MH,pos)
        MH_after=refFA.fetch(chrom,pos-1+l_indel,pos-1+l_indel+l_MH)
        if MH_before!=MH_after:
            if printlogs:
                print("1bp+MH-indel-MH+1bp:" + MH_before + "|" + indel[1:] + "|" + MH_after)
            isMH=False
            return(l_MH-1)

def MH2prob_delNHEJ(refFA,chrom,pos,indel,l_MH,win_size=2000):
    """Calculate likelihood that a given pattern occurs via NHEJ
    """
    if l_MH==0:
        return(1)
    else:
        winregion=refFA.fetch(chrom,pos-1-win_size,pos+win_size)
        MH=refFA.fetch(chrom,pos-l_MH,pos)
        prob_MH=float(count_overlapping(winregion, MH))/(win_size*2-l_MH+1)
        prob_MH=( (1-prob_MH) ** (len(indel)-1) ) * prob_MH #likelihood
        prob_MH=(1- ( (1-prob_MH) ** (len(indel)-1) )) * prob_MH #CDF of geometric
        return(float(prob_MH))

def create_variant_id(chrom,original_pos):
        """Define name (id) for a set of potential realignments of the same original indel
        """
        return(chrom+"_"+str(original_pos).split(".")[0])


#fastafile='GWHAAEV00000000.1.genome.fasta.gz'
#vcffile="/home/labs/alevy/fabrizio/workspace/guy/simulations/realigned/test13_realigned.vcf"
#vcffile2="/home/labs/alevy/fabrizio/workspace/guy/simulations/realigned/test10_realigned.vcf"

print("Import and process input files")
refFA=FastaFile(fastafile)
#xvcf=pd.read_csv(vcffile, sep='\t')
xvcf_0=pd.read_csv(vcffile, sep='\t',dtype = {'#chr': str, 'pos': int, 'REF': str, 'ALT': str,'original_pos':str})
xvcf_0['variant_id']=xvcf_0.apply(lambda row : create_variant_id(row['#chr'],row['original_pos']), axis = 1)
unique_variants_1=np.unique(np.array(xvcf_0['variant_id']))
if vcffile2:
        xvcf2_0=pd.read_csv(vcffile2, sep='\t',dtype = {'#chr': str, 'pos': int, 'REF': str, 'ALT': str,'original_pos':str})
        xvcf2_0['variant_id']=xvcf2_0.apply(lambda row : create_variant_id(row['#chr'],row['original_pos']), axis = 1)
        unique_variants_2=np.unique(np.array(xvcf2_0['variant_id']))

file1 = open(outputfile, 'w')
file1.write("prob_NHEJ\n")
file1.close()



for iboostrap in range(0,nbootstraps+1):
    if nbootstraps>0:
        sampled_vars=np.random.choice(unique_variants_1, size=nvariants_vcf1,replace=False)
        xvcf=xvcf_0[xvcf_0['variant_id'].isin(sampled_vars)]
        print("total number variants in file1 is: " + str(xvcf.shape[0]))
        if vcffile2:
            sampled_vars=np.random.choice(unique_variants_2, size=nvariants_vcf2,replace=False)
            xvcf2_b=xvcf2_0[xvcf2_0['variant_id'].isin(sampled_vars)]
            xvcf = pd.concat([xvcf, xvcf2_b], axis=0)
            print("total number variants in file1+file2 is: " + str(xvcf.shape[0]))
    else:
        xvcf=xvcf_0

    n_variants=len(np.unique(np.array(xvcf['variant_id'])))
    print("total number variants is: " + str(n_variants) + "in vcf with:" + str(xvcf.shape))

    print("Start EM")
    prob_NHEJ=0.9
    xvcf['max_length_MH']=xvcf.apply(lambda row : find_longest_MH(refFA,row['#chr'],row['pos'], row['REF']), axis = 1)
    xvcf['prob_NHEJ']=xvcf.apply(lambda row : MH2prob_delNHEJ(refFA,row['#chr'],row['pos'], row['REF'], row['max_length_MH']), axis = 1)
    xvcf=xvcf.reset_index()
    not_converged=True
    while not_converged:
        #Expectation step
        #Given the current estimated proportion of MMEJ, calculate the probability that each alignment is the real one
        xvcf['Estep']=xvcf['prob_NHEJ']*prob_NHEJ+(1-xvcf['prob_NHEJ'])*(1-prob_NHEJ)
        Estep=xvcf.groupby('variant_id',sort=False).apply(lambda x: x['Estep']/sum(x['Estep']))
        xvcf['Estep']=Estep.reset_index()['Estep']
        #Maximization step
        prob_NHEJ_t=sum(xvcf['Estep']*xvcf['prob_NHEJ'])/n_variants
        distance_from_convergence=(prob_NHEJ_t-prob_NHEJ) ** 2
        prob_NHEJ=prob_NHEJ_t
        print(prob_NHEJ)
        if distance_from_convergence < convergence_threshold:
            not_converged=False
            file1 = open(outputfile, 'a')
            file1.write(str(prob_NHEJ)+"\n")
            file1.close()
            if printlogs:
                xvcf.to_csv(outputfile+"_logs.vcf'", sep="\t",index=None)

