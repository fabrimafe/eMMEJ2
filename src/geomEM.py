#import pandas as pd
#from pysam import FastaFile
#from Bio import pairwise2
#from Bio.pairwise2 import format_alignment
#import random
#import re as re
#import numpy as np


def count_overlapping(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count


def find_longest_MH(refFA,chrom,pos,indel):
    """Find the longest possible MH motif around indel
    """
    printlogs=False
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


