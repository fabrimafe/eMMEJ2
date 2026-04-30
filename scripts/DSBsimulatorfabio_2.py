import random
from pysam import FastaFile
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import pandas as pd
import argparse

# Construct an argument parser
print("Parse arguments")
all_args = argparse.ArgumentParser()

# Add arguments to the parser
all_args.add_argument("-r", "--ref", required=True,
   help="reference genome in fasta format")
all_args.add_argument("-o", "--out", required=False,
   help="output file")
all_args.add_argument("-c", "--chrom", required=False,
   help="chrom of position of DSB; necessary of nsims=1, otherwise ignored")
all_args.add_argument("-p", "--pos", required=False,type=int,
   help="position of DSB; necessary if nsims=1; otherwise ignored")
all_args.add_argument("-n", "--nsims", required=False,default=100,
      help="number of random DSB to simulate")
all_args.add_argument("-m", "--mechanism", required=True,default="deletion_MMEJ",
   help="mechanism of DSB-repair. One of: deletion_NHEJ,deletion_MMEJ,SDloopout,SDsnapback,insertion,MMEJtrans")
all_args.add_argument("-l", "--MHlength", required=False,type=int,default=2,
   help="length of microhomology motif ( or minimum indel size ) ")
all_args.add_argument("-s", "--SDlength", required=False,type=int,default=5,
      help="length of SD motif/distance second MH in MMEJtrans")
all_args.add_argument("-d", "--maxdistance", required=False,type=int,default=10,
      help="maximum distance of MH motif or maximum indel length")
all_args.add_argument("-z", "--Z_length", required=False,type=int,default=10,
        help="maximum distance between MH and SD")

args = vars(all_args.parse_args())

def seq2complement(xseq):
    """Calculate complement of a sequence and NOT reverse complement"""
    x_complement = ""
    for i in xseq:
        if i == "T" or i == "t":
            x = "A"
        elif i == "A" or i == "a":
            x = "T"
        elif i == "G" or i == "g":
            x = "C"
        elif i == "C" or i == "c":
            x = "G"
        else:
            x = "N"
        x_complement = x_complement + x
    return(x_complement)

def seq2revcomplement(xseq):
    """Calculate reverse complement of a DNA sequence"""
    return(seq2complement(xseq)[::-1])

def normalize_vcf(anc, der):
    
    while len(anc) > 1 and len(der) > 1 and anc[-1] == der[-1]:
        anc = anc[:-1]
        der = der[:-1]
    
    
    while len(anc) > 1 and len(der) > 1 and anc[0] == der[0]:
        anc = anc[1:]
        der = der[1:]
    
    return anc, der

def fa2deletion(refFA,chrom,pos,min_length,max_length):
    """Function to create deletions on a fasta file giving outuput in vcf format. First and last 2kbp of chromosomes should not be used.
    Arguments are:
    refFA: a reference fasta file in pysam format
    chrom: the chromosome to be mutated
    pos: the start (0-based) position of the mutation, i.e. if you want to delete position 121 and 122 you need to set 120
    min_length: minimum length of indel
    max_length: minimum length of indel
    Output has indel_length in addition to vcf-like annotation.
    """
    indel_length=random.sample(range(min_length,max_length+1),1)[0]
    indel_seq=refFA.fetch(chrom,pos,pos+indel_length)
#    seq_indel120bp=refFA.fetch(chr,pos-120,pos+120+indel_length)
#    seq_ref2kbp=refFA.fetch(chr,pos-2000,pos+2000)
#    seq_ref120bp=len(refFA.fetch(chr,pos-120,pos+120))
    ancestral_state_vcf=refFA.fetch(chrom,pos-1,pos+indel_length)
    derived_state_vcf=refFA.fetch(chrom,pos-1,pos)
    indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
    if (indelNs>-3):
        print("N found in indel")
        return('error')
    else:
        return(chrom,pos,ancestral_state_vcf,derived_state_vcf,indel_length)

def fa2deletion_MMEJ(refFA,chrom,pos,MH_length,ma_distance_MMEJ):
    """Function to create deletions arising from MMEJ on a fasta file giving outuput in vcf format. First and last 2kbp of chromosomes should not be used.
    Arguments are:
    refFA: a reference fasta file in pysam format
    chr: the chromosome to be mutated
    pos: the start (0-based) position of the DSB, i.e. if you want to delete position 121 and 122 you need to set 120
    MH_length: length of Microhomology motif
    max_distance_MMEJ: max distance between cut site and Microhomology motif
    """
    MH_seq=refFA.fetch(chrom,pos-MH_length,pos)
    seq_window=refFA.fetch(chrom,pos,pos+MH_length+max_distance_MMEJ)
    i_MH=seq_window.find(MH_seq)
    if i_MH < 0:
        print("no MH possible in range")
        return("error")
    else:
        indel_length=i_MH+MH_length
        indel_seq=refFA.fetch(chrom,pos,pos+i_MH+MH_length)
        #seq_indel120bp=refFA.fetch(chr,pos-120,pos-1)+refFA.fetch(chr,pos+indel_length-1,pos+indel_length-1+120) #pos+120+indel_length)
        #seq_ref2kbp=refFA.fetch(chr,pos-2000,pos+2000)
        #seq_ref120bp=len(refFA.fetch(chr,pos-120,pos+120))
        print("MH motif: "+MH_seq)
        print("MH|DEL|MH:")
        print(refFA.fetch(chrom,pos-MH_length,pos)+"|"+refFA.fetch(chrom,pos,pos+indel_length-MH_length)+"|"+refFA.fetch(chrom,pos+indel_length-MH_length,pos+indel_length) )
        #return(chr,pos,"DEL",indel_seq,"-",indel_length,seq_indel120bp,121,seq_ref2kbp,seq_ref120bp)
        ancestral_state_vcf=refFA.fetch(chrom,pos-1,pos+indel_length)
        derived_state_vcf=refFA.fetch(chrom,pos-1,pos)
        indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
        if (indelNs>-3):
            print("N found in indel")
            return('error')
        else:
            return(chrom,pos,ancestral_state_vcf,derived_state_vcf,indel_length)

def fa2insertion(refFA,chrom,pos,min_length,max_length):
    """Function to create insertions on a fasta file giving outuput in vcf format. First and last 2kbp of chromosomes should not be used.
    Insertion sequences are generated by transposition, copying a region in the nearby 2000bp.
    Arguments are:
    refFA: a reference fasta file in pysam format
    chr: the chromosome to be mutated
    pos: the start (0-based) position of the mutation, i.e. if you want to delete position 121 and 122 you need to set 120
    min_length: minimum length of indel
    max_length: minimum length of indel
    """
    indel_length=random.sample(range(min_length,max_length+1),1)[0]
    insertion_from=random.sample(range(500,1900),1)[0]
    indel_seq=refFA.fetch(chrom,insertion_from,insertion_from+indel_length)
#    seq_indel120bp=refFA.fetch(chr,pos-120,pos) + indel_seq + refFA.fetch(chr,pos+1,pos+120)
#    seq_ref2kbp=refFA.fetch(chr,pos-2000,pos+2000)
#    seq_ref120bp=len(refFA.fetch(chr,pos-120,pos+120))
#    print("original context (+/-1bp): "+refFA.fetch(chrom,pos-1,pos+1)) #print around
#    print("mutated context: "+refFA.fetch(chrom,pos-1,pos) + indel_seq + refFA.fetch(chrom,pos,pos+1)) #print insertion and around
    ancestral_state_vcf=refFA.fetch(chrom,pos-1,pos)
    derived_state_vcf=refFA.fetch(chrom,pos-1,pos)+indel_seq
#    return(chr,pos,indel_seq,indel_length,seq_indel120bp,121,seq_ref2kbp,seq_ref120bp)
    indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
    print(indelNs)
    if (indelNs>-3 or indel_length==0 or ancestral_state_vcf=="N" or ancestral_state_vcf=="-" or ancestral_state_vcf=="n"):
        print("N or NULL found in indel")
        return('error')
    else:
        return(chrom,pos,ancestral_state_vcf,derived_state_vcf,indel_length)


def fa2insertion_snapback(refFA,chrom,pos,MH_length,SD_motif_length,max_distance,pre_SD_length=0):
    """Function to create snapback insertions on a fasta file giving outuput in vcf-like entries. First and last 2kbp of chromosomes should not be used.
    Arguments are:
    refFA: a reference fasta file in pysam format
    chr: the chromosome to be mutated
    pos: the start (0-based) position of the DSB
    MH_length: length of Microhomology motif
    SD_motif_length: length of SD motif
    max_distance: starts from end deletion of step1
    pre_SD_length: standard snapback have SD starting right after Microhomology motif. By setting this to >0 one can change that.
    Output has indel_length in addition to vcf-like annotation.
    """
    SD_motif=refFA.fetch(chrom,pos+pre_SD_length,pos+SD_motif_length+pre_SD_length)
    MH_seq=refFA.fetch(chrom,pos-MH_length,pos)
    seq_window=refFA.fetch(chrom,pos+SD_motif_length+pre_SD_length,pos+pre_SD_length+max_distance-1)
    seq_window_fabio=refFA.fetch(chrom,pos-MH_length,pos+max_distance-1)
    #seq_window_fabio is only for the print
    #1° MH is always in position 0 in seq_window_fabio
    SD_revc=seq2revcomplement(SD_motif)
    i_SD=seq_window.find(SD_revc)+pre_SD_length
    print("seq window from 1° MH : "+seq_window_fabio)
    print("MH seq: "+MH_seq)
    print("SD motif: "+SD_motif)
    print("position 2nd SD motif from 1° SD motif, (snapback_loop): "+str(i_SD))
    #print("MH+SD+..+SDrevc: "+refFA.fetch(chr,pos-MH_length,pos+i_SD+SD_motif_length))
    if i_SD < 0:
        print("no SD possible in range")
        return("error")
    else:
        i_SD=i_SD+SD_motif_length
        #here i_SD is distance between the last nt 1° SD and  2° SD 
        seq_window=refFA.fetch(chrom,pos+i_SD+SD_motif_length,pos+max_distance)
        MH_revc=seq2revcomplement(MH_seq)
        i_MH=seq_window.find(MH_revc)
        if i_MH < 0:
            print("no MH possible in range")
            return("error")
        else:
            i_MH=i_MH+i_SD+SD_motif_length
            #i_MH distance between 1° SD and 2° MH
            indel_length=i_MH-i_SD-SD_motif_length
            #indel_length defined as the distance between 2° SD and 2° MH 
            print("distance between SDrevc and MHrevc (i_MH) == length of insertion : "+str(indel_length))
            print("MH|SD|snapback-loop|SDrevc|INSrevc|MHrevc: ")#refFA.fetch(chr,pos-MH_length,pos+i_MH+MH_length))
            indel_seq_annotated=refFA.fetch(chrom,pos-MH_length,pos)+"|"+refFA.fetch(chrom,pos,pos+SD_motif_length) +"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_SD)+"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length)+"|"+refFA.fetch(chrom,pos+i_SD+SD_motif_length,pos+i_MH) +"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)
            indel_seq=seq2revcomplement(refFA.fetch(chrom,pos+i_SD+SD_motif_length,pos+i_MH))
            indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
            seq_after_insertion=refFA.fetch(chrom,pos-MH_length,pos)+"|"+str(indel_seq)+"|"+refFA.fetch(chrom,pos,pos+SD_motif_length) +"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_SD)+"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length)+"|"+refFA.fetch(chrom,pos+i_SD+SD_motif_length,pos+i_MH) +"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)
            ancestral_state=refFA.fetch(chrom,pos-1,pos)
            derived_state=refFA.fetch(chrom,pos-1,pos)+indel_seq
            ancestral_state_vcf=ancestral_state.upper()
            derived_state_vcf=derived_state.upper()
            if (indelNs>-3 or indel_length==0 or ancestral_state_vcf=="N" or ancestral_state_vcf=="-" or ancestral_state_vcf=="n"):
                print("N found in indel")
                return('error')
            else:
                print(indel_seq_annotated)
                print("after insertion:")
                print("MH|INS|SD|snapback-loop|SDrevc|INSrevc|MHrevc")
                print(seq_after_insertion)
                return(chrom,pos,ancestral_state_vcf,derived_state_vcf,indel_length)



def fa2deletion_snapback(refFA,chrom,pos,MH_length,SD_motif_length,max_distance,pre_SD_length=0):
    """Function to create snapback deletions on a fasta file giving outuput in vcf-like entries.
    It happens when Z>0 and X=0.
    First and last 2kbp of chromosomes should not be used.
    Arguments are:
    refFA: a reference fasta file in pysam format
    chr: the chromosome to be mutated
    pos: the start (0-based) position of the DSB
    MH_length: length of Microhomology motif
    SD_motif_length: length of SD motif
    max_distance: starts from end deletion of step1
    pre_SD_length: standard snapback have SD starting right after Microhomology motif. By setting this to >0 one can change that.
    here pos indicate the position of rhe MH|Sd that mediates the deletion
    Output has indel_length in addition to vcf-like annotation.
    """
    SD_motif=refFA.fetch(chrom,pos+pre_SD_length,pos+SD_motif_length+pre_SD_length)
    MH_seq=refFA.fetch(chrom,pos-MH_length,pos)
    seq_window=refFA.fetch(chrom,pos+SD_motif_length+pre_SD_length,pos+pre_SD_length+max_distance-1)
    seq_window_fabio=refFA.fetch(chrom,pos-MH_length,pos+max_distance-1)
    #seq_window_fabio is only for the print
    #1° MH is always in position 0 in seq_window_fabio
    SD_revc=seq2revcomplement(SD_motif)
    i_SD=seq_window.find(SD_revc)+pre_SD_length
    print("seq window: "+seq_window_fabio)
    print("MH seq: "+MH_seq)
    print("SD motif: "+SD_motif)
    print("distance between SD motif and SDrevc motif, (snapback-loop): "+str(i_SD))
    #print("MH+SD+..+SDrevc: "+refFA.fetch(chr,pos-MH_length,pos+i_SD+SD_motif_length))
    if i_SD < 0:
        print("no SD possible in range")
        return("error")
    else:
        i_SD=i_SD+SD_motif_length
        seq_window=refFA.fetch(chrom,pos+i_SD+SD_motif_length,pos+max_distance)
        MH_revc=seq2revcomplement(MH_seq)
        i_MH=seq_window.find(MH_revc)
        if i_MH < 0:
            print("no MH possible in range")
            return("error")
        else:
            i_MH=i_MH+i_SD+SD_motif_length
            indel_length=i_MH-i_SD-SD_motif_length
            print("distance between SDrevc and MHrevc (i_MH) == length deletion: "+str(indel_length))
            print("MH|SD|snapback-loop|SDrevc|DEL|MHrevc: ")#refFA.fetch(chr,pos-MH_length,pos+i_MH+MH_length))
            indel_seq_annotated=refFA.fetch(chrom,pos-MH_length,pos)+"|"+refFA.fetch(chrom,pos,pos+SD_motif_length) +"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_SD)+"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length)+"|"+refFA.fetch(chrom,pos+i_SD+SD_motif_length,pos+i_MH) +"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)
            indel_seq=seq2revcomplement(refFA.fetch(chrom,pos+i_SD+SD_motif_length,pos+i_MH))
            indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
            seq_after_deletion=refFA.fetch(chrom,pos-MH_length,pos)+"|"+refFA.fetch(chrom,pos,pos+SD_motif_length) +"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_SD)+"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length)+"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)
            derived_state=refFA.fetch(chrom,pos-1,pos)
            ancestral_state=refFA.fetch(chrom,pos-1,pos)+indel_seq
            ancestral_state_vcf=ancestral_state.upper()
            derived_state_vcf=derived_state.upper()
            if (indelNs>-3 or indel_length==0 or ancestral_state_vcf=="N" or ancestral_state_vcf=="-" or ancestral_state_vcf=="n"):
                print("N found in indel")
                return('error')
            else:
                print(indel_seq_annotated)
                print("after deletion")
                print("MH|SD|snapback-loop|SDrevc|MHrevc: ")
                print(seq_after_deletion)
                return(chrom,pos,ancestral_state_vcf,derived_state_vcf,indel_length)



def fa2SD_direct_insertion(refFA,chrom,pos,MH_length,SD_motif_length,max_distance,pre_SD_length=0):
    """Function to create scars of a SD-direct insertion, no matter if trans or loop-out.
    It happens when Z=0, X>0.
    DSB repair on a fasta file giving outuput in vcf-like entries. First and last 2kbp of chromosomes should not be used.
    Arguments are:
    refFA: a reference fasta file in pysam format
    chr: the chromosome to be mutated
    pos: the start (0-based) position of the DSB
    MH_length: length of Microhomology motif
    SD_motif_length: length of SD motif
    max_distance: starts from end deletion of step1
    pre_SD_length: standard SD-loopouts have SD starting right after Microhomology motif. By setting this to >0 one can change that.
    Output has indel_length in addition to vcf-like annotation.
    """
    MH_seq=refFA.fetch(chrom,pos-MH_length,pos)
    SD_motif=refFA.fetch(chrom,pos+pre_SD_length,pos+pre_SD_length+SD_motif_length)
    seq_window=refFA.fetch(chrom,pos+SD_motif_length+pre_SD_length,pos+max_distance-1)
    seq_window_fabio=refFA.fetch(chrom,pos-MH_length,pos+max_distance - 1)
    i_MH=seq_window.find(MH_seq) #=position of MH from SD in seq window
    print("seq window from 1° MH: "+seq_window_fabio)
    print("MH seq: "+MH_seq)
    print("SD motif: "+SD_motif)
    print("distance between 1°SD and 2° MH,spacer seq ,(i_MH): "+str(i_MH))
#    print("MH++SD+..+MH: "+refFA.fetch(chr,pos-MH_length,pos+i_MH+MH_length))
#the 1° MH is always at 0 pos in seq_window_fabio
    if i_MH < 0:
        print("no MH possible in range")
        return ("error")
    else:
        i_MH=i_MH+SD_motif_length+pre_SD_length
        #here i_MH is the pos of the last nt of 2° SD from 2° MH
        seq_window=refFA.fetch(chrom,pos+i_MH+MH_length,pos+max_distance)
        #up here seq_window is from 2°MH to max distance, it's used to find 2°SD
        i_SD=seq_window.find(SD_motif)
        if i_SD < 0:
            print("no SD possible in range")
            return "error" 
        else:
            i_SD=i_SD+i_MH+MH_length
             #up here i_sd is the position of 2° SD from the last nt of 1° MH
            indel_length=i_SD-i_MH-MH_length
            #indel= distance from 2° SD and the 1° SD - distance from 2° SD and 2°MH - length MH
            print("distance between 2° MH and 2° SD(i_SD) == length insertion: "+str(indel_length))
            print("MH|SD|spacer seq between two patterns|MH|INS|SD: ") #+refFA.fetch(chr,pos-MH_length,pos+i_MH+MH_length+i_SD+SD_motif_length))
            indel_seq_annotated=refFA.fetch(chrom,pos-MH_length,pos,)+"|"+refFA.fetch(chrom,pos,pos+SD_motif_length) +"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_MH)+"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)+"|"+refFA.fetch(chrom,pos+i_MH+MH_length,pos+i_SD) +"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length)
            
            seq_after_insertion=refFA.fetch(chrom,pos-MH_length,pos,)+"|"+refFA.fetch(chrom,pos+i_MH+MH_length,pos+i_SD) +"|"+refFA.fetch(chrom,pos,pos+SD_motif_length)+"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_MH)+"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)+"|"+refFA.fetch(chrom,pos+i_MH+MH_length,pos+i_SD) +"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length)

            indel_seq=refFA.fetch(chrom,pos+i_MH+MH_length,pos+i_SD)
            #indel_seq = from pos, + distance between 2° SD and 2° MH + length of MH to 2°SD from the the last nt of 2° MH
            ancestral_state=refFA.fetch(chrom,pos-1,pos)
            derived_state=refFA.fetch(chrom,pos-1,pos)+indel_seq
            ancestral_state_vcf=ancestral_state.upper()
            derived_state_vcf=derived_state.upper()

            indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
            if (indelNs>-3 or indel_length==0 or ancestral_state_vcf=="N" or ancestral_state_vcf=="-" or ancestral_state_vcf=="n"):
                print("N found in indel")
                return('error')
            else:
                print(indel_seq_annotated)
                print("after insertion:")
                print("MH|INS|SD|spacer seq between two patterns|MH|INS|SD: ")
                print(seq_after_insertion)
                return(chrom,pos,ancestral_state_vcf,derived_state_vcf,indel_length)

def fa2SD_direct_deletion(refFA,chrom,pos,MH_length,SD_motif_length,max_distance,pre_SD_length=0):
    """Function to create scars of a SD-direct_deletion, no matter if trans or loop-out.
    It occur when Z>0 and x = 0
    DSB repair on a fasta file giving outuput in vcf-like entries.
    First and last 2kbp of chromosomes should not be used. 
    MH_seq = seq of MH
    SD_motif= seq of SD
    i_MH=position of MH from the start of seq window
    Arguments are:
    refFA: a reference fasta file in pysam format
    chr: the chromosome to be mutated
    pos: the start (0-based) position of the DSB
    MH_length: length of Microhomology motif
    SD_motif_length: length of SD motif
    max_distance: starts from end deletion of step1
    pre_SD_length: standard SD-loopouts have SD starting right after Microhomology motif. By setting this to >0 one can change that.
    Output has indel_length in addition to vcf-like annotation.
    i converted anc into der e der rather than insertion, here the pos is the position of the seqeunce that mediates the deletion of sthe seq between MH2 and SD2
    """
    MH_seq=refFA.fetch(chrom,pos-MH_length,pos)
    SD_motif=refFA.fetch(chrom,pos+pre_SD_length,pos+pre_SD_length+SD_motif_length)
    seq_window=refFA.fetch(chrom,pos+SD_motif_length+pre_SD_length,pos+max_distance-1) #from SD to maxdistance
    seq_window_fabio=refFA.fetch(chrom,pos-MH_length,pos+max_distance - 1) 
    SD_motif_fabio=refFA.fetch(chrom,pos+pre_SD_length,pos+pre_SD_length+SD_motif_length)
    i_MH=seq_window.find(MH_seq) #=position of MH from SD in seq window
    #pos of the sequence MH|SD that mediates the deletion of sequence where occur the DSB
    print("seq window from 1° MH to max distance: "+seq_window_fabio)
    print("MH seq: "+MH_seq)
    print("SD motif: "+SD_motif)
    print("position of 2°MH from 1°SD,(spacer seq),(i_MH): "+str(i_MH)) 
    #the 1° MH is always at 0 pos in seq_window_fabio
#    print("MH++SD+..+MH: "+refFA.fetch(chr,pos-MH_length,pos+i_MH+MH_length))
    if i_MH < 0:
        print("no MH possible in range")
        return("error")
    else:
        i_MH=i_MH+SD_motif_length+pre_SD_length
        #here i_MH is the pos of the last nt of 2° SD from 2° MH
        seq_window=refFA.fetch(chrom,pos+i_MH+MH_length,pos+max_distance)
        #up here seq_window is from 2°MH to max distance, it's used to find 2°SD 
        i_SD=seq_window.find(SD_motif) #position of 2° SD from the end of 2° MH
        if i_SD < 0:
            print("no SD possible in range")
            return("error") 
        else:
            i_SD=i_SD+i_MH+MH_length
            #up here i_sd is the position of 2° SD from the last nt of 1° MH 
            indel_length=i_SD-i_MH-MH_length
            #indel= distance from 2° SD and the 1° SD - distance from 2° SD and 2°MH - length MH
            print("distance between 2°MH and 2° SD == length deletion : "+str(indel_length))
            print("MH|SD|spacer seq between two patterns|MH|DEL|SD: ") #+refFA.fetch(chr,pos-MH_length,pos+i_MH+MH_length+i_SD+SD_motif_length))
            indel_seq_annotated=refFA.fetch(chrom,pos-MH_length,pos,)+"|"+refFA.fetch(chrom,pos,pos+SD_motif_length) +"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_MH)+"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)+"|"+refFA.fetch(chrom,pos+i_MH+MH_length,pos+i_SD) +"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length)
            indel_seq=refFA.fetch(chrom,pos+i_MH+MH_length,pos+i_SD)
            #indel_seq = from pos, + distance between 2° SD and 2° MH + length of MH to 2°SD from the the last nt of 2° MH
            seq_after_deletion=refFA.fetch(chrom,pos-MH_length,pos,)+"|"+refFA.fetch(chrom,pos,pos+SD_motif_length) +"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_MH)+"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)+"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length) 
            derived_state=refFA.fetch(chrom,pos-1,pos)
            derived_state_vcf=derived_state.upper()
            ancestral_state=refFA.fetch(chrom,pos,pos)+indel_seq
            ancestral_state_vcf=ancestral_state.upper()
            indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
            if (indelNs>-3 or indel_length==0 or ancestral_state_vcf=="N" or ancestral_state_vcf=="-" or ancestral_state_vcf=="n"):
                print("N found in indel")
                return('error')
            else:
                print(indel_seq_annotated)
                print("after deletion:")
                print("MH|SD|spacer seq between two patterns|MH|SD: ")
                print(seq_after_deletion)
                return(chrom,pos,ancestral_state_vcf,derived_state_vcf,indel_length)

def fa2SD_direct_substitution(refFA,chrom,pos,MH_length,SD_motif_length,max_distance,Z_length): 
    """Function to create scars of a SD-direct_substitution, no matter if trans or loop-out.
    It occur when Z>0 and X > 0
    DSB repair on a fasta file giving outuput in vcf-like entries.
    First and last 2kbp of chromosomes should not be used. 
    MH_seq = seq of MH
    SD_motif= seq of SD
    Z_seq = seq between MH and SD(replaced by X)

    i_MH=position of MH from the start of seq window
    Arguments are:
    refFA: a reference fasta file in pysam format
    chr: the chromosome to be mutated
    pos: the start (0-based) position of the DSB
    MH_length: length of Microhomology motif
    SD_motif_length: length of SD motif
    max_distance: starts from end deletion of step1
    Z_length: length of the seqeunce between MH and SD, is the replaced sequence. 
    By setting this to >0 one can change that.
    Output has indel_length in addition to vcf-like annotation.
    """

    #fabio: devo ben definire la lunghezza dell'indel come dice gemini, o dato che prendo le prima 4 colonne non serve fare ciò?

    MH_seq=refFA.fetch(chrom,pos-MH_length,pos) 
    SD_motif=refFA.fetch(chrom,pos+Z_length,pos+Z_length+SD_motif_length) 
    #Z_length=pre_SD_length
    Z_seq=refFA.fetch(chrom,pos,pos+Z_length)
    seq_window=refFA.fetch(chrom,pos+Z_length+SD_motif_length,pos+max_distance-1) #from SD to maxdistance
    seq_window_fabio=refFA.fetch(chrom,pos-MH_length,pos+max_distance - 1) 
    #SD_motif_=refFA.fetch(chrom,pos+pre_SD_length,pos+pre_SD_length+SD_motif_length)
    i_MH=seq_window.find(MH_seq) #=position of 2°MH from 1°SD in seq window
    print("seq window from 1° MH: "+seq_window_fabio)
    print("MH seq: "+MH_seq)
    print("SD motif: "+SD_motif)
    print("Z_seq: " + Z_seq)
    print("distance between 1°SD and 2° MH == Y (i_MH): "+str(i_MH))
#    print("MH++SD+..+MH: "+refFA.fetch(chr,pos-MH_length,pos+i_MH+MH_length))
#the 1° MH is always at 0 pos in seq_window_fabio
    if i_MH < 0:
        print("no MH possible in range")
        return("error")
    else:
        i_MH=i_MH+SD_motif_length+Z_length
        #here i_MH is the pos of the last nt of 2° SD from 2° MH
        seq_window=refFA.fetch(chrom,pos+i_MH+MH_length,pos+max_distance)
        #up here seq_window is from 2°MH to max distance, it's used to find 2°SD 
        i_SD=seq_window.find(SD_motif) #position of 2° SD from the end of 2° MH
        print("position of 2° SD from the end of 2° MH:  "+str(i_SD))
        if i_SD < 0:
            print("no SD possible in range")
            return("error") 

        else:
            i_SD=i_SD+i_MH+MH_length
            #up here i_sd is the position of 2° SD from the last nt of 1° MH 
            indel_length=i_SD-i_MH-MH_length
            #indel= distance from 2° SD and the 1° SD - distance from 2° SD and 2°MH - length MH
            print("distance between 2°MH print and 2° SD == length of X : "+str(indel_length))
            print("MH|Z|SD|spacer seq between two patterns|MH|X|SD: ") #+refFA.fetch(chr,pos-MH_length,pos+i_MH+MH_length+i_SD+SD_motif_length))
 #da qui in giù devi sistemarlo eh

            indel_seq_annotated=refFA.fetch(chrom,pos-MH_length,pos)+"|"+refFA.fetch(chrom,pos,pos+Z_length)+"|"+refFA.fetch(chrom,pos+Z_length,pos+Z_length+SD_motif_length) +"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_MH)+"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)+"|"+refFA.fetch(chrom,pos+i_MH+MH_length,pos+i_SD) +"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length)
            #X seq that replaces Z
            indel_seq=refFA.fetch(chrom,pos+i_MH+MH_length,pos+i_SD)# indel_seq = DER
            #Z that is replaced
            Z_ANC_seq=refFA.fetch(chrom,pos,pos+Z_length) #Z_ANC_seq = ANC
            nt_anchor=refFA.fetch(chrom,pos-1,pos)
            seq_after_deletion=refFA.fetch(chrom,pos-MH_length,pos,)+"|"+refFA.fetch(chrom,pos+i_MH+MH_length,pos+i_SD)+"|"+refFA.fetch(chrom,pos+Z_length,pos+Z_length+SD_motif_length) +"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_MH)+"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)+"|"+refFA.fetch(chrom,pos+i_MH+MH_length,pos+i_SD) +"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length)

            
            anc, der = normalize_vcf(Z_ANC_seq,indel_seq)
            #print(f"POS={pos}, ANC={anc}, DER={der}")
            
            #i put the nt anchor first cause the vcf
            derived_state_vcf=nt_anchor.upper()+der.upper()
            ancestral_state_vcf=nt_anchor.upper()+anc.upper()

            indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
            if (indelNs>-3 or indel_length==0 or ancestral_state_vcf=="N" or ancestral_state_vcf=="-" or ancestral_state_vcf=="n"):
                print("N found in indel")
                return('error')
            else:
                print(indel_seq_annotated)
                print("after substitution:")
                print("MH|SUB(X)|SD|spacer seq between two patterns|MH|X|SD: ")
                print(seq_after_deletion)
                #chrom,pos,ref,alt
                return(chrom,pos,ancestral_state_vcf,derived_state_vcf,indel_length)
#devi fare questa
def fa2SD_snapback_substitution(refFA,chrom,pos,MH_length,SD_motif_length,max_distance,Z_length): 
    """Function to create scars of a SD_inverted_substitution, thanks to snapback mechanism.
    It occur when Z>0 and X > 0
    DSB repair on a fasta file giving outuput in vcf-like entries.
    First and last 2kbp of chromosomes should not be used. 
    MH_seq = seq of MH
    SD_motif= seq of SD
    Z_seq = seq between MH and SD (replaced by X)

    i_MH=position of MH from the start of seq window
    Arguments are:
    refFA: a reference fasta file in pysam format
    chr: the chromosome to be mutated
    pos: the start (0-based) position of the DSB
    MH_length: length of Microhomology motif
    SD_motif_length: length of SD motif
    max_distance: starts from end deletion of step1
    pre_SD_length: standard SD-loopouts have SD starting right after Microhomology motif.
    By setting this to >0 one can change that.
    Output has indel_length in addition to vcf-like annotation.
    """

    #fabio: 
    MH_seq=refFA.fetch(chrom,pos-MH_length,pos) 
    SD_motif=refFA.fetch(chrom,pos+Z_length,pos+Z_length+SD_motif_length) 
    Z_seq=refFA.fetch(chrom,pos,pos+Z_length)
    seq_window=refFA.fetch(chrom,pos+Z_length+SD_motif_length,pos+max_distance-1) #from SD to maxdistance
    seq_window_fabio=refFA.fetch(chrom,pos-MH_length,pos+max_distance - 1) 
    
    SD_revc=seq2revcomplement(SD_motif)
    i_SD=seq_window.find(SD_revc)+Z_length 
    print("seq window from 1° MH: "+seq_window_fabio)
    print("MH seq: "+MH_seq)
    print("SD motif: "+SD_motif)
    print("Z_seq: " + Z_seq)
    print("distance between 1°SD and 2° SD == Y (i_MH): "+str(i_SD))
#    print("MH++SD+..+MH: "+refFA.fetch(chr,pos-MH_length,pos+i_MH+MH_length))
#the 1° MH is always at 0 pos in seq_window_fabio
    if i_SD < 0:
        print("no SD possible in range")
        return("error")
    else:
        i_SD=i_SD+SD_motif_length
        #here i_MH is the pos of the last nt of 2° SD from 1° SD
        seq_window=refFA.fetch(chrom,pos+i_SD+SD_motif_length,pos+max_distance)
        #up here seq_window is from 2°SD to max distance, it's used to find 2°MH or mhrevc
        MH_revc=seq2revcomplement(MH_seq)
        i_MH=seq_window.find(MH_revc) #position of 2° SD from the end of 2° MH
        if i_MH < 0:
            print("no SD possible in range")
            return("error") 

        else:
            i_MH=i_MH+i_SD+SD_motif_length
            #i_MH distance between 2° SD adn 2°MH  
            indel_length=i_MH-i_SD-SD_motif_length
            #indel= distance from 2° SD and the 1° SD - distance from 2° SD and 2°MH - length MH
            print("distance between 2°MH print and 2° SD == length substitution : "+str(indel_length))
            print("MH|Z|SD|spacer seq between two patterns|SDrevc|Xrevc|MHrevc: ") #+refFA.fetch(chr,pos-MH_length,pos+i_MH+MH_length+i_SD+SD_motif_length))
 #da qui in giù devi sistemarlo eh
            
            indel_seq_annotated=refFA.fetch(chrom,pos-MH_length,pos)+"|"+refFA.fetch(chrom,pos,pos+Z_length)+"|"+refFA.fetch(chrom,pos+Z_length,pos+Z_length+SD_motif_length) +"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_SD)+"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length)+"|"+refFA.fetch(chrom,pos+i_SD+SD_motif_length,pos+i_MH) +"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)
            
            #X revc
            X_revc=refFA.fetch(chrom,pos+i_SD+SD_motif_length,pos+i_MH)
            #but the revc of Xrevc replaces Z
            indel_seq=seq2revcomplement(X_revc)
            
            #Z that is replaced
            Z_ANC_seq=refFA.fetch(chrom,pos,pos+Z_length)
            
            nt_anchor=refFA.fetch(chrom,pos-1,pos)

            seq_after_substitution=refFA.fetch(chrom,pos-MH_length,pos,)+"|"+str(indel_seq)+"|"+refFA.fetch(chrom,pos+Z_length,pos+Z_length+SD_motif_length)+"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_SD)+"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length)+"|"+refFA.fetch(chrom,pos+i_SD+SD_motif_length,pos+i_MH) +"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)
            
            anc, der = normalize_vcf(Z_ANC_seq,indel_seq)

            #derived_state_vcf=nt_anchor.upper()+indel_seq.upper()# i put the anchor for the happyness of vcf
            #ancestral_state_vcf=nt_anchor.upper()+Z_ANC_seq.upper()

            ancestral_state_vcf=nt_anchor.upper()+anc.upper()
            derived_state_vcf=nt_anchor.upper()+der.upper()

            indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
            if (indelNs>-3 or indel_length==0 or ancestral_state_vcf=="N" or ancestral_state_vcf=="-" or ancestral_state_vcf=="n"):
                print("N found in indel")
                return('error')
            else:
                print(indel_seq_annotated)
                print("after substitution:")
                print("MH|SUB(X)|SD|spacer seq between two patterns|MHrevc|Xrevc|SDrevc: ")
                print(seq_after_substitution)
                #chrom,pos,ref,alt
                return(chrom,pos,ancestral_state_vcf,derived_state_vcf,indel_length)




#fastafile='GWHAAEV00000000.1.genome.fasta.gz'
fastafile=args['ref']
outfile=args['out']
chrom=args['chrom']
pos=args['pos']
nsims=int(args['nsims'])
xmechanism=args['mechanism']
xMHlength=args['MHlength']
xSDlength=args['SDlength']
xmaxdistance=args['maxdistance']
#inserito da fabio
xZ_length=args['Z_length']

refFA=FastaFile(fastafile)

myvcf=list()
while len(myvcf)<nsims:
    if (chrom==None):
        xchr=random.choices(refFA.references,weights=refFA.lengths,k=1)[0]
    else:
        xchr=chrom
    if (pos==None):
        xpos=random.sample(range(2000,refFA.lengths[refFA.references.index(xchr)]-2000),1)[0]
    else:
        xpos=int(pos)
    #print("create indel")
    print(xchr)
    print(xpos)
    if xmechanism=="deletion_MMEJ":
        res=fa2deletion_MMEJ(refFA,xchr,xpos,xMHlength,xmaxdistance)
    elif xmechanism=="deletion_NHEJ":
        res=fa2deletion(refFA,xchr,xpos,xMHlength,xmaxdistance)
    elif xmechanism=="SD_inverted_insertion":
        res=fa2insertion_snapback(refFA,xchr,xpos,MH_length=xMHlength,SD_motif_length=xSDlength,max_distance=xmaxdistance)
    elif xmechanism=="SD_inverted_deletion":
        res=fa2deletion_snapback(refFA,xchr,xpos,MH_length=xMHlength,SD_motif_length=xSDlength,max_distance=xmaxdistance)
    elif xmechanism=="insertion":
        res=fa2insertion(refFA,xchr,xpos,xMHlength,xmaxdistance)
    elif xmechanism=="SD_direct_insertion":
        res=fa2SD_direct_insertion (refFA,xchr,xpos,xMHlength,xSDlength,xmaxdistance)
    elif xmechanism=="SD_direct_deletion":
        res=fa2SD_direct_deletion (refFA,xchr,xpos,xMHlength,xSDlength,xmaxdistance)
    elif xmechanism=="SD_direct_substitution":
        res=fa2SD_direct_substitution (refFA,xchr,xpos,xMHlength,xSDlength,xmaxdistance,xZ_length)
    elif xmechanism=="SD_inverted_substitution":
        res=fa2SD_snapback_substitution (refFA,xchr,xpos,xMHlength,xSDlength,xmaxdistance,xZ_length)

    if res != "error":
        myvcf.append(res)
        print(res)
    else: 
        print("error")
        break    # questo era il motivo per cui andava in loop

myvcf = pd.DataFrame(myvcf, columns=['#CHR', 'POS','REF','ALT',"INDEL_LENGTH"])
#myvcf = myvcf[['#CHR', 'POS','REF','ALT']]
if outfile!=None:
    myvcf.to_csv(outfile, sep="\t",index=None)

#chr="GWHAAEV00000001.1"
#pos=5000
#MH_length1=2
#MH_length2=2
#max_distance_MMEJ1=100
#distance_MMEJ2=10

