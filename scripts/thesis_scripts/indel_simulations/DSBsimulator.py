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

def fa2deletion_MMEJ(refFA,chrom,pos,MH_length,max_distance_MMEJ):
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

def fa2insertion_SDloopout(refFA,chrom,pos,MH_length,SD_motif_length,max_distance,pre_SD_length=0):
    """Function to create scars of a SD-loop-out DSB repair on a fasta file giving outuput in vcf-like entries. First and last 2kbp of chromosomes should not be used.
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
    i_MH=seq_window.find(MH_seq)
    print("seq window: "+seq_window)
    print("MH seq: "+MH_seq)
    print("SD motif: "+SD_motif)
    print("position MH (i_MH): "+str(i_MH))
#    print("MH++SD+..+MH: "+refFA.fetch(chr,pos-MH_length,pos+i_MH+MH_length))
    if i_MH < 0:
        print("no MH possible in range")
        return("error")
    else:
        i_MH=i_MH+SD_motif_length+pre_SD_length
        seq_window=refFA.fetch(chrom,pos+i_MH+MH_length,pos+max_distance)
        i_SD=seq_window.find(SD_motif)
        if i_SD < 0:
            print("no SD possible in range")
            return("error")
        else:
            i_SD=i_SD+i_MH+MH_length
            indel_length=i_SD-i_MH-MH_length
            print("position SD (i_SD) == length insertion: "+str(indel_length))
            print("MH|SD|loop|MH|INS|SD: ") #+refFA.fetch(chr,pos-MH_length,pos+i_MH+MH_length+i_SD+SD_motif_length))
            indel_seq_annotated=refFA.fetch(chrom,pos-MH_length,pos,)+"|"+refFA.fetch(chrom,pos,pos+SD_motif_length) +"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_MH)+"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)+"|"+refFA.fetch(chrom,pos+i_MH+MH_length,pos+i_SD) +"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length)
            indel_seq=refFA.fetch(chrom,pos+i_MH+MH_length,pos+i_SD)
            print(indel_seq_annotated)
            ancestral_state_vcf=refFA.fetch(chrom,pos-1,pos)
            derived_state_vcf=refFA.fetch(chrom,pos-1,pos)+indel_seq
            indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
            if (indelNs>-3 or indel_length==0 or ancestral_state_vcf=="N" or ancestral_state_vcf=="-" or ancestral_state_vcf=="n"):
                print("N found in indel")
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
    SD_revc=seq2revcomplement(SD_motif)
    i_SD=seq_window.find(SD_revc)+pre_SD_length
    print("seq window: "+seq_window)
    print("MH seq: "+MH_seq)
    print("SD motif: "+SD_motif)
    print("position 2nd SD motif (i_SD): "+str(i_SD))
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
            print("position MH (i_MH) == length insertion: "+str(indel_length))
            print("MH|SD|snapback-loop|SDrevc|INSrevc|MHrevc: ")#refFA.fetch(chr,pos-MH_length,pos+i_MH+MH_length))
            indel_seq_annotated=refFA.fetch(chrom,pos-MH_length,pos)+"|"+refFA.fetch(chrom,pos,pos+SD_motif_length) +"|"+refFA.fetch(chrom,pos+SD_motif_length,pos+i_SD)+"|"+refFA.fetch(chrom,pos+i_SD,pos+i_SD+SD_motif_length)+"|"+refFA.fetch(chrom,pos+i_SD+SD_motif_length,pos+i_MH) +"|"+refFA.fetch(chrom,pos+i_MH,pos+i_MH+MH_length)
            indel_seq=seq2revcomplement(refFA.fetch(chrom,pos+i_SD+SD_motif_length,pos+i_MH))
            indelNs=indel_seq.find('N') + indel_seq.find('n') + indel_seq.find('-')
            ancestral_state_vcf=refFA.fetch(chrom,pos-1,pos)
            derived_state_vcf=refFA.fetch(chrom,pos-1,pos)+indel_seq
            if (indelNs>-3 or indel_length==0 or ancestral_state_vcf=="N" or ancestral_state_vcf=="-" or ancestral_state_vcf=="n"):
                print("N found in indel")
                return('error')
            else:
                print(indel_seq_annotated)
                return(chrom,pos,ancestral_state_vcf,derived_state_vcf,indel_length)

def fa2insertion_MMEJtrans(refFA,chrom,pos,MH_length1,MH_length2,max_distance_MMEJ,distance_MMEJ2,MMEJ2_beforedeletion=0):
    """Function to create deletions on a fasta file giving outuput in vcf-like format. First and last 2kbp of chromosomes should not be used.
    Arguments are:
    refFA: a reference fasta file in pysam format
    chr: the chromosome to be mutated
    pos: the start (0-based) position of the mutation, i.e. if you want to delete position 121 and 122 you need to set 120
    MH_length: length of Microhomology motif
    max_distance_MMEJ: max distance between cut site and Microhomology motif
    distance_MMEJ2: starts from end deletion of step1
    The output shows:
    chr: chromosome
    pos: position
    the ancestral_state;
    the derived_state after repair
    the size of the deletion after the first MMEJ step;
    the size of the insertion caused by template switching;
    an alignment file of the 20bp region surrounding the indel in Biopython format
    """
    print("Step1: first annealing and MMEJ DSB repair")
    step1=fa2deletion_MMEJ(refFA,chrom,pos,MH_length1,max_distance_MMEJ)
    #step2: a) the DNA is partially repaired (elongated) on one strand and then dissociated.
    #       b) a microhomology motif is found on this elongated strand.
    #       c) a new MMEJ repair occurs at this second motif
    if step1=="error":
        print("no MH1 possible in range")
        return("error")
    else:
        deletion1_length=step1[4]
        print("Step2: template switching and second MMEJ DSB repair")
        MH_seq2=refFA.fetch(chrom,pos+deletion1_length+distance_MMEJ2,pos+deletion1_length+distance_MMEJ2+MH_length2)
        seq_window=refFA.fetch(chrom,pos+deletion1_length-MMEJ2_beforedeletion,pos+deletion1_length+distance_MMEJ2)
        i_MH2=seq_window.find(MH_seq2) #this is necessary to find position on reverse strand
        if i_MH2 < 0:
            print("no MH2 possible in range")
            return("error")
        else:
            insertion_length=distance_MMEJ2+MH_length2
    #note that if max_distance_MMEJ<MH_length1+distance_MMEJ2, MH2 will be found on reverse at positions where for sure it must be single strand, i.e. before deletion1_length+MH_length1+distance_MMEJ2
    #in these cases part of deletion1:
            insertion_seq=refFA.fetch(chrom,pos+deletion1_length,pos+deletion1_length+distance_MMEJ2)
            deletion_seq=step1[2]
            #print("DEL1: ",deletion_seq)
            #print("INS: ",insertion_seq)
            ancestral_state=refFA.fetch(chrom,pos,pos+deletion1_length)
            derived_state=refFA.fetch(chrom,pos-1,pos)+refFA.fetch(chrom,pos+deletion1_length,pos+deletion1_length-MMEJ2_beforedeletion+i_MH2+insertion_length)
            derived_state_vcf=refFA.fetch(chrom,pos-1,pos)+derived_state+refFA.fetch(chrom,pos+deletion1_length,pos+deletion1_length+insertion_length)
            #print("DEL1|not-INS/beforeMH2-on-reverse|MH2onrev|template.onreverse-for-INS")
            #print(refFA.fetch(chrom,pos-1,pos+deletion1_length)+"|"+refFA.fetch(chrom,pos+deletion1_length,pos+deletion1_length+i_MH2)+"|"+refFA.fetch(chrom,pos+deletion1_length+i_MH2,pos+deletion1_length+i_MH2+MH_length2)+"|"+refFA.fetch(chrom,pos+deletion1_length+i_MH2+MH_length2)
            #,pos+deletion1_length-MMEJ2_beforedeletion+i_MH2)
            ancestral_state_vcf=refFA.fetch(chrom,pos-1,pos+deletion1_length+insertion_length)
            derived_state_around=refFA.fetch(chrom,pos-20,pos)+derived_state+refFA.fetch(chrom,pos+deletion1_length,pos+deletion1_length+insertion_length+20)
            ancestral_state_around=refFA.fetch(chrom,pos-20,pos+deletion1_length+insertion_length+20)
            alignments = pairwise2.align.globalms(ancestral_state_around, derived_state_around ,5, -1, -0.5, -0.1)
            for a in alignments:
                print(format_alignment(*a))
            return(chrom,pos,ancestral_state_vcf,derived_state_vcf,deletion1_length,insertion_length-MH_length2,alignments)


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
    elif xmechanism=="SDloopout":
        res=fa2insertion_SDloopout(refFA,xchr,xpos,xMHlength,xSDlength,xmaxdistance)
    elif xmechanism=="SDsnapback":
        res=fa2insertion_snapback(refFA,xchr,xpos,MH_length=xMHlength,SD_motif_length=xSDlength,max_distance=xmaxdistance)
    elif xmechanism=="insertion":
        res=fa2insertion(refFA,xchr,xpos,xMHlength,xmaxdistance)
    elif xmechanism=="MMEJtrans":
        res=fa2insertion_MMEJtrans(refFA,chrom,pos,xMHlength,xMHlength,xmaxdistance,distance_MMEJ2=xSDlength)
    if res != "error":
        myvcf.append(res)
        print(res)

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

