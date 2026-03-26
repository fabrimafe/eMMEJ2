import pandas as pd
from pysam import FastaFile
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import argparse

# Construct an argument parser
all_args = argparse.ArgumentParser()

# Add arguments to the parser
all_args.add_argument("-v", "--vcf", required=True,
   help="input vcf")
all_args.add_argument("-r", "--ref", required=True,
   help="reference genome in fasta format")
all_args.add_argument("-o", "--outputfile", required=True,
      help="output file in vcf-like format")

args = vars(all_args.parse_args())


#xdata=pd.read_csv('test6.vcf', sep='\t') #, names=["chr", "pos","REF","ALT"])
#xdata=pd.read_csv('test8.vcf', sep='\t') #, names=["chr", "pos","REF","ALT"])
xdata=pd.read_csv(f'{args["vcf"]}.vcf' , sep='\t')
xdata.columns=["chr","pos","REF","ALT", 'ground_true']

fastafile=args['ref'] #'GWHAAEV00000000.1.genome.fasta.gz'
size_window=2000
#refFA=DSBsimulate.FastaFile(fastafile)
refFA=FastaFile(fastafile)

def alignments2vcf(xalignment,chr='CHR',pos=0,starting_base="N"):
    """convert a single alignment to vcf-like format, adding the original starting position as a tag for the alignment (useful in case of many alignments)"""
    genotypes = []
    inindel=0
    xindel=''
    xpos=pos
    ref_genotype="N"
    for a, b in zip(xalignment[0]+'n',xalignment[1]+'n'):
        if a==b and inindel==1:
            genotypes.append([chr,pos_indel,ref_genotype_var,xindel,pos])
            inindel=0
        if a != b and a!="-" and b!="-":
            genotypes.append([chr,xpos,a,b,pos])
        if a != b and a=='-' and b!='-':
            if inindel==0:
                ref_genotype_var=ref_genotype
                xindel=ref_genotype_var
                pos_indel=xpos-1
            inindel=1
            xindel=xindel+b
        if a != b and a!='-' and b=='-':
            if inindel==0:
                pos_indel=xpos
                xindel=ref_genotype
                ref_genotype_var=ref_genotype
                pos_indel=xpos-1
            inindel=1
            ref_genotype_var=ref_genotype_var+a
        ref_genotype=a
        xpos=xpos+1
    return(genotypes)

def vcf2realignedvcfs(refA,chr,pos,REF,ALT,length_around):
    """perform indel realignment for a single variant in vcf format listing all possible equally-best alignments. Then it applies alignments2vcf for each of these, returning all possible
    indel calls for a single variant. Output in vcf format + an identifier tag that indicates original_position.index, where the index specifies the index of the alignment of that position (useful when a single indel can be split
    in several sub-indels)"""
    if REF!=ALT and REF!="N" and ALT!="N" and not pd.isna(REF) and not pd.isna(ALT):
        starting_base=refFA.fetch(chr,pos-2-length_around,pos)
        seqREF=refFA.fetch(chr,pos-1-length_around,pos-1)+REF+refFA.fetch(chr,pos+len(REF)-1,pos+length_around+len(REF))
        seqALT=refFA.fetch(chr,pos-1-length_around,pos-1)+ALT+refFA.fetch(chr,pos+len(REF)-1,pos+length_around+len(REF))
        alignments = pairwise2.align.localms(seqREF, seqALT ,5, -1, -0.5, -0.1)
        scores=[]
        for a in alignments:
            # print(format_alignment(*a))
            scores.append(a[2])
        vcfs=[]
        al_counter=0
        for a in alignments:
            if a[2]==max(scores):
                al_counter=al_counter+1
                if a[0][0]!="-" and a[0][len(a[0])-1]!="-":
                    vcfs.append(alignments2vcf(a,chr,pos-1-length_around)[0]) #,starting_base
                    #fix tag to have the original position of the vcf
                    vcfs[len(vcfs)-1][4]=str(vcfs[len(vcfs)-1][4]+length_around+1)+"."+str(al_counter)
                    #for ia in vcfs:
                    #    a[4]=pos
    else:
            vcfs=[]
    return(vcfs)

def flatten_2list(list_of_lists):
    flat_list=[]
    for item in list_of_lists:
        for item2 in item:
            flat_list.append(item2)
    return(flat_list)

        #alignments=pd.DataFrame(alignments, columns = ['seqREF','seqALT','score',"x","length"])
        #alignments=alignments[alignments["score"] == max(alignments["score"])]
        #return(alignments)
#i=0
#for i in range(10000):
#    print(i)
#    vcfs=vcf2realignedvcfs(refFA,xdata['chr'][i],xdata['pos'][i], xdata['REF'][i],xdata['ALT'][i],20)

#i=661
#refFA
#chr=xdata['chr'][i]
#pos=xdata['pos'][i]
#REF=xdata['REF'][i]
#ALT=xdata['ALT'][i]
#length_around=20
#pd.DataFrame(vcfs,columns = ['#chr','pos','REF',"ALT","original_pos"])

df=xdata.apply(lambda row : vcf2realignedvcfs(refFA,row['chr'],row['pos'], row['REF'],row['ALT'],40), axis = 1)
df=flatten_2list(df.tolist())
df=pd.DataFrame(df,columns = ['#chr','pos','REF',"ALT","original_pos"])

df.loc[:, 'pos'] = df.loc[:, 'pos'] + 1
df['N'] = df.apply(lambda x: 'N' in x['REF'], axis=1)
df = df.loc[(df['N'] == False), :]
df.drop(columns=['N'], inplace=True)
df.rename(columns={'#chr': 'CHR', 'pos': 'POS'}, inplace=True)
df.to_csv(args['outputfile'], sep="\t",index=None)

