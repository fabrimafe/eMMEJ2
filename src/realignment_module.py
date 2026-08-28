import pandas as pd
from pysam import FastaFile
#from Bio import pairwise2
#from Bio.pairwise2 import format_alignment
import Bio
from Bio import Align
import numpy as np
from Bio import pairwise2
import sys

sys.path.append("/home/fabio/git/eMMEJ2/scripts")
from sliding_vcf import sliding_vcf
#import argparse

#from /home/fabio/git/eMMEJ2/scripts/DSBsimulatorfabio_2.py import sliding_vcf 

def alignments2vcf(xalignment, chrom='CHR', pos=0, starting_base="N"):
    """convert a single alignment to vcf-like format, adding the
        original starting position as a tag for the
        alignment (useful in case of many alignments)"""
    genotypes = []
    inindel=0
    xindel=''
    xpos=pos
    ref_genotype="N"
    for a, b in zip(xalignment[0]+'n',xalignment[1]+'n'):
        if a==b and inindel==1:
            genotypes.append([chrom,pos_indel,ref_genotype_var,xindel,pos])
            inindel=0
        if a != b and a!="-" and b!="-":
            genotypes.append([chrom,xpos,a,b,pos])
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

def vcf2realignedvcfs_pairwise2(refFA,chrom,pos,REF,ALT,length_around):
    """perform indel realignment for a single variant in vcf format
    listing all possible equally-best alignments. Then it applies
    alignments2vcf for each of these, returning all possible
    indel calls for a single variant. Output in vcf format + an identifier
    tag that indicates original_position.index, where the index specifies
    the index of the alignment of that position (useful when a
    single indel can be split in several sub-indels).
    It uses Bio.pairwise2 was deprecated in Biopython Release 1.80.
"""
    if REF!=ALT and REF!="N" and ALT!="N" and not pd.isna(REF) and not pd.isna(ALT):
        starting_base=refFA.fetch(chrom,pos-1-length_around,pos)
        
        seqREF=refFA.fetch(chrom,pos-1-length_around,pos-1)+REF+refFA.fetch(chrom,pos+len(REF)-1,pos+length_around+len(REF))
        seqALT=refFA.fetch(chrom,pos-1-length_around,pos-1)+ALT+refFA.fetch(chrom,pos+len(REF)-1,pos+length_around+len(REF))
        alignments = pairwise2.align.localms(seqREF, seqALT ,5, -1, -0.5, -0.1)
        scores=[]
        for a in alignments:
            scores.append(a[2])
        print(len(alignments))
        vcfs=[]
        al_counter=0
        for a in alignments:
            if a[2]==max(scores):
                al_counter=al_counter+1
                if a[0][0]!="-" and a[0][len(a[0])-1]!="-":
                    vcfs.append(alignments2vcf(a,chrom,pos-1-length_around)[0]) #,starting_base
                    #fix tag to have the original position of the vcf
                    vcfs[len(vcfs)-1][4]=str(vcfs[len(vcfs)-1][4]+length_around+1)+"."+str(al_counter)
    else:
            vcfs=[]
    return(vcfs)

def flatten_2list(list_of_lists):
    flat_list=[]
    for item in list_of_lists:
        for item2 in item:
            flat_list.append(item2)
    return(flat_list)


def DSB_background_generator(refFA, chrom, pos, ref, alt, indel_type, max_distance):
    """
    creating 'background' variants shifting the pos upstream and downsteam 


        refFA - genoma di riferimento
        chrom - cromosoma (str)
        pos -  posizione originale (1-based, VCF)
        ref - REF originale (stringa, con anchor) senza ancora nelle sub
        alt - ALT originale (stringa, con anchor) senza ancora nelle sub
        indel_type - 'DEL' o 'INS' o 'SUB'
        max_distance - numero massimo di basi di shift (int)

    Return:
    chrom, new_pos,variant_id, new_anc, new_der,pos,indel_type
    """
    results = []
    ref_len = len(ref)
    alt_len = len(alt)

    for n in range(0, max_distance + 1):
        for direction in (-1, 1):          # -1 = sinistra, +1 = destra
            new_pos = pos + direction * n
            new_variant_id = f"{chrom}_{pos}"
            
            if indel_type == 'DEL':
                # ri-pesco l'intera regione (anchor + basi delete) dal riferimento
                # pysam.fetch è 0-based half-open -> new_pos-1 come start
                new_anc = (refFA.fetch(chrom, new_pos - 1, new_pos - 1 + ref_len)).upper()
                new_der = new_anc[0].upper()        # solo l'anchor

            elif indel_type == 'INS':
                # dalle inserzione anc è una nuova ancora che dipende dalla nuova pos 
                new_anc = (refFA.fetch(chrom, new_pos - 1, new_pos)).upper()
                inserted_seq = alt[1:]      # la parte inserita resta invariata
                new_der = (new_anc + inserted_seq).upper() #e poi le sommo dando vita a new_der
            elif indel_type == 'SUB':
                anc = (refFA.fetch(chrom,new_pos, new_pos + ref_len)).upper() #come nelle delezioni cambia in base alla new_pos
                der=alt.upper() #e come nelle ins rimane invariata ma senza ancora sta volta
                new_anc,new_der,NT_removed_suffix, NT_removed_prefix, total_NT_removed = sliding_vcf (anc,der)
                new_pos=new_pos+NT_removed_prefix
            else:
                raise ValueError(f"indel_type non riconosciuto: {indel_type}")

            results.append((chrom, new_pos,new_variant_id, new_anc, new_der,pos,indel_type))

    return results

def de_collapse(refFA, chrom, pos, ref, alt,de_collapse_distance):
    """
    de_collapse extend to the left and to the right ANC and DER, it adds a number of nucleotides that can't be > de_collapse_distance
    if de_collapse_distance = 2, we will have:
    0.ANC(OR DER).0
    0.ANC(OR DER).1
    0.ANC(OR DER).2
    2.ANC(OR DER).0
    2.ANC(OR DER).1
    1.ANC(OR DER).0

    """


    results = [] #apro la lista
    new_variant_id = f"{chrom}_{pos}"
    
    ref_len = len(ref)
    anc = ref
    der = alt

    for dist in range(0, de_collapse_distance + 1):
           for left in range(0, dist + 1):
                
                right = dist - left 
                
                left_ex = (refFA.fetch(chrom, pos - 1 - left, pos - 1)).upper()
                right_ex = (refFA.fetch(chrom, pos - 1 + ref_len, pos - 1 + ref_len + right)).upper()

                new_anc = left_ex + anc + right_ex
                new_der = left_ex + der + right_ex
                new_pos = pos - left

                results.append((chrom, new_pos, new_variant_id, new_anc, new_der, pos))
    
    return results
    



def vcf2realignedvcfs(refFA,chrom,pos,REF,ALT,length_around):
    """perform indel realignment for a single variant in vcf format
    listing all possible equally-best alignments. Then it applies
    alignments2vcf for each of these, returning all possible
    indcalls for a single variant. Output in vcf format + an identifier
    tag that indicates original_position.index, where the index specifies
    the index of the alignment of that position (useful when a
    single indel can be split in several sub-indels)"""
    if REF!=ALT and REF!="N" and ALT!="N" and not pd.isna(REF) and not pd.isna(ALT):
        starting_base=refFA.fetch(chrom,pos-2-length_around,pos)
        seqREF=refFA.fetch(chrom,pos-1-length_around,pos-1)+REF+refFA.fetch(chrom,pos+len(REF)-1,pos+length_around+len(REF))
        
        seqALT=refFA.fetch(chrom,pos-1-length_around,pos-1)+ALT+refFA.fetch(chrom,pos+len(REF)-1,pos+length_around+len(REF))
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score= 1.000000
        aligner.mismatch_score= -1.000000
        aligner.target_internal_open_gap_score= -2
        aligner.target_internal_extend_gap_score= -0.1
        aligner.query_internal_open_gap_score= -2
        aligner.query_internal_extend_gap_score= -0.1
        alignments = aligner.align(seqREF, seqALT)
        scores=[]
        #NB: only top scoring alignments shown with PairwiseAligner (biopython v>1.80). 
        for a in alignments:
            scores.append(a.score)
        #    print(a). 
        #print(len(alignments))
        vcfs=[]
        al_counter=0
        for a in alignments:
            if a.score==max(scores):
                al_counter=al_counter+1
                if a[0][0]!="-" and a[0][len(a[0])-1]!="-":
                    vcfs.append(alignments2vcf(a,chrom,pos-1-length_around)[0]) #,starting_base
                    #fix tag to have the original position of the vcf
                    vcfs[len(vcfs)-1][4]=str(vcfs[len(vcfs)-1][4]+length_around+1)+"."+str(al_counter)
    else:
            vcfs=[]
    return(vcfs)

def flatten_2list(list_of_lists):
    flat_list=[]
    for item in list_of_lists:
        for item2 in item:
            flat_list.append(item2)
    return(flat_list)

def set_ancestral_state_indel(df: pd.DataFrame):
    """
    Set columns for ancestral and derived alleles based on
    a binary column (ANCESTRAL) indicating whether the ancestral state is:
    REF allele - if ANCESTRAL==0
    ALT allele - if ANCESTRAL==1
    Args:
        df (pd.DataFrame) : A dataframe with the following
            columns: 'CHR', 'POS', 'REF', 'ALT', 'ANCESTRAL'

    Implementation:
    it uses numpy select, which vectorialize if-conditions (specified in condlist) and
    respective consequences (specified in choicelist)
    """
    condl = [df['ANCESTRAL'] == 0, df['ANCESTRAL'] == 1]
    choicel = [df['REF'], df['ALT']]
    df['ANC'] = np.select(condlist=condl,  choicelist=choicel)
    choicel = [df['ALT'], df['REF']]
    df['DER'] = np.select(condlist=condl,  choicelist=choicel)

