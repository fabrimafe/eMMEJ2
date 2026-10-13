""" Code Glossary:
    up -> upstream
    down -> downstream
    ref -> reference sequence
    mutant -> reference sequence with the indel embedded inside
    cand -> candidate
    dsb -> double strand break
    TS -> template_switch (relevant to insertions)
    inv_comp -> invrs and complementary sequence
"""
# importing libraries
import sys
import re
import regex as reg
# import typing
 
import pandas as pd

import numpy as np
from pysam import FastaFile
from helper import type_checker

sys.path.append('src')
from pysam_getfasta import *

def checkNaN(str):
  return str != str

class emMEJrealignment:
    def __init__(self, 
                ANC: str, DER: str,
                indel_type: str, flip: bool, include_context: bool,
                extension: bool, windowsize: int, refFA: str, MH_lengths: list, 
                chrom: str, pos_on_chr:int):
        self.indel_position = windowsize #fabri: OK, this is weird but try to understand that.
        self.indel_type = indel_type
        self.windowsize = windowsize
        self.flip = flip        
        self.include_context = include_context
        self.extension = extension
        self.refFA=refFA
        self.chrom=chrom
        self.pos_on_chr=pos_on_chr
       
        self.MH_lengths=MH_lengths 
#        if MH_lengths[0] == 0:
#            self.MH_lengths = "None"

        if (self.indel_type == 'INS'):  # cropping indel sequence to the right size.
            #print("ins")
            self.indel_length = len(DER) - 1
            self.INDEL =  DER[len(ANC):] #fabri: is this for substitutions? Then why before length-1?
            self.ANC = ANC
            self.windowsize = windowsize
        if (self.indel_type == 'DEL'): 
            #print("del")
            self.indel_length = len(ANC) - 1
            self.ANC = ANC[len(DER):]
            self.INDEL = ANC[len(DER):]
            self.windowsize = windowsize + self.indel_length
            self.indel_position = self.windowsize
        if (self.indel_type == 'SUB'):
            #print("sub")
            if len(ANC) >= len (DER):
                self.indel_length = len(ANC) - len (DER)
                self.INDEL = DER
                self.ANC = ANC
                self.DER = DER
                                
            else: 
                self.indel_length = len(DER) - len(ANC)
                self.INDEL = DER
                self.ANC = ANC
                self.DER = DER       
                
            self.windowsize = windowsize
        #print(f"DEBUG: indel_type = {repr(self.indel_type)}")
        #print(ANC)
        #print(DER)
        #print(self.pos_on_chr)
        self.ref_seq = get_ref_context(refFA=refFA, chrom=chrom,indel_pos=pos_on_chr,
                            context_window_size=self.windowsize,indel_seq=self.INDEL)
        self.ref_seq = self.ref_seq.upper()
        
        # INSERTIONS:
        #getting sequence contexts
        if (self.indel_type == 'INS'):
            self.mutant_sequence = get_mutant_context(
                        reference_contex=self.ref_seq,
                        alt=DER,
                        ref=ANC,
                        window_size=(windowsize))
            # Flipping sequences (looking upward, i.e. in 5' direction from indel, without complement)
            if self.flip:
                self.ref_seq = self.ref_seq[::-1]
                self.mutant_sequence = self.mutant_sequence[::-1]
                self.ANC = self.ANC[::-1]
                self.INDEL = self.INDEL[::-1]
                self.indel_position = self.indel_position + self.indel_length +1
        
        # DELETIONS    
        if (self.indel_type == 'DEL'):
            if self.flip:
                self.ANC = self.ANC[::-1]
                self.indel_position = self.indel_position - self.indel_length +1
                #cioè self.indel_position = a window size, in caso di deletion è windows size - la lunghezza della indel + 1
                self.ref_seq = self.ref_seq[::-1]
            #fabri: the mutant context is the alternative?
            self.mutant_sequence = get_mutant_context(
                        reference_contex=self.ref_seq,
                        alt=self.INDEL, 
                        ref=self.ANC, 
                        window_size=windowsize) 

        #SUBSTITUTION
        if (self.indel_type == 'SUB'):
            if len(ANC) >= len (DER):
                self.mutant_sequence = get_mutant_context(
                        reference_contex=self.ref_seq,
                        alt=self.INDEL,
                        ref=self.ANC,
                        window_size=windowsize)
            if len(ANC) < len (DER):
                self.mutant_sequence = get_mutant_context(
                        reference_contex=self.ref_seq,
                        alt=DER,
                        ref=ANC,
                        window_size=(windowsize))



        if self.include_context:
            self.context = {
                'ref_genome_context':self.ref_seq, 
                'mutant_sequence': self.mutant_sequence}

        
        self.tab = str.maketrans("ACTG", "TGAC") #fabri: define for complement in reverse_complement_converter. It would be nicer within reverse_complement_converter. I guess put here for faster. But then perhaps better out of the module not to do it for every row?
        self.microhomology = self.microhomology_detection()
        self.ex_data = self.export_data()


    def reverse_complement_converter(self, seq: str):
        """
        A fuction that calculate reverse complements
        """
        return seq.translate(self.tab)[::-1] #reverse

    def get_seq_updownstream(self):
        """
        get_seq_updownstream will get self as an input and will return
        the upstream and downstream sequences of both the reference and
        the mutant as an output.
        """       
        self.DSB_up = self.mutant_sequence[:(self.indel_position)] 
        self.DSB_down = self.mutant_sequence[(self.indel_position):] 
        self.ref_genome_up = self.ref_seq[:(self.indel_position)] 
        self.ref_genome_down = self.ref_seq[(self.indel_position):]


    def first_MH_finder(self):
        """ For deletion MMEJ: create a list of all possible MH motifs in ANC, i.e. deleted region:
        e.g. ANC=ACT -> [T,CT,ACT]
        Mutate:
            self.mmej_motifs (list): a list of all possible MH motifs    
        """
        mmej_motifs = pd.Series([str(self.ANC)[i:] for i in 
                              range((self.indel_length-1), -1, -1)])
        #if MH_length != "None":
              #mmej_motifs = mmej_motifs[mmej_motifs.str.len() == MH_length]
        self.mmej_motifs = mmej_motifs.to_list()
        
    def first_MH_indicator(self): 
        """
        Indicate if the 1st condition which is existence of MH sequences
        in both sides of the DSB
        Mutate:
            self.microhomology (str) : mutate to potential_INS_MMEJ
                if the 1st condition is met.
        """
        if len(self.mmej_motifs) > 0:
            self.microhomology = 'potential_INS_MMEJ'              
        
    """
    This function will detect microhomology in both insertions and deletions cases.
    it will take both reference and the mutant sequences
    as an input and will give True/False flag as an output,
    in cases of deletion mmej it will also give mmej_marked_on_ref.
    There are two conditions that must exist in order to detect
    microhomology in deletions indels:
    1. the microhomology sequence must exist both rigth before the
    indel and downstream (when looking at the reference sequence).
    2. the distance between the two microhomology sequences is equal to
    the deletion's length
    optional:
    3. the sequence around the 2nd microhomology sequence in the reference
    wont exist in the mutant's sequence
    In case of insertions, the pathway must go through the deletion step
    but with an alternative 2nd condition:
    2alt. The function will look for the following pattern -> repeat
    - microhomology - repeat.
    3. The length of the repeat must be smaller then the indel_len.
    """

    def find_subrepeat(self,x):
        """
        Function that finds the shortest (which is also the most frequent) perfectly-repeated submotif in a sequence. 
        When no submotif is repeated, this corresponds to the actual sequence. 
        The number or repeats is returned along the subrepeat.
        """
        subrepeat=x
        nrepeats_old=1
        for i in range(1,len(x)+1):
            if len(x)%i==0:
                 nrepeats=int(len(x)/i)
                 if x[:i]*nrepeats==x and nrepeats>=nrepeats_old:
                     nrepeats_old=nrepeats
                     subrepeat=x[:i]
        return([subrepeat,nrepeats_old])                 

    def microhomology_detection(self):
        self.pol_slip
        # get seqs upstream and downstream
        self.get_seq_updownstream() 
        self.microhomology = False
        self.snap_out_dict = {}
        self.loop_out_dict = {}
        self.pol_slip_dict = {}
        self.del_out_dict = {}

        if self.indel_type == 'INS': 
            self.sd_snap_back_MMEJ()
            self.sd_direct_insertions()
            self.pol_slip(xINDEL=self.INDEL)

        if self.indel_type == 'DEL': 
            self.deletion_MMEJ()
            self.pol_slip(xINDEL=self.ANC)
            self.sd_direct_deletions()
            self.sd_inverted_deletions()
        
        if self.indel_type == 'SUB':
            self.sd_direct_substitution()
            self.sd_inverted_substitution()

    def deletion_MMEJ(self):
    
        MMEJ_motif_found,MMEJ_motif_len = list(), list()

        ###1st step is to bring MH1 and P1 from the ref seq 
        for MH_lengths in self.MH_lengths:
            MH = self.ref_genome_up [ - MH_lengths :]  
            #P1 = self.ref_genome_down [(self.indel_length):(self.indel_length+MH_lengths)]
            #MH1_P1 = MH1 + P1
            print("MH: " + MH)
            print("self.INDEL :" + self.INDEL)
            print(len(self.INDEL))
            print(len(MH))
            
            #print("ex: " + str(self.extension))
            #searching the pattern in self.DSB_down
        
            if MH in self.DSB_down[(len(self.INDEL)) - len (MH):(len(self.INDEL))]:#1 nt is the necessary for the formation of the loop
                MMEJ_F = True

            else:
                MMEJ_F = False
                print ("MMEJ_F not found")
            
            if MMEJ_F:
                print("MMEJ_F found!")
                
                if not self.extension:
                    print("OFF")
                    print(self.windowsize)
                    for MH_lengths in self.MH_lengths:
                        search_window = self.INDEL
                        found, length = MMEJ_found(MH, search_window)

    

                    MMEJ_motif_found.append(found)

                    MMEJ_motif_len.append(length)
                    

                    _d = {
                        'MMEJ': MMEJ_motif_found,
                        'MMEJ_motif_len': MMEJ_motif_len
                        }
        
                    self.del_out_dict = _d
                if self.extension:
                    print("ON")
                    window = self.ref_genome_down[:len(self.INDEL)]

                    for MH_length in self.MH_lengths:
                        MH = self.ref_genome_up[-MH_length:]

                    # estendo a sinistra finché il motivo continua a essere trovato
                        if MH in window:
                            for i in range(1, len(self.INDEL) + 1):
                                candidate = self.ref_genome_up[-(MH_length + i):]
                                if len(candidate) < MH_length + i:   # finita la sequenza
                                    break
                                if candidate in window:
                                    MH = candidate
                                else:
                                    break

                        print("MH esteso(MMEJ): " + MH)
                        found, length = MMEJ_found(MH, window)
                        MMEJ_motif_found.append(found)
                        MMEJ_motif_len.append(length)

                    self.del_out_dict = {
                        'MMEJ': MMEJ_motif_found,
                        'MMEJ_motif_len': MMEJ_motif_len,
                        }

    def pol_slip(self,xINDEL):

        """
        A function that looks for the pattern of a polymerase slippage, i.e. a submotif preceding the indel that is repeated tp form the indel

        """
        subrepeat=self.find_subrepeat(xINDEL)
        nrepeats_upstream=0
        if self.ref_genome_up[(len(self.ref_genome_up)-len(subrepeat[0])):]==subrepeat[0]: 
            nrepeats_upstream=1
            nrepeats_downstream=0
            isrepeated=True
            while isrepeated:
                 nrepeats_downstream=nrepeats_downstream+1
                 #upstream_repeat=self.DSB_up[(len(self.DSB_up)-len(subrepeat[0])*nrepeats_upstream):]
                 #NB: variable is called upstream_repeat because first implemented looking upstream. Currently downstream 
                 downstream_repeat=self.ref_genome_down[:len(subrepeat[0])*nrepeats_downstream]
                 if downstream_repeat!=(subrepeat[0]*nrepeats_downstream):
                     isrepeated=False
                     nrepeats_downstream=nrepeats_downstream-1 
        if nrepeats_upstream>0:
            #print(xINDEL+" - " +subrepeat[0] + "x" + str(nrepeats_upstream) + " - " + self.DSB_up+ " - " + self.DSB_down + " ref: " + self.ref_genome_up + " - " + self.ref_genome_down)
            pol_slip = True
            pol_slip_submotif = subrepeat[0]
            pol_slippage_repeatsIndel = subrepeat[1]
            pol_slippage_repeatsDownstream = nrepeats_downstream #NB:repeats upstream includes the indel, thus goes from 1 (only the indel) to many
            pol_slip_pos=get_motifs_pos(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, motif=subrepeat[0], windowsize=len(subrepeat[0])*nrepeats_downstream)
           #print(pol_slip_pos)
        else:
            pol_slip = False
            pol_slip_submotif = '-'
            pol_slippage_repeatsIndel = subrepeat[1]
            pol_slippage_repeatsDownstream = 0
            pol_slip_pos='-'

        _d = { 'pol_slip': pol_slip ,'pol_slip_submotif': pol_slip_submotif ,
                    'pol_slippage_repeatsIndel': pol_slippage_repeatsIndel, 
                    'pol_slippage_repeatsDownstream':pol_slippage_repeatsDownstream,
                    'pol_slip_pos':pol_slip_pos
                    } 
                
        self.pol_slip_dict = _d
    
    def sd_inverted_deletions(self):
        """
        sd_inverted_deletions will detect sd_inverted_deletions by looking for the following pattern:
        [MH1->DEL->P1] -> random seq (with length>0) -> ([MH2-P2]) MH2-P2=revc of MH1-P1
        This pattern is based on the one represented in mfabio's drawings
        Args: 
            self attributes of the EMMEJrealignment object
        Returns:
            ID = Inverted Deletion

            'SD_inverted_deletion': True or false 
            'SD_ID_repeat_pat_len': SD_ID_repeat_pat's length
            'SD_ID_motif_pos' : position of the pattern from the position of the DSB (pos)
                   
            No extension: MH2 is the revc of MH found downstream
                          P2 is the revc of MH found downstream
            Extension: MH1 is MH not extended
                       MH1_1 is MH during the extension
                       MH1_2 is MH with max extension
                       
                       P1 is P not extended
                       P1_1 is P during extension
                       P1_2 is P with max extension
        """
        #creating lists
        rep_pat_len, rep_pat_2, ID_motif_found, SD_ID_motif_pos_3, SD_ID_motif_freq_large, SD_ID_motif_freq_small = list(), list(), list(), list(), list(), list()               
        ###1st step is to bring MH1 and P1 from the ref seq, define them and create the revc
        for MH_lengths in self.MH_lengths:
            MH1 = self.ref_genome_up [ - MH_lengths :]  
            P1 = self.ref_genome_down [(self.indel_length):(self.indel_length+MH_lengths)]
            MH1_P1 = MH1 + P1
            MH1_P1revc = self.reverse_complement_converter(seq = MH1_P1) #make the revc)
            #print("MH1_P1: " + MH1_P1)
            #print("self.INDEL :" + self.INDEL)
            #print("MH1: "+ MH1)
            #print("P1: "+ P1)
            #print("MH1_P1revc: "+ MH1_P1revc)
            #print("ex: " + str(self.extension))
        
            #searchind MH1_P1revc in self.DSB_down
            #if MH1_P1revc in self.DSB_down[(len(self.INDEL) + len(P1) + 1):]:
            if MH1_P1revc in self.ref_genome_down[(len(self.INDEL) + len(P1) + 1):]:
                SD_inverted_deletion = True

            else:
                SD_inverted_deletion = False
                print("BR not found")
            
            if SD_inverted_deletion:
                print("BR found!!!")
            
                
                if not self.extension:
                    #print("OFF")
                    rep_pat = MH1_P1revc
                    rep_pat_pos = self.mutant_sequence[
                    (self.indel_position+len(P1)):].index(rep_pat)
                 
                    rep_pat_pos_2 = self.mutant_sequence[
                        (self.indel_position + len(self.INDEL) + len (P1)):].index(rep_pat)
                    #print("rep_pat_pos: "+str (rep_pat_pos)) #distance from P1motif to P2revc motif
                    #print("rep_pat_pos_2: "+str (rep_pat_pos_2))#distance from pos to P2revc motif
        
                    #define MH2 and P2, MH2 = revc of MH1, P2 = revc of P1
                    P2 = self.mutant_sequence[self.indel_position + len(P1) + rep_pat_pos : 
                                self.indel_position + len(P1) + rep_pat_pos + len(MH1)]
                
                    MH2 = self.mutant_sequence[self.indel_position + len(P1) + rep_pat_pos + len(MH1):
                        self.indel_position + len(P1) + rep_pat_pos + len(MH1) + len(P1)]
                
                    #print("MH2: "+MH2)
                    #print("P2: "+P2)
        
                   
                    # +1 and - 1 because 1 nt is necessary due to the fomration of the loop in snap_back mechanism
                    up_start= len(self.INDEL) + len(P1) + 1
                    down_start = - len(MH1) - 1
                    #print(up_start)
                    #print(down_start)

                    for MH_lengths in self.MH_lengths:
                        pos=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)

                    SD_ID_motif_pos_3.append(pos)
                    ID_motif_found.append(len(pos) > 0)
                    #rep_pat_2.append(rep_pat)
                    rep_pat_len.append(len(rep_pat))

                    # setting output dict
                    _d = {
                        'SD_inverted_deletion': ID_motif_found,
                        
                        'SD_ID_repeat_pat_len':rep_pat_len,
                        
                        'SD_ID_motif_pos' : SD_ID_motif_pos_3,
                        }
                    self.snap_out_dict = _d
                    
                if self.extension:
                    # ---------- PATTERN DI PARTENZA ----------
                    L = MH_lengths        # lunghezza MH iniziale
                    MH1_rc = self.reverse_complement_converter(seq=MH1)
                    P1_rc = self.reverse_complement_converter(seq=P1)
                    search_pattern = P1_rc + MH1_rc

                    after_P = self.ref_genome_down[len(self.INDEL)+len(P1):] # tutto ciò che sta dopo dopo indel e dopo P1
                    distance_Y = after_P.find(search_pattern, 1)  # parte da 1: salta Y = 0
                    if distance_Y == -1:
                        raise ValueError("pattern non trovato con Y >= 1")
                    
                    
                    Y_len = distance_Y# Y = after_P[:idx]
                    print("Y_len: "+ str(Y_len))
                    assert Y_len >= 1, "Y deve essere almeno 1 nucleotide"

                    # inizio di rc(MH) in after_P (fisso: non dipende dalle estensioni di MH)
                    base_mh = len(P1_rc) + distance_Y
                    print("base_mh: "+str(base_mh))

                    # ---------- ESTENSIONE DI MH ----------
                    # MH cresce a sinistra in ref_genome_up, rc(MH) cresce a destra in after_P.
                    # Non dipende da Y: dipende solo da quanto sta a monte e a valle.
                    MH1_2 = MH1
                    MH1_2_rc = MH1_rc

                    for i in range(1, len(self.ref_genome_up) - L + 1):
                        MH1_i = self.ref_genome_up[-(L + i):] # MH esteso a sinistra di i
                        MH1_i_rc = self.reverse_complement_converter(seq=MH1_i)

                        # rc(MH esteso) = rc(MH) + rc(basi aggiunte): deve combaciare dopo l'INDEL
                        if after_P[base_mh:].startswith(MH1_i_rc):
                            MH1_2 = MH1_i                             # ultima estensione valida
                            MH1_2_rc = MH1_i_rc
                        else:
                            break
                    # ---------- ESTENSIONE DI P ----------
                    # con P esteso di j: Y = after_P[j : idx - j], lunghezza idx - 2j >= 1
                    P1_2 = P1
                    Y_2 = after_P[:distance_Y]
                    #II_pattern = MH1_2 + P1 + Y_2 + P1_rc + self.INDEL + MH1_2_rc

                    j_max = (distance_Y - 1) // 2

                    for j in range(1, j_max + 1):
                        # P non può crescere oltre quello che c'è a valle
                        if len(self.INDEL)+len(P1) + j > len(self.ref_genome_down):
                            break

                        P1_j = self.ref_genome_down[len(self.INDEL):len(self.INDEL)+len(P1) + j]     # P esteso a destra di j
                        P1_j_rc = self.reverse_complement_converter(seq=P1_j)
                        candidato = P1_j_rc + MH1_2_rc      # usa già l'rc(MH) esteso

                        # rc(P esteso) inizia j basi prima di rc(P1) originale
                        if after_P[distance_Y - j:].startswith(candidato): # se in after_p dopo distance_Y - j non inizia con i P1revc allora breakka
                            P1_2 = P1_j                               # ultima estensione valida
                            
                        else:
                            break

                    j_final = len(P1_2) - len(P1)
                    i_final = len(MH1_2) - L
                    P1_2_rc = self.reverse_complement_converter(seq=P1_2)
                    Y_2 = after_P[j_final: distance_Y - j_final]
                    ID_pattern = P1_2_rc +  MH1_2_rc

                    print("MH1_2:", MH1_2, "(+%d)" % i_final)
                    print("P1_2:", P1_2, "(+%d)" % j_final)
                    print("Y finale:", Y_2, "len:", len(Y_2))
                    print("ID_pattern:", ID_pattern)
                   
                    rep_pat = ID_pattern
                      
                    up_start= len(self.INDEL) + len(P1_2) + 1
                    down_start = - len(MH1_2) 
                    #print("rep_pat: "+rep_pat)
                    #print("up_start: "+str(up_start))

                    for MH_lengths in self.MH_lengths:
                          pos=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)

                    SD_ID_motif_pos_3.append(pos)
                    ID_motif_found.append(len(pos) > 0)
                    #rep_pat_2.append(rep_pat)
                    rep_pat_len.append(len(rep_pat))
  
                    # setting output dict
                    _d = {
                        'SD_inverted_deletion': ID_motif_found,
                        'SD_ID_repeat_pat_len':rep_pat_len,
                        'SD_ID_motif_pos' : SD_ID_motif_pos_3,
                         }
                    self.snap_out_dict = _d
                              
        

    
    def sd_snap_back_MMEJ(self):
        """    
        sd_snap_back_MMEJ will detect SD-Snap back MMEJ insertions  by looking for the following pattern:
        [MH2->INS->P2] -> random seq (with length>0) -> invert_and_complement([MH2->INS->P2])
        This pattern is based on the one represented in fabio's drawings:
        
        Args: 
            self attributes of the EMMEJrealignment object
        Returns:
            II = Inverted Insertion
            'SD_inverted_insertion': True or false 
            'SD_II_repeat_pat_len': SD_II_repeat_pat's length
            'SD_II_motif_pos' : position of the pattern from the position of the DSB (pos)
                   
            No extension: MH2 is the revc of MH found downstream
                          P2 is the revc of MH found downstream
            Extension: MH1 is MH not extended
                       MH1_1 is MH during the extension
                       MH1_2 is MH with max extension
                       
                       P1 is P not extended
                       P1_1 is P during extension
                       P1_2 is P with max extension 

            
        """
        SD_II_last_dimer,rep_pat_len,rep_pat_2, II_motif_found,SD_II_motif_pos_3, SD_II_motif_freq_large, SD_II_motif_freq_small = list(), list(),list(),list(),list(), list(), list()       #creation of the lists               
        
        ###1st step is to bring MH1 and P1 from the ref seq 
        for MH_lengths in self.MH_lengths:
            MH1 = self.ref_genome_up [ - MH_lengths :]  
            P1 = self.ref_genome_down [ : + MH_lengths ]
            #defining the II_pattern 
            DI_pattern = MH1 + self.INDEL + P1
            II_pattern = self.reverse_complement_converter(seq = DI_pattern)
            #print("DI_pattern: " + DI_pattern)
            #print ("II-pattern_before elongation: " + II_pattern)
            #print("self.INDEL :" + self.INDEL)
            #print("MH1: "+ MH1)
            #print("P1: "+ P1)
            #print("ex: " + str(self.extension))
            #searching II_pattern in downstream in the genome
        
            if II_pattern in self.ref_genome_down[(len(P1) + 1):]:#1 nt is the necessary for the formation of the loop
                SD_inverted_insertion = True

            else:
                SD_inverted_insertion = False
                #print ("BZR not found")
            
            if SD_inverted_insertion:
                #print("BZR found!")
                 # Creating mmej_marked
                
                if not self.extension:
                    #in not elongation we define MH2 and P2 according to the distance of the pattern from the pos 
                    #print("OFF")
                    rep_pat = II_pattern
                    rep_pat_pos_2 = self.mutant_sequence[
                    (self.indel_position):].index(rep_pat)
                
                    rep_pat_pos = self.mutant_sequence[
                        (self.indel_position + len(self.INDEL) + len (P1)):].index(rep_pat)
                    #print("rep_pat_pos: "+str (rep_pat_pos)) #distance between P  e Prevc (Y seq)
                    #print("rep_pat_pos_2: "+str (rep_pat_pos_2))#distance between pos e Prevc

                    #define MH2 e P2
                    P2 = self.mutant_sequence[
                            self.indel_position + len(self.INDEL) + len(P1) + rep_pat_pos : 
                                self.indel_position + len(self.INDEL) + len(P1) + rep_pat_pos + len(MH1)]
                
                    MH2  = self.mutant_sequence[
                            self.indel_position + len(self.INDEL) + len(P1) + 
                            rep_pat_pos + len(MH1)+ len(self.INDEL) :
                            self.indel_position + len(self.INDEL) + len(P1) +
                            rep_pat_pos + len(MH1) + len(self.INDEL) + len(P1)]
                               
                    up_start= len(P2) + 1
                    down_start = - len(MH2) 
                    #print("rep_pat: "+rep_pat)
                    #print("up_start: "+str(up_start))
                    #print("down_start : "+str (down_start))
                    for MH_lengths in self.MH_lengths:
                        pos=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)

                    SD_II_motif_pos_3.append(pos)
                    
                    II_motif_found.append(len(pos) > 0)
                    
                    rep_pat_len.append(len(rep_pat) - len(self.INDEL))
                                          
                      # setting output dict
                    _d = {
                          'SD_inverted_insertion': II_motif_found,
                          'SD_II_repeat_pat_len':rep_pat_len,
                          'SD_II_motif_pos' : SD_II_motif_pos_3,
                           }
                    self.snap_out_dict = _d

                if self.extension:
                    # ---------- PATTERN DI PARTENZA ----------
                    L = MH_lengths        # lunghezza MH iniziale
                    MH1_rc = self.reverse_complement_converter(seq=MH1)
                    P1_rc = self.reverse_complement_converter(seq=P1)
                    self.INDEL_revc = self.reverse_complement_converter(seq=self.INDEL)
                    search_pattern = P1_rc + self.INDEL_revc + MH1_rc # rc(P1) + INDEL + rc(MH1)

                    after_P = self.ref_genome_down[len(P1):] # tutto ciò che sta dopo P1
                    distance_Y = after_P.find(search_pattern,1)# posizione di rc(P1) in after_P
                    if distance_Y == -1:
                        raise ValueError("pattern non trovato con Y >= 1")

                    Y_len = distance_Y# Y = after_P[:idx]
                    print("Y_len: "+ str(Y_len))
                    assert Y_len >= 1, "Y deve essere almeno 1 nucleotide"

                    # inizio di rc(MH) in after_P (fisso: non dipende dalle estensioni di MH)
                    base_mh = distance_Y + len(P1_rc) + len(self.INDEL)

                    # ---------- ESTENSIONE DI MH ----------
                    # MH cresce a sinistra in ref_genome_up, rc(MH) cresce a destra in after_P.
                    # Non dipende da Y: dipende solo da quanto sta a monte e a valle.
                    MH1_2 = MH1
                    MH1_2_rc = MH1_rc

                    for i in range(1, len(self.ref_genome_up) - L + 1):
                        MH1_i = self.ref_genome_up[-(L + i):] # MH esteso a sinistra di i
                        MH1_i_rc = self.reverse_complement_converter(seq=MH1_i)

                        # rc(MH esteso) = rc(MH) + rc(basi aggiunte): deve combaciare dopo l'INDEL
                        if after_P[base_mh:].startswith(MH1_i_rc):
                            MH1_2 = MH1_i                             # ultima estensione valida
                            MH1_2_rc = MH1_i_rc
                        else:
                            break

                    # ---------- ESTENSIONE DI P ----------
                    # con P esteso di j: Y = after_P[j : idx - j], lunghezza idx - 2j >= 1
                    P1_2 = P1
                    Y_2 = after_P[:distance_Y]
                    #II_pattern = MH1_2 + P1 + Y_2 + P1_rc + self.INDEL + MH1_2_rc

                    j_max = (distance_Y - 1) // 2

                    for j in range(1, j_max + 1):
                        # P non può crescere oltre quello che c'è a valle
                        if len(P1) + j > len(self.ref_genome_down):
                            break

                        P1_j = self.ref_genome_down[:len(P1) + j]     # P esteso a destra di j
                        P1_j_rc = self.reverse_complement_converter(seq=P1_j)
                        self.INDEL_rc = self.reverse_complement_converter(seq=self.INDEL)
                        candidato = P1_j_rc + self.INDEL_rc + MH1_2_rc      # usa già l'rc(MH) esteso

                        # rc(P esteso) inizia j basi prima di rc(P1) originale
                        if after_P[distance_Y - j:].startswith(candidato):
                            P1_2 = P1_j                               # ultima estensione valida
                        else:
                            break

                    j_final = len(P1_2) - len(P1)
                    i_final = len(MH1_2) - L
                    Y_2 = after_P[j_final: distance_Y - j_final]
                    P1_2_rc = self.reverse_complement_converter(seq=P1_2)
                    # pattern finale: solo rc(P esteso) + rc(INDEL) + rc(MH esteso)
                    II_pattern = P1_2_rc + self.INDEL_revc + MH1_2_rc


                    print("MH1_2:", MH1_2, "(+%d)" % i_final)
                    print("P1_2:", P1_2, "(+%d)" % j_final)
                    print("Y finale:", Y_2, "len:", len(Y_2))
                    print("II_pattern:", II_pattern)
                    rep_pat = II_pattern
                    print ("final pattern found: " + rep_pat)
                     
                    up_start= len(P1_2) + 1
                    down_start = - len(MH1_2) 
                    #print("rep_pat: "+rep_pat)
                    #print("up_start: "+str(up_start))
                    #print("down_start : "+str (down_start))
                    for MH_lengths in self.MH_lengths:
                      pos=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)

                    SD_II_motif_pos_3.append(pos)
                    
                    II_motif_found.append(len(pos) > 0)
                    
                    rep_pat_len.append(len(rep_pat) - len(self.INDEL))
                    _d = {
                       'SD_inverted_insertion': II_motif_found,
                       
                       'SD_II_repeat_pat_len': rep_pat_len,
                       
                       'SD_II_motif_pos' : SD_II_motif_pos_3,
                        }
                    self.snap_out_dict = _d



    def sd_direct_deletions(self):
        """
        
        sd_direct_deletion will detect sd direct deletions by looking for the following pattern:
        [MH1->DEL->P1] -> random seq (with length>0) -> ([MH2-P2]) 
        This pattern is based on the one represented in fabio's drawings
        Args: 
            self attributes of the EMMEJrealignment object
        Returns:
            DD = Inverted Deletion

            'SD_direct_deletion': True or false 
            'SD_DD_repeat_pat_len': SD_DD_repeat_pat's length
            'SD_DD_motif_pos' : position of the pattern from the position of the DSB (pos)
                   
            No extension: MH2 is the revc of MH found downstream
                          P2 is the revc of MH found downstream
            Extension: MH1 is MH not extended
                       MH1_1 is MH during the extension
                       MH1_2 is MH with max extension
                       
                       P1 is P not extended
                       P1_1 is P during extension
                       P1_2 is P with max extension
       
            
        """  
        rep_pat_len,rep_pat_2, DD_motif_found, SD_DD_motif_pos_3, SD_DD_motif_freq_large, SD_DD_motif_freq_small = list(), list(), list(), list(), list(), list()               
        ###1st step is to bring MH1 and P1 from the ref seq 
        for MH_lengths in self.MH_lengths:
            MH1 = self.ref_genome_up [ - MH_lengths :]  
            P1 = self.ref_genome_down [(self.indel_length):(self.indel_length+MH_lengths)]
            MH1_P1 = MH1 + P1
            #searching the pattern in self.DSB_down
        
            if MH1_P1 in self.DSB_down[(len(self.INDEL) + len(P1) + 0):]:#1 nt is the necessary for the formation of the loop
                SD_direct = True

            else:
                SD_direct = False
                #print ("BR not found")
            
            if SD_direct:
                #print("BR found!")
            
                if not self.extension:
                    #print("OFF")
                    rep_pat = MH1_P1
                    rep_pat_pos = self.mutant_sequence[
                    (self.indel_position+len(P1)):].index(rep_pat)
                
                
                
                    rep_pat_pos_2 = self.mutant_sequence[
                        (self.indel_position + len(self.INDEL) + len (P1)):].index(rep_pat)
                    #print("rep_pat_pos: "+str (rep_pat_pos)) #distance from P1 to MH2
                    #print("rep_pat_pos_2: "+str (rep_pat_pos_2))#distance from pos to MH2

                    #define MH2 e P2 according to distance from P1 to MH2
                    MH2 = self.mutant_sequence[self.indel_position + len(P1) + rep_pat_pos : 
                                self.indel_position + len(P1) + rep_pat_pos + len(MH1)]
                
                    P2 = self.mutant_sequence[self.indel_position + len(P1) + rep_pat_pos + len(MH1):
                        self.indel_position + len(P1) + rep_pat_pos + len(MH1) + len(P1)]
                
                    up_start = len(self.INDEL) + len (P1) 
                    down_start = - len(MH1)
                    #print(up_start)
                    #print(down_start)
                    #print(self.pos_on_chr)
                    #print(rep_pat)
                      
                    
                    for MH_lengths in self.MH_lengths:
                        pos=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)
                    
                    SD_DD_motif_pos_3.append(pos)
                   
                    DD_motif_found.append(len(pos) > 0)
                   
                    rep_pat_len.append(len(rep_pat))


                    _d = {
                        'SD_direct_deletion':DD_motif_found,
                        'SD_DD_repeat_pat_len': rep_pat_len,                        
                        'SD_DD_motif_pos' : SD_DD_motif_pos_3,
                        }
                    self.loop_out_dict = _d

                if self.extension:
                    print("ON")
                    rep_pat = MH1 + P1 # MH1  + P1, non esteso
                    
                    after_P = self.ref_genome_down[len(self.INDEL)+len(P1):]
                    Y_space = after_P.index(rep_pat)# distanza fra P1 e MH2
                    # valori di partenza: nessuna estensione
                    MH1_2 = MH1
                    MH1_2_P1 = rep_pat
                    tmp_pos = Y_space

                    for i in range(1, Y_space + 1):
                        # MH non può crescere oltre quello che c'è a monte
                        if MH_lengths + i > len(self.ref_genome_up):
                            break
                        MH1_1 = self.ref_genome_up[-(MH_lengths + i):]
                        MH1_1_P1 = MH1_1 + P1

                        # MH2 esteso inizia i basi prima nella reference
                        if after_P[Y_space - i:].startswith(MH1_1_P1): 
                            #lo cerca fino a quando after P dopo Y_space - i inizia con il pattern allungato
                            MH1_2 = MH1_1                     # ultima estensione valida
                            MH1_2_P1 = MH1_1_P1
                            tmp_pos = Y_space - i             # spazio residuo per estendere P
                        else:
                            break
                    # ---------- ESTENSIONE DI P ----------
                    # confinata in after_P[:tmp_pos] (prima dell'inizio di MH2 esteso)
                    # livello 1: fuori dal for di MH

                    # valori di partenza: nessuna estensione di P
                    P1_2 = P1
                    DD_pattern = MH1_2_P1    # MH1_2 + P1
                    for j in range(1, tmp_pos + 1):                   # livello 1: stesso del for di MH
                        # P non può crescere oltre quello che c'è a valle
                        if len(self.INDEL)+len(P1) + j > len(self.ref_genome_down):
                            break
                        P1_1 = self.ref_genome_down[len(self.INDEL):len(self.INDEL)+len(P1) + j]
                        MH1_2_P1_1 = MH1_2 + P1_1

                    # MH2 esteso inizia a tmp_pos basi dopo la fine di P1
                        if after_P[tmp_pos:].startswith(MH1_2_P1_1):
                            P1_2 = P1_1                               # ultima estensione valida
                            DD_pattern = MH1_2_P1_1
                        else:
                            break

                    #print("P1_2: " + P1_2)
                    #print("DD_pattern: " + DD_pattern)
                    
                    rep_pat = DD_pattern
                    #rep_pat_pos = self.mutant_sequence[
                    #     (self.indel_position+len(P1_2)):].index(rep_pat)
                    up_start = len(self.INDEL) + len (P1_2) 
                    down_start = - len(MH1_2)  
                    #print(up_start)
                    #print(down_start)
                    #print("pos on chrom: "+self.pos_on_chr)
              
                    for MH_lengths in self.MH_lengths:
                        pos=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)

                    SD_DD_motif_pos_3.append(pos)
                    DD_motif_found.append(len(pos) > 0)
                    #rep_pat_2.append(rep_pat)
                    rep_pat_len.append(len(rep_pat))
              
                    _d = {
                        'SD_direct_deletion': DD_motif_found,
                        'SD_DD_repeat_pat_len':rep_pat_len,
                        'SD_DD_motif_pos' : SD_DD_motif_pos_3
                           }
                    self.loop_out_dict = _d
            
        
    def sd_direct_insertions(self):
        """
        ####### CORRECTION ###########
        sd_direct_insertions will detect sd loop out MMEJ by looking for the following pattern:
        [MH1->INS->P1] -> random seq -> [MH2->INS->P2]
        This pattern is based on the one represented in the fabios's notebook.
        
        Procedure works by first expanding MH2, then P2. This does not guarantee finding all, nor longest.
        Not a very good approach.

        Args: 
            self attributes of the EMMEJrealignment object

        returns:
            SD_direct_insertion: True/False
            'SD_DI_repeat_pat_len': length of SD_DI_repeat_pat,
            'SD_DI_motif_pos' : positions of BZR motif 
                    
            
        """
        rep_pat_len,rep_pat_2, DI_motif_found ,SD_DI_motif_pos_3, SD_DI_motif_freq_large, SD_DI_motif_freq_small = list(), list(), list(), list(), list(), list()

        ###1st step is to bring MH1 and P1 from the ref seq 
        for MH_lengths in self.MH_lengths:
            MH1 = self.ref_genome_up [ - MH_lengths :]  
            #P1 = self.ref_genome_down [(self.indel_length):(self.indel_length+MH_lengths)]
            P1 = self.ref_genome_down [ : + MH_lengths ]
            #MH1_P1 = MH1 + P1
            DI_pattern = MH1 + self.INDEL + P1 
            #print("DI_pattern: " + DI_pattern)
            #print("MH1_P1: " + MH1_P1)
            #print("self.INDEL :" + self.INDEL)
            #print("MH1: "+ MH1)
            #print("P1: "+ P1)
            #print("ex: " + str(self.extension))
        
            if DI_pattern in self.ref_genome_down[(len(P1)) : ]:#1 nt is the necessary for the formation of the loop
                SD_direct_insertion = True

            else:
                SD_direct_insertion = False
                print ("BZR not found")
            
            if SD_direct_insertion:
                print("BZR found!")
        
                # Creating mmej_marked
                if not self.extension:
                    #MH1 = non allungato, MH1_1 in fase di allungamento, MH1_2 =  MH1 a fine allungamento 
                    print("OFF")
                    rep_pat = DI_pattern
                    rep_pat_pos_2 = self.mutant_sequence[
                    (self.indel_position):].index(rep_pat)
                
                    rep_pat_pos = self.mutant_sequence[
                        (self.indel_position + len(self.INDEL) + len (P1)):].index(rep_pat)
                    #print("rep_pat_pos: "+str (rep_pat_pos)) #distance between P1 e MH2 (Y seq)
                    #print("rep_pat_pos_2: "+str (rep_pat_pos_2))#distance between pos e MH2

                    #definisci MH2 e P2
                    MH2 = self.mutant_sequence[
                            self.indel_position + len(self.INDEL) + len(P1) + rep_pat_pos : 
                                self.indel_position + len(self.INDEL) + len(P1) + rep_pat_pos + len(MH1)]
                
                    P2 = self.mutant_sequence[
                            self.indel_position + len(self.INDEL) + len(P1) + 
                            rep_pat_pos + len(MH1)+ len(self.INDEL) :
                            self.indel_position + len(self.INDEL) + len(P1) +
                            rep_pat_pos + len(MH1) + len(self.INDEL) + len(P1)]
                
                    up_start= len(P2)
                    down_start = - len(MH2) 
                    #print("rep_pat: "+rep_pat)
                    #print("up_start: "+str(up_start))
                    #print("down_start : "+str (down_start))
                    
                    for MH_lengths in self.MH_lengths:
                        pos=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)

                    SD_DI_motif_pos_3.append(pos)
                    DI_motif_found.append(len(pos) > 0)
                    #rep_pat_2.append(rep_pat)
                    rep_pat_len.append(len(rep_pat)-len(self.INDEL))

                          
                    #setting output dict
                    _d = {
                        'SD_direct_insertion': DI_motif_found, 
                        'SD_DI_repeat_pat_len':rep_pat_len,
                        'SD_DI_motif_pos': SD_DI_motif_pos_3
                        }
                    self.loop_out_dict = _d 
                 
                if self.extension:
                    print("ON")
                    rep_pat = DI_pattern                      # MH1 + INDEL + P1, non esteso
                    after_P = self.ref_genome_down[len(P1):]
                    Y_space = after_P.index(rep_pat)          # distanza fra P1 e MH2
                    #print(Y_space)

                    # valori di partenza: nessuna estensione
                    MH1_2 = MH1
                    MH1_2_indel_P1 = rep_pat
                    tmp_pos = Y_space

                    for i in range(1, Y_space + 1):
                        # MH non può crescere oltre quello che c'è a monte
                        if MH_lengths + i > len(self.ref_genome_up):
                            break
                        MH1_1 = self.ref_genome_up[-(MH_lengths + i):]
                        MH1_1_indel_P1 = MH1_1 + self.INDEL + P1

                        # MH2 esteso inizia i basi prima nella reference
                        if after_P[Y_space - i:].startswith(MH1_1_indel_P1):
                            MH1_2 = MH1_1                     # ultima estensione valida
                            MH1_2_indel_P1 = MH1_1_indel_P1
                            tmp_pos = Y_space - i             # spazio residuo per estendere P
                        else:
                            break

                    # ---------- ESTENSIONE DI P ----------
                    # confinata in after_P[:tmp_pos] (prima dell'inizio di MH2 esteso)
                    #print("tmp_pos: " + str(tmp_pos))                 # livello 1: fuori dal for di MH

                    # valori di partenza: nessuna estensione di P
                    P1_2 = P1
                    DI_pattern = MH1_2_indel_P1                       # MH1_2 + INDEL + P1

                    for j in range(1, tmp_pos + 1):                   # livello 1: stesso del for di MH
                        # P non può crescere oltre quello che c'è a valle
                        if len(P1) + j > len(self.ref_genome_down):
                            break
                        P1_1 = self.ref_genome_down[:len(P1) + j]
                        MH1_2_indel_P1_1 = MH1_2 + self.INDEL + P1_1

                    # MH2 esteso inizia a tmp_pos basi dopo la fine di P1
                        if after_P[tmp_pos:].startswith(MH1_2_indel_P1_1):
                            P1_2 = P1_1                               # ultima estensione valida
                            DI_pattern = MH1_2_indel_P1_1
                        else:
                            break

                    #print("P1_2: " + P1_2)
                    #print("DI_pattern: " + DI_pattern)
                    """
                      # elongation of P2, confined in extention_space[:*starting position of MH2_2*]
                      print("tmp_pos: "+str(tmp_pos))
                      
                      for i in range(1, tmp_pos):
                          P1_1 = self.ref_genome_down[:( + MH_lengths + i)]
                          MH1_2_indel_P1_1 = MH1_2 + self.INDEL + P1_1

                          if MH1_2_indel_P1_1 in self.ref_genome_down [ len(P1):]:
                              continue
                          else:
                              P1_2 = P1_1[:-1]
                              #print("P1_2: "+P1_2)
                              DI_pattern = MH1_2 + self.INDEL + P1_2
                              #print("DI_pattern : "+ DI_pattern)
                              break
                      else:
                          DI_pattern=MH1_2_indel_P1 #se P2 non viene allungato allora il BR rimane con R non allungato
                          P1_2 = P1_1
                    """
                    
                    rep_pat = DI_pattern
                    #print(self.pos_on_chr)
                    #print(MH1_2)
                    rep_pat_pos_3 = self.mutant_sequence[
                         (self.indel_position+self.indel_length+len(P1_2)):].index(rep_pat)
                    #print("rep_pat_pos_3: "+str (rep_pat_pos_3)) #distanza fra P1 e MH2
                      
                    up_start= len(P1_2)
                    down_start = - len(MH1_2) 
                    #print("rep_pat: "+rep_pat)
                    #print("up_start: "+str(up_start))
                    #print("down_start : "+str (down_start))
                    for MH_lengths in self.MH_lengths:
                         pos=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)

                    SD_DI_motif_pos_3.append(pos)
                    DI_motif_found.append(len(pos) > 0)
                    #rep_pat_2.append(rep_pat)
                    rep_pat_len.append(len(rep_pat) - len(self.INDEL))

                    #setting output dict
                    _d = {
                        'SD_direct_insertion': DI_motif_found, 
                        'SD_DI_repeat_pat_len':rep_pat_len,
                        'SD_DI_motif_pos': SD_DI_motif_pos_3,
                        }
                    self.loop_out_dict = _d 

    def sd_direct_substitution(self):
        #print("self.INDEL :"+self.INDEL)
        #print("self.ANC:"+self.ANC)
        #print("self.DER:"+self.DER)
        """
        """

        rep_pat_len,rep_pat_2, DS_motif_found,SD_DS_motif_pos_3, SD_DS_motif_freq_large, SD_DS_motif_freq_small = list(), list(), list(), list(), list(), list()               
        ###1st step is to bring MH1 and P1 from the ref seq 
        for MH_lengths in self.MH_lengths:
            MH1 = self.ref_genome_up [ - MH_lengths :]  
            P1 = self.ref_genome_down [ len(self.ANC) : len(self.ANC) + MH_lengths ]
            DS_pattern = MH1 + self.INDEL + P1 
            #print("DS_pattern: " + DS_pattern)
            #print("MH1_P1: " + MH1_P1)
            #print("self.INDEL :" + self.INDEL)
            #print("MH1: "+ MH1)
            #print("P1: "+ P1)
            #print("ex: " + str(self.extension))
            
            if DS_pattern in self.ref_genome_down [len(self.ANC) + len(P1):]:
                SD_direct_substitution = True
            else:
                SD_direct_substitution = False 
                #print("BZR direct not found")

            if SD_direct_substitution:
                #print("BZR direct found")
                
                if not self.extension:
                    #MH1 = non allungato, MH1_1 in fase di allungamento, MH1_2 =  MH1 a fine allungamento 
                    print("OFF")
                    rep_pat = DS_pattern
                    rep_pat_pos_2 = self.mutant_sequence[
                    (self.indel_position):].index(rep_pat)
                    #print(self.indel_position)
                    #print(self.INDEL)
                    #print(P1)
                    #print("distance between pos e MH2: "+str(rep_pat_pos_2))
                    #print(self.pos_on_chr)
                    rep_pat_pos = self.mutant_sequence[
                        (self.indel_position + len(self.INDEL) + len (P1)):].index(rep_pat)
                    #print("rep_pat_pos: "+str (rep_pat_pos)) #distance between P1 e MH2 (Y seq)
                    #print("rep_pat_pos_2: "+str (rep_pat_pos_2))#distance between pos e MH2

                    #definisci MH2 e P2
                    MH2 = self.mutant_sequence[
                            self.indel_position + len(self.INDEL) + len(P1) + rep_pat_pos : 
                                self.indel_position + len(self.INDEL) + len(P1) + rep_pat_pos + len(MH1)]
                
                    P2 = self.mutant_sequence[
                            self.indel_position + len(self.INDEL) + len(P1) + 
                            rep_pat_pos + len(MH1)+ len(self.INDEL) :
                            self.indel_position + len(self.INDEL) + len(P1) +
                            rep_pat_pos + len(MH1) + len(self.INDEL) + len(P1)]
                
                    #print("MH2: "+MH2)
                    #print("P2: "+P2)
                    
                    up_start = len(self.ANC) + len (P2)  
                    down_start = - len(MH2) 
                    #print(up_start)
                    #print(down_start)
                    
                    for MH_lengths in self.MH_lengths:
                        pos=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)

                    SD_DS_motif_pos_3.append(pos)
                    DS_motif_found.append(len(pos) > 0)
                    #rep_pat_2.append(rep_pat)
                    rep_pat_len.append(len(rep_pat) - len(self.INDEL))

                          
                    #setting output dict
                    _d = {
                        'SD_direct_substitution': DS_motif_found, 
                        
                        'SD_DS_repeat_pat_len': rep_pat_len,
                        
                        'SD_DS_motif_pos': SD_DS_motif_pos_3
                         }
                    self.loop_out_dict = _d 
                
                if self.extension:
                    print("ON")
                    rep_pat = DS_pattern                      # MH1 + INDEL + P1, non esteso
                    after_P = self.ref_genome_down[len(self.ANC)+ len(P1):]
                    Y_space = after_P.index(rep_pat)          # distanza fra P1 e MH2
                    print(Y_space)
                    print(rep_pat)
                    print(self.INDEL)

                    # valori di partenza: nessuna estensione
                    MH1_2 = MH1
                    MH1_2_indel_P1 = rep_pat
                    tmp_pos = Y_space

                    for i in range(1, Y_space + 1):
                        # MH non può crescere oltre quello che c'è a monte
                        if MH_lengths + i > len(self.ref_genome_up):
                            break
                        MH1_1 = self.ref_genome_up[-(MH_lengths + i):]
                        MH1_1_indel_P1 = MH1_1 + self.INDEL + P1

                        # MH2 esteso inizia i basi prima nella reference
                        if after_P[Y_space - i:].startswith(MH1_1_indel_P1):
                            MH1_2 = MH1_1                     # ultima estensione valida
                            MH1_2_indel_P1 = MH1_1_indel_P1
                            tmp_pos = Y_space - i             # spazio residuo per estendere P
                        else:
                            break
                    # ---------- ESTENSIONE DI P ----------
                    # confinata in after_P[:tmp_pos] (prima dell'inizio di MH2 esteso)
                    #print("tmp_pos: " + str(tmp_pos))                 # livello 1: fuori dal for di MH

                    # valori di partenza: nessuna estensione di P
                    P1_2 = P1
                    DS_pattern = MH1_2_indel_P1                       # MH1_2 + INDEL + P1

                    for j in range(1, tmp_pos + 1):                   # livello 1: stesso del for di MH
                        # P non può crescere oltre quello che c'è a valle
                        if len(P1) + j > len(self.ref_genome_down):
                            break
                        P1_1 = self.ref_genome_down[len(self.ANC):len(self.ANC)+len(P1) + j]
                        MH1_2_indel_P1_1 = MH1_2 + self.INDEL + P1_1

                    # MH2 esteso inizia a tmp_pos basi dopo la fine di P1
                        if after_P[tmp_pos:].startswith(MH1_2_indel_P1_1):
                            P1_2 = P1_1                               # ultima estensione valida
                            DS_pattern = MH1_2_indel_P1_1
                        else:
                            break
                    rep_pat = DS_pattern
                    print("rep_pat: "+rep_pat)
    
                    
                    up_start = len(self.ANC) + len (P1_2)  
                    down_start = - len(MH1_2) 
                    #print(up_start)
                    #print(down_start)
                                        
                    for MH_lengths in self.MH_lengths:
                        pos=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)

                    SD_DS_motif_pos_3.append(pos)
                    DS_motif_found.append(len(pos) > 0)
                    #rep_pat_2.append(rep_pat)
                    rep_pat_len.append(len(rep_pat) - len(self.INDEL))
                    
                    #setting output dict
                    _d = {
                        'SD_direct_substitution': DS_motif_found, 
                    
                        'SD_DS_repeat_pat_len': rep_pat_len,
                        
                        'SD_DS_motif_pos': SD_DS_motif_pos_3,
                         }
                    self.loop_out_dict = _d 

    def sd_inverted_substitution(self):
        #print("self.ANC:"+self.ANC)
        #print("self.DER:"+self.DER)
        """
        """
        
        rep_pat_len, rep_pat_2, IS_motif_found, SD_IS_motif_pos_3, SD_IS_motif_freq_large, SD_IS_motif_freq_small = list(),list(),list(),list(), list(), list()               
        ###1st step is to bring MH1 and P1 from the ref seq 
        for MH_lengths in self.MH_lengths:
            MH1 = self.ref_genome_up [ - MH_lengths :]  
            P1 = self.ref_genome_down [ len(self.ANC) : len(self.ANC) + MH_lengths ]
            #create the IS_pattern and search it
            IS_pattern = self.reverse_complement_converter(seq = MH1 + self.INDEL + P1)
            #print("IS_pattern: " + IS_pattern)
            #print("self.INDEL :" + self.INDEL)
            #print("MH1: "+ MH1)
            #print("P1: "+ P1)
            #print("ex: " + str(self.extension))
            #search the pattern after 1° P
            if IS_pattern in self.ref_genome_down [len(self.ANC) + len(P1) +1:]:
                SD_inverted_substitution = True
            else:
                SD_inverted_substitution = False 
                print("RZBrevc not found")

            if SD_inverted_substitution:
                print("RZBrevc found")
               
                if not self.extension:
                    #MH1/P1 = not extended, MH2/P2 = MH and P revc not extended 
                    print("OFF")
                    rep_pat = IS_pattern
                    #rep_pat_pos_2 is the distance of the pattern from DSB in the mutant context
                    rep_pat_pos_2 = self.mutant_sequence[
                    (self.indel_position):].index(rep_pat)
                    #rep_pat_pos is the distance of the pattern from 1° P in the mutant context
                    rep_pat_pos = self.mutant_sequence[
                        (self.indel_position + len(self.INDEL) + len (P1)):].index(rep_pat)
                    
                    #define MH2 e P2
                    Prevc = self.mutant_sequence[
                            self.indel_position + len(self.INDEL) + len(P1) + rep_pat_pos : 
                                self.indel_position + len(self.INDEL) + len(P1) + rep_pat_pos + len(MH1)]
                
                    MHrevc = self.mutant_sequence[
                            self.indel_position + len(self.INDEL) + len(P1) + 
                            rep_pat_pos + len(MH1)+ len(self.INDEL) :
                            self.indel_position + len(self.INDEL) + len(P1) +
                            rep_pat_pos + len(MH1) + len(self.INDEL) + len(P1)]
                
                    #print("MHrevc: " + MHrevc)
                    #print("Prevc: " + Prevc)

                    up_start = len(self.ANC) + len (Prevc)  
                    down_start = - len(MHrevc) 
                    #print(up_start)
                    #print(down_start)
                                       
                    for MH_lengths in self.MH_lengths:
                        pos=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)

                    SD_IS_motif_pos_3.append(pos)
                    IS_motif_found.append(len(pos) > 0)
                    #rep_pat_2.append(rep_pat)
                    rep_pat_len.append(len(rep_pat) - len(self.INDEL))
                    #setting output dict
                    _d = {
                        'SD_inverted_substitution': IS_motif_found, 
                       
                        'SD_IS_repeat_pat_len': rep_pat_len,
                       
                        'SD_IS_motif_pos': SD_IS_motif_pos_3,
                          }
                    self.snap_out_dict = _d 
                    """
                if self.extension:
                    print("ON")
                    rep_pat = IS_pattern
                    rep_pat_pos_2 = self.mutant_sequence[
                            (self.indel_position):].index(rep_pat)#distance from pos to Prevc
                    max_elongations_p1 = rep_pat_pos_2 // 2
                    MH1_1 = MH1
                    P1_1 = P1
                      

                    for i in range (1, rep_pat_pos_2):
                      ###elongation of MH, creating the pattern and search it downstream 
                         MH1_1 = self.ref_genome_up [(- MH_lengths -i):] 
                         #MH1_1_indel_P1 = MH1_1 + self.INDEL + P1
                         P1revc_indel_MH1_1revc = self.reverse_complement_converter(seq = MH1_1 + self.INDEL + P1) 
                         print("self.INDEL :" + self.INDEL)
                         print("MH1_1: "+ MH1_1)
                         print("P1: "+ P1)
                         if P1revc_indel_MH1_1revc in self.ref_genome_down [ len(self.ANC) + len(P1):]:
                             tmp_pos = self.ref_genome_down[ len(self.ANC) + (len(P1)):].index(P1revc_indel_MH1_1revc)
                             #tmp_pos is the distance of the BZR after elongation 
                             #print("tmp_pos: "+str(tmp_pos))

                         else:
                             MH1_2 = self.ref_genome_up[(- MH_lengths -i +1):]
                             print("elongations finished")
                             print("MH1_2: "+str(MH1_2))
                             P1revc_indel_MH1_1revc = self.reverse_complement_converter(seq = MH1_2 + self.INDEL + P1) 
                             print("P1revc_indel_MH1_1revc : " + P1revc_indel_MH1_1revc)
                             tmp_pos_2 = self.ref_genome_down[
                                     len(self.ANC) + len(P1):].index(P1revc_indel_MH1_1revc)
                             max_elongations_p2 = tmp_pos_2 // 2
                             #print("tmp_pos_2: " +str(tmp_pos_2))
                             break
                    else:
                        P1revc_indel_MH1_1revc = rep_pat
                        MH1_2 = self.reverse_complement_converter(seq = MH1_1)
                        max_elongations_p2 = max_elongations_p1


                    for i in range(1, max_elongations_p2):
                        P1_1 = self.ref_genome_down[
                                len(self.ANC) : len(self.ANC) + MH_lengths+i]
                        print("P1_1:"+P1_1)
                        #MH1_2_indel_P1_1 = MH1_2 + self.INDEL + P1_1
                        #print(MH1_2_indel_P1_1)
                        P1_1revc_indel_MH1_2revc = self.reverse_complement_converter(seq = MH1_2 + self.INDEL + P1_1)

                        if P1_1revc_indel_MH1_2revc in self.ref_genome_down [ len(self.ANC) + len(P1):]:
                            continue
                        else:
                            P1_2 = P1_1[:-1]
                            print("P1_2: "+P1_2)
                            IS_pattern = self.reverse_complement_converter(seq = MH1_2 + self.INDEL + P1_2)
                            #DS_pattern = MH1_2 + self.INDEL + P1_2
                            print("P1_2revc_indel_MH1_2revc : "+ IS_pattern)
                            break
                    else:
                        
                        IS_pattern=self.reverse_complement_converter(seq = MH1_2 + self.INDEL + P1)
                        #if P2 doesn't get extended RZB remains with R not extendend P2 non viene allungato allora il BR rimane con R non allungato
                        P1_2 = self.reverse_complement_converter(seq = P1_1)
                    """
                if self.extension:
                # ---------- PATTERN DI PARTENZA ----------
                    L = MH_lengths        # lunghezza MH iniziale
                    MH1_rc = self.reverse_complement_converter(seq=MH1)
                    P1_rc = self.reverse_complement_converter(seq=P1)
                    self.INDEL_revc = self.reverse_complement_converter(seq=self.INDEL)

                    after_P = self.ref_genome_down[len(self.ANC)+ len(P1):] # tutto ciò che sta dopo P1
                    search_pattern = P1_rc + self.INDEL_revc + MH1_rc   # rc(P1) + INDEL + rc(MH1)
                    print("search_pattern: "+search_pattern)
                    distance_Y = after_P.find(search_pattern, 1)   # parte da 1: salta Y = 0
                    if distance_Y == -1:
                        raise ValueError("pattern non trovato con Y >= 1")      # posizione di rc(P1) in after_P
                    Y_len = distance_Y# Y = after_P[:idx]
                    print("Y_len: "+ str(Y_len))
                    assert Y_len >= 1, "Y deve essere almeno 1 nucleotide"

                    # inizio di rc(MH) in after_P (fisso: non dipende dalle estensioni di MH)
                    base_mh = distance_Y + len(P1_rc) + len(self.INDEL)
                    # ---------- ESTENSIONE DI MH ----------
                    # MH cresce a sinistra in ref_genome_up, rc(MH) cresce a destra in after_P.
                    # Non dipende da Y: dipende solo da quanto sta a monte e a valle.
                    MH1_2 = MH1
                    MH1_2_rc = MH1_rc

                    for i in range(1, len(self.ref_genome_up) - L + 1):
                        MH1_i = self.ref_genome_up[-(L + i):] # MH esteso a sinistra di i
                        MH1_i_rc = self.reverse_complement_converter(seq=MH1_i)

                        # rc(MH esteso) = rc(MH) + rc(basi aggiunte): deve combaciare dopo l'INDEL
                        if after_P[base_mh:].startswith(MH1_i_rc):
                            MH1_2 = MH1_i                             # ultima estensione valida
                            MH1_2_rc = MH1_i_rc
                        else:
                            break
                    # ---------- ESTENSIONE DI P ----------
                    # con P esteso di j: Y = after_P[j : idx - j], lunghezza idx - 2j >= 1
                    P1_2 = P1
                    Y_2 = after_P[:distance_Y]
                    #II_pattern = MH1_2 + P1 + Y_2 + P1_rc + self.INDEL + MH1_2_rc

                    j_max = (distance_Y - 1) // 2

                    for j in range(1, j_max + 1):
                        # P non può crescere oltre quello che c'è a valle
                        if len(P1) + j > len(self.ref_genome_down):
                            break

                        P1_j = self.ref_genome_down[len(self.ANC):len(self.ANC)+len(P1) + j]     # P esteso a destra di j
                        P1_j_rc = self.reverse_complement_converter(seq=P1_j)
                        self.INDEL_rc = self.reverse_complement_converter(seq=self.INDEL)
                        candidato = P1_j_rc + self.INDEL_rc + MH1_2_rc      # usa già l'rc(MH) esteso

                        # rc(P esteso) inizia j basi prima di rc(P1) originale
                        if after_P[distance_Y - j:].startswith(candidato): # se in after_p dopo distance_Y - j non inizia con i P1revc allora breakka
                            P1_2 = P1_j                               # ultima estensione valida
                            
                        else:
                            break

                    j_final = len(P1_2) - len(P1)
                    i_final = len(MH1_2) - L
                    Y_2 = after_P[j_final: distance_Y - j_final]
                    P1_2_rc = self.reverse_complement_converter(seq = P1_2)
                    IS_pattern = P1_2_rc + self.INDEL_revc + MH1_2_rc

                    print("MH1_2:", MH1_2, "(+%d)" % i_final)
                    print("P1_2:", P1_2, "(+%d)" % j_final)
                    print("Y finale:", Y_2, "len:", len(Y_2))
                    print("IS_pattern:", IS_pattern)
                    rep_pat = IS_pattern
                    
                    up_start = len(self.ANC) + len (P1_2_rc)
                    down_start = - len(MH1_2_rc) 
                    #print(up_start)
                    #print(down_start)
                                        
                    for MH_lengths in self.MH_lengths:
                        pos=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)

                    SD_IS_motif_pos_3.append(pos)
                    IS_motif_found.append(len(pos) > 0)
                    #rep_pat_2.append(rep_pat)
                    rep_pat_len.append(len(rep_pat))
       
                    #setting output dict
                    _d = {
                        'SD_inverted_substitution': IS_motif_found, 
                        
                        'SD_IS_repeat_pat_len': rep_pat_len,
                        
                        'SD_IS_motif_pos': SD_IS_motif_pos_3,
                        }
                    self.snap_out_dict = _d 


    def export_data(self):
        """
        Encapsulate all the relevant values into a dataframe object with standard
        format for all MMEJ pathways.
        Args:
            self emMEJrealignment inputs.
        Returns:
            ex_df (pd.DataFrame): dataframe object with all the relevant values
        """
        _colnames = [
            # deletions
            'MMEJ', 'MMEJ_motif_len',
                       
            # SD_inverted_insertion
            'SD_inverted_insertion',
            'SD_II_repeat_pat_len', 
            'SD_II_motif_pos',
            
            #SD-direct_deletion
            'SD_direct_deletion', 
            'SD_DD_repeat_pat_len',
            'SD_DD_motif_pos'

            # SD_direct_insertion
            'SD_direct_insertion',
            'SD_DI_repeat_pat_len',
            'SD_DI_motif_pos',
                        
            # SD_inverted_deletion
            'SD_inverted_deletion',
            'SD_ID_repeat_pat_len',
            'SD_ID_motif_pos',

            #SD_direct_substitution
            'SD_direct_substitution',
            
            'SD_DS_repeat_pat_len',
            
            'SD_DS_motif_pos',
            
            #SD_inverted_substitution
             'SD_inverted_substitution',
            
             'SD_IS_repeat_pat_len',
             
             'SD_IS_motif_pos',
                                      
            # polymerase slippage
            #'pol_slip' , 'pol_slip_submotif' , 'pol_slippage_repeatsIndel',
            #'pol_slippage_repeatsDownstream', 'pos_slip_pos'
            ]

        _Dtypes={
            'MMEJ':str,'MMEJ_motif_len':str,
           
            #SD_inverted_insertion
            'SD_inverted_insertion':str,
            'SD_II_repeat_pat_len':str,
            'SD_II_motif_pos':str,
            
            
            #SD-direct_deletion
            'SD_direct_deletion':str,
            'SD_DD_repeat_pat_len':str,
            'SD_DD_motif_pos':str,
            
            # SD_direct_insertion
            'SD_direct_insertion':str,
            'SD_DI_repeat_pat_len':str,
            'SD_DI_motif_pos':str, 

                      
            # SD_inverted_deletion
            'SD_inverted_deletion':str,
            'SD_ID_repeat_pat_len':str,
            'SD_ID_motif_pos':str, 
            
            #SD_direct_substitution
            'SD_direct_substitution':str,
         
            'SD_DS_repeat_pat_len':str,
            
            'SD_DS_motif_pos':str,
            
            #SD_inverted_substitution
            'SD_inverted_substitution':str,
            
            'SD_IS_repeat_pat_len':str,
            
            'SD_IS_motif_pos':str,
              }
            # polymerase slippage
            #'pol_slip':str, 'pol_slip_submotif':str, 'pol_slippage_times':str,
            #'pol_slippage_last_dimer':str, 'pol_slip_motif_len':str}
        #,dtype=_Dtypes)
        if self.include_context:
            _colnames = _colnames + ['ref_genome_context', 'mutant_sequence']
        ex_df = pd.DataFrame(columns=_colnames )   
        ex_df.loc[0,:] = np.nan
        # merging dicts from all paths
        d = {**self.del_out_dict, 
            **self.snap_out_dict, **self.loop_out_dict, **self.pol_slip_dict}
        if self.include_context: d = {**d, **self.context}
        """ 
        # setting the right values to match corresponding columns
        for k, v in d.items():
            if isinstance(v, (list, tuple)) and len(v) == 1:
                v = v[0]
            ex_df.at[0, k] = v
        return ex_df
        """
        def _flatten(x):
            #Appiattisce liste/tuple/set/array annidati in una lista piatta
            if isinstance(x, (list, tuple, set, np.ndarray, pd.Series)):
                out = []
                for i in x:
                    out.extend(_flatten(i))
                return out
            return [x]

            # setting the right values to match corresponding columns
        for k, v in d.items():
            if k == "SD_DD_motif_pos" and isinstance(v, (list, tuple)):
                # un gruppo per ogni mhl, separati da ";"  ->  "12,15,30;8,12"
                v = ";".join(",".join(map(str, _flatten(g))) for g in v)
            else:
                if isinstance(v, (list, tuple)) and len(v) == 1:
                    v = v[0]
                if isinstance(v, (list, tuple, set, np.ndarray, pd.Series)):
                    v = ",".join(map(str, _flatten(v)))

            ex_df.at[0, k] = v
        return ex_df
