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
            print("ins")
            self.indel_length = len(DER) - 1
            self.INDEL =  DER[len(ANC):] #fabri: is this for substitutions? Then why before length-1?
            self.ANC = ANC
            self.windowsize = windowsize
        if (self.indel_type == 'DEL'): 
            print("del")
            self.indel_length = len(ANC) - 1
            self.ANC = ANC[len(DER):]
            self.INDEL = ANC[len(DER):]
            self.windowsize = windowsize + self.indel_length
            self.indel_position = self.windowsize
        if (self.indel_type == 'SUB'):
            print("sub")
            if len(ANC) > len (DER):
                self.indel_length = len(ANC) - len (DER)
                self.INDEL = DER
                self.ANC = ANC
                self.DER = DER
                #self.INDEL = DER[1:] # self.INDEL = self.DER 
                #self.ANC = ANC[1:] #delete the anchor
                #self.DER = DER[1:]
                
            else: 
                self.indel_length = len(DER) - len(ANC)
                self.INDEL = DER
                self.ANC = ANC
                self.DER = DER
                #self.INDEL = DER[1:] # self.INDEL = self.DER
                #self.ANC = ANC[1:] #delete the anchor
                #self.DER = DER[1:]
                
            self.windowsize = windowsize

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
            if len(ANC) > len (DER):
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
        """
        Finds the longest subsequence (a motif) that exist both before and after
        the DSB, and ends at the last base of the deletion sequence
        """
        
        # set default values for mmej_cand
        del_mmej = False
        del_mmejl = list()
        del_mmej_cand, del_mmej_cand_len ,del_mmej_marked= list(), list(), list() #np.nan
        del_last_dimer, mmej_cand = list() , list()
        distance_to_indel_position, del_mmej_marked_on_ref = list(), list()
        del_mmej_motif_pos, del_mmej_freq_large, del_mmej_freq_small = list(), list(), list()
        # 1st condition:
        self.first_MH_finder() # finding all possible MH motifs
        self.first_MH_indicator()

        # 2nd condition
        # reverse integrating over all motifs and look whether it exists before the break. If not breaks and longest is considered as the potential MH
        mmej_cand_position = self.indel_length + 1 #fabri: why is this within the loop? should be above for efficiency
        imhl=-1
        for mh_motif in self.mmej_motifs[::-1]:
              if mh_motif == self.ref_genome_up[-len(mh_motif):]:
                    if ( 0 in self.MH_lengths and len(mh_motif)>max(self.MH_lengths) and imhl==-1) or (len(mh_motif) in self.MH_lengths):
                          imhl=imhl+1
                          del_mmej = True
                          del_mmejl.append(True)
                          mmej_cand.append(mh_motif)
                          del_mmej_cand.append(mh_motif)
                    #also, ugly and misleading name for variable. This is better as simply indel_length+1, or indel_length_plus1
                          distance_to_indel_position.append(self.indel_length - len(mmej_cand[imhl]))
                          del_mmej_marked.append(self.del_mutant_mmej_marked(mmej_cand[imhl]))
                          del_mmej_marked_on_ref.append(self.del_mmej_marked_on_reference(mmej_cand[imhl],mmej_cand_position))
                          #define last dimer, because used later for 2nd order Markov Chain                
                          if len(mmej_cand[imhl]) >= 2: del_last_dimer.append(mmej_cand[imhl][-2:])
                          else: del_last_dimer.append(self.ref_genome_up[-2] + mmej_cand[imhl])
                          del_last_dimer[imhl]= del_last_dimer[imhl].upper()
                          del_mmej_cand_len.append(int(len(del_mmej_cand[imhl])))
                          temp=get_motifs_freqs(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, large_window=1000,
                               small_window=self.windowsize,motif=mmej_cand[imhl],indel_type='DEL')
                          del_mmej_motif_pos.append(temp[0])
                          del_mmej_freq_small.append(temp[1])
                          del_mmej_freq_large.append(temp[2])
        if not del_mmej:
              for i in self.MH_lengths:
                    del_mmej_cand.append('')
                    del_mmejl.append(False)
              #      del_mmej_motif_pos.append([])
              distance_to_indel_position=np.nan
              del_mmej_cand_len=del_mmej_cand
              del_mmej_marked=del_mmej_cand
              del_mmej_marked_on_refa=del_mmej_cand
              del_mmej_freq_small=del_mmej_cand
              del_mmej_freq_large=del_mmej_cand
              del_mmej_motif_pos=del_mmej_cand
              #del_mmej_cand,mmej_cand= np.nan, np.nan
              self.distance_to_indel_position=np.nan
        else:
              for zz in range(len(self.MH_lengths)-len(mmej_cand)):
                    mmej_cand.insert(0, '')     
                    del_mmej_cand.insert(0,'')
                    del_mmej_cand_len.insert(0,'')
                    del_mmej_marked.insert(0,'')
                    del_mmej_marked_on_ref.insert(0,'')
                    del_mmej_motif_pos.insert(0,'')
#                    del_mmej_motif_pos.insert(0,[])
                    del_mmej_freq_small.insert(0,'')
                    del_mmej_freq_large.insert(0,'')  
                    del_mmejl.insert(0,False)  
              #      del_last_dimer.insert(0,[])

        del_mmej_cand=str(del_mmej_cand)
        del_mmej_cand_len=str(del_mmej_cand_len)
        del_mmej_marked=str(del_mmej_marked)
        del_mmej_marked_on_ref=str(del_mmej_marked_on_ref)
        del_last_dimer=str(del_last_dimer)
        del_mmej_motif_pos=str(del_mmej_motif_pos)
        del_mmej_freq_small=str(del_mmej_freq_small)
        del_mmej_freq_large=str(del_mmej_freq_large)
        del_mmejl=str(del_mmejl)
#        print(mmej_cand)
#        if not checkNaN(mmej_cand):
#            del_mmej_motif_pos = np.array([m.end() for m in #positions are the end positions of the motif in the sequence (e.g. TA if sequence is TATATA haso indices 2,4,6 (1based))
#            reg.finditer(mmej_cand, self.ref_genome_down, overlapped=False)])
#        else:
#            del_mmej_motif_pos = ""
        _d = {
            'del_mmej': del_mmej, 'del_mmejl':del_mmejl, 'del_mmej_cand': del_mmej_cand ,
            'del_mmej_marked': del_mmej_marked ,'del_mmej_marked_on_ref':del_mmej_marked_on_ref,
            'del_last_dimer':del_last_dimer, 'del_mmej_cand_len':del_mmej_cand_len,
            'del_mmej_motif_pos':del_mmej_motif_pos,'del_mmej_freq_small':del_mmej_freq_small,
            'del_mmej_freq_large':del_mmej_freq_large 
            }
        
        self.del_out_dict = _d         
    
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
            'SD_ID_mutant_pattern': the genomic context after the repair mechanism
            'SD_ID_MHrevc': the revc sequence of MH (MH2)
            'SD_ID_Prevc': the revc sequence of P (P2)
            'SD_ID_repeat_pat': the 2nd pattern found downstream the DSB
            'SD_ID_repeat_pat_len': SD_ID_repeat_pat's length
            'SD_ID_last_dimer': SD_ID_repeat_pat's last dimer useful for hidden markov
            'SD_ID_dist_between_reps': distance between the patterns (Y sequence)
            'SD_ID_motif_pos' : position of the pattern from the position of the DSB (pos)
            'SD_ID_motif_freq_small' : frequence of the pattern in small window
            'SD_ID_motif_freq_large' : frequence of the pattern in big window       
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
        SD_ID_motif_pos, SD_ID_motif_freq_large, SD_ID_motif_freq_small = list(), list(), list()               
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
            if MH1_P1revc in self.DSB_down[(len(self.INDEL) + len(P1) + 1):]:
                SD_inverted_deletion = True

            else:
                SD_inverted_deletion = False
                #print("BR not found")
            
            if SD_inverted_deletion:
                #print("BR found!!!")
            
                
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
        
                    #sequence between the pattern (from P1 to MH2)
                    mmej_marked_inter_reps_seq = f'{self.DSB_down[len(P2) + self.indel_length :(len(P2) + self.indel_length + rep_pat_pos_2)]}'
                    #SD_inverted_deletion_mutant_pattern is useful to illustrate the genomic context in the output file
                    SD_inverted_deletion_mutant_pattern = self.inverted_deletions_mutant_pattern_generator(P1=P1, MH1=MH1,rep_pat=rep_pat)
                    
                    #print(mmej_marked_inter_reps_seq)
                    #print(SD_inverted_deletion_mutant_pattern)

                    #defining SD_ID_motif_pos,SD_ID_motif_freq_small,SD_ID_motif_freq_large
                    temp_SD_inverted_deletion=get_motifs_freqs(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, large_window=1000,small_window=self.windowsize,motif=rep_pat,indel_type='DEL')
                    
                    #print(temp_SD_inverted_deletion)
                    
                    #SD_ID_motif_pos.append(temp_SD_inverted_deletion[0])
                    SD_ID_motif_freq_small.append(temp_SD_inverted_deletion[1])
                    SD_ID_motif_freq_large.append(temp_SD_inverted_deletion[2])
                    
                    # +1 and - 1 because 1 nt is necessary due to the fomration of the loop in snap_back mechanism
                    up_start= len(self.INDEL) + len(P1) + len (rep_pat) + 1
                    down_start = - len(MH1) - len (rep_pat) - 1
                    print(up_start)
                    print(down_start)

                    #SD_ID_motif_pos_2=get_motifs_pos_2(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)
                    
                    #print(SD_ID_motif_pos)
                    SD_ID_motif_pos_3=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)

                   #       
                    SD_inverted_deletion_P2 = P2
                    SD_inverted_deletion_MH2 = MH2
                    SD_inverted_deletion_repeat_pat = rep_pat
                    dist_between_reps = len(mmej_marked_inter_reps_seq)
                    SD_inverted_deletion_repeat_pat_len=len(SD_inverted_deletion_repeat_pat)
                    SD_inverted_deletion_last_dimer = SD_inverted_deletion_repeat_pat[-2:]
                
                    # setting output dict
                    _d = {
                        'SD_inverted_deletion': SD_inverted_deletion,
                        'SD_ID_mutant_pattern':SD_inverted_deletion_mutant_pattern,
                        'SD_ID_MHrevc': SD_inverted_deletion_MH2,
                        'SD_ID_Prevc': SD_inverted_deletion_P2, 
                        'SD_ID_repeat_pat': SD_inverted_deletion_repeat_pat,
                        'SD_ID_repeat_pat_len':SD_inverted_deletion_repeat_pat_len,
                        'SD_ID_last_dimer': SD_inverted_deletion_last_dimer,
                        'SD_ID_dist_between_reps': dist_between_reps,
                        'SD_ID_motif_pos' : SD_ID_motif_pos_3,
                        'SD_ID_motif_freq_small' : SD_ID_motif_freq_small,
                        'SD_ID_motif_freq_large' : SD_ID_motif_freq_large
                        }
                    self.snap_out_dict = _d
                    
              
                if self.extension:
                      print("ON")
                      rep_pat = MH1_P1revc # defining rep_pat and P1_1 for the extarnal else
                      P1_1 = self.reverse_complement_converter(seq = P1)
                      print("rep_pat: "+rep_pat)
                      rep_pat_pos_2 = self.mutant_sequence[
                              (self.indel_position + len(self.INDEL) + len (P1))+1:].index(rep_pat)
                      max_elongation_p = rep_pat_pos_2 // 2

                            #distance from the end of SD and start of SDrevc
                      print("distance from the end of SD and start of SDrevc: " + str(rep_pat_pos_2))
        
                      for i in range (0, self.windowsize):
                      ###elongation of MH, from 0 to windowsize beacuse we can extend Mh at the max, without limits, according to the pattern 
                          MH1_1 = self.ref_genome_up [(- MH_lengths -i):] 
                          #P1_1 = self.ref_genome_down [(self.indel_length):(self.indel_length+MH_lengths+i)]
                          MH1_1_P1 = MH1_1 + P1
                          #define the revc in order to find that
                          MH1_1_P1revc = self.reverse_complement_converter(seq = MH1_1_P1)
                          #print("MH1_1_P1: "+MH1_1_P1)
                          #print("MH1_1_P1revc: "+MH1_1_P1revc)
                          #print("self.INDEL :" + self.INDEL)
                          #print("MH1_1: "+ MH1_1)
                          #print("P1: "+ P1)
                          #searching MH1_1_P1revc in self.ref_genome_down
                          if MH1_1_P1revc in self.ref_genome_down [(len(self.INDEL) + len(P1)):]:
                              tmp_pos = self.ref_genome_down[(len(self.INDEL) + len(P1)):].index(MH1_1_P1revc)
                              #tmp_pos è la posizione del secondo template
                              #print("tmp_pos: "+str(tmp_pos))
        
                          else:
                              MH1_2 = self.ref_genome_up[(- MH_lengths -i +1):]
                              print("elongations finished")
                              #print("MH1_2: "+str(MH1_2))
                              MH1_2_P1 = MH1_2 + P1
                              #print("MH1_2_P1_2: " + MH1_2_P1)
                              MH1_2_P1revc = self.reverse_complement_converter(seq = MH1_2_P1)
                              #print("MH1_2_P1_2revc: " + MH1_2_P1revc)
                              tmp_pos = self.ref_genome_down[
                                      (len(self.INDEL) + len(P1)):].index(MH1_2_P1revc)
                              #print("tmp_pos" +str(tmp_pos))
                              break
                          #if tmp_pos < (len(self.INDEL) + len(P1_2)): break
                          #why i should break the loop?
                      

                      # elongation of P, from 1 to max_elongation_p becasue 1 is the 1 nt necessary fro the formation of the loop,
                      # and max_elongation_p is the distance from 1° and 2° SD / 2,
                      for i in range(1, max_elongation_p):
                          P1_1 = self.ref_genome_down[(self.indel_length):(self.indel_length+MH_lengths+i)]
                          MH1_2_P1_1 = MH1_2 + P1_1
                          MH1_2_P1_1revc = self.reverse_complement_converter(seq = MH1_2_P1_1)
                                                    
                          if MH1_2_P1_1revc in self.ref_genome_down[(len(self.INDEL) + len(P1_1)):]:
                              continue
                          else:
                              P1_2 = P1_1[:-1]
                              #print("P1_2: "+P1_2)
                              MH1_2_P1_2 = MH1_2 + P1_2
                              #print("MH1_2_P1_2: "+MH1_2_P1_2)
                              MH1_2_P1_2revc = self.reverse_complement_converter(seq = MH1_2_P1_2)
                              #print("MH1_2_P1_2revc: "+MH1_2_P1_2revc)
                              break
                      else:
                          MH1_2_P1_2revc = MH1_2_P1revc #se P2 non viene allungato allora il P-MH rimane con P non allungato, 
                          P1_2 = P1_1  #P non si allunga se la distanza fra i due P centrali è 1
                      
                      rep_pat = MH1_2_P1_2revc
                      rep_pat_pos = self.mutant_sequence[
                          (self.indel_position+len(P1_2)):].index(rep_pat)
                      #print("rep_pat_pos"+str (rep_pat_pos)) #distanza fra P1 e MH2
                      #seq between the pattern
                      mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P1_2) + self.indel_length):(len(P1_2) + rep_pat_pos)]}'
                        #seq from P1 to MH2
                      #print(mmej_marked_inter_reps_seq)
                      #SD_inverted_deletion_mutant_pattern is useful to illustrate the genomic context in the output file
                      
                      SD_inverted_deletion_mutant_pattern = self.inverted_deletions_mutant_pattern_generator(P1=P1_2, MH1=MH1_2,rep_pat=rep_pat)
                      #print(SD_inverted_deletion_mutant_pattern)
                      
                      #define SD_ID_motif_pos,SD_ID_motif_freq_small, SD_ID_motif_freq_large
                      temp_SD_inverted_deletion=get_motifs_freqs(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, large_window=1000,small_window=self.windowsize,motif=rep_pat,indel_type='DEL')
                      #print(temp_SD_inverted_deletion)
                            
                    
                      
                      #SD_ID_motif_pos.append(temp_SD_inverted_deletion[0])
                      SD_ID_motif_freq_small.append(temp_SD_inverted_deletion[1])
                      SD_ID_motif_freq_large.append(temp_SD_inverted_deletion[2])

                      #-1 and + 1 because 1 nt is necessary due to the formation of the loop
                      up_start= len(self.INDEL) + len(P1_2) + len(rep_pat) + 1
                      down_start = - len(MH1_2) - len(rep_pat)
                      print("rep_pat: "+rep_pat)
                      print("up_start: "+str(up_start))

                      print("down_start : "+str (down_start))

                      #SD_ID_motif_pos_2=get_motifs_pos_2(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)
                      
                      SD_ID_motif_pos_3=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)

                      SD_inverted_deletion_P2 = self.reverse_complement_converter(seq = P1_2)
                      SD_inverted_deletion_MH2 = self.reverse_complement_converter(seq = MH1_2)
                      SD_inverted_deletion_repeat_pat = rep_pat
                      dist_between_reps = len(mmej_marked_inter_reps_seq)
                      SD_inverted_deletion_repeat_pat_len=len(SD_inverted_deletion_repeat_pat)
                      SD_inverted_deletion_last_dimer = SD_inverted_deletion_repeat_pat[-2:]
                      
                  
                      
                      # setting output dict
                      _d = {
                          'SD_inverted_deletion': SD_inverted_deletion,
                          'SD_ID_mutant_pattern':SD_inverted_deletion_mutant_pattern,
                          'SD_ID_MHrevc': SD_inverted_deletion_MH2,
                          'SD_ID_Prevc': SD_inverted_deletion_P2, 
                          'SD_ID_repeat_pat': SD_inverted_deletion_repeat_pat,
                          'SD_ID_repeat_pat_len':SD_inverted_deletion_repeat_pat_len,
                          'SD_ID_last_dimer': SD_inverted_deletion_last_dimer,
                          'SD_ID_dist_between_reps': dist_between_reps,
                          'SD_ID_motif_pos' : SD_ID_motif_pos_3,
                          'SD_ID_motif_freq_small' : SD_ID_motif_freq_small,
                          'SD_ID_motif_freq_large' : SD_ID_motif_freq_large
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
            'SD_II_mutant_pattern': the genomic context after the repair mechanism
            'SD_II_MHrevc': the revc sequence of MH (MH2)
            'SD_II_Prevc': the revc sequence of P (P2)
            'SD_II_repeat_pat': the 2nd pattern found downstream the DSB
            'SD_II_repeat_pat_len': SD_II_repeat_pat's length
            'SD_II_last_dimer': SD_II_repeat_pat's last dimer useful for hidden markov
            'SD_II_dist_between_reps': distance between the patterns (Y sequence)
            'SD_II_motif_pos' : position of the pattern from the position of the DSB (pos)
            'SD_II_motif_freq_small' : frequence of the pattern in small window
            'SD_II_motif_freq_large' : frequence of the pattern in big window       
            No extension: MH2 is the revc of MH found downstream
                          P2 is the revc of MH found downstream
            Extension: MH1 is MH not extended
                       MH1_1 is MH during the extension
                       MH1_2 is MH with max extension
                       
                       P1 is P not extended
                       P1_1 is P during extension
                       P1_2 is P with max extension 

            
        """
        SD_II_motif_pos, SD_II_motif_freq_large, SD_II_motif_freq_small = list(), list(), list() #creation of the lists               
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
                    print("OFF")
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
                
                    #print("Prevc: "+P2)
                    #print("MHrevc: "+MH2)

                    # creating Crearing mmej_marked
                    mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P2) + self.indel_length):(self.indel_length +len(P2) + rep_pat_pos)]}'
                    print("snapback - loop :"+mmej_marked_inter_reps_seq)
        
                    # SD_I_Insertion_mutant_pattern is useful to illustrate the genomic context in the output file
                    SD_I_Insertion_mutant_pattern = self.snap_mutant_pattern_generator(P1=P1,P2=P2,MH1=MH1,
                            MH2=MH2,rep_pat=rep_pat)
                    #print (SD_I_Insertion_mutant_pattern)
                                 
                    
                    #defining SD_II_motif_pos, SD_II_motif_freq_small,SD_II_motif_freq_large 
                    temp_SD_inverted_insertion=get_motifs_freqs(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, large_window=1000,small_window=self.windowsize,motif=rep_pat,indel_type='DEL')
                    #print (temp_SD_inverted_insertion)
          
                    SD_II_motif_pos.append(temp_SD_inverted_insertion[0])
                    SD_II_motif_freq_small.append(temp_SD_inverted_insertion[1])
                    SD_II_motif_freq_large.append(temp_SD_inverted_insertion[2])
            
                    SD_inverted_insertion_P2 = self.reverse_complement_converter(seq = P2)
                    SD_inverted_insertion_MH2 = self.reverse_complement_converter(seq = MH2)
                    SD_inverted_insertion_repeat_pat = rep_pat
                    dist_between_reps = len(mmej_marked_inter_reps_seq)
                    SD_inverted_insertion_repeat_pat_len=len(SD_inverted_insertion_repeat_pat)
                    SD_inverted_insertion_last_dimer = SD_inverted_insertion_repeat_pat[-2:]
                      
                  
                      
                      # setting output dict
                    _d = {
                          'SD_inverted_insertion': SD_inverted_insertion,
                          'SD_II_mutant_pattern':SD_I_Insertion_mutant_pattern,
                          'SD_II_MHrevc': SD_inverted_insertion_MH2,
                          'SD_II_Prevc': SD_inverted_insertion_P2, 
                          'SD_II_repeat_pat': SD_inverted_insertion_repeat_pat,
                          'SD_II_repeat_pat_len':SD_inverted_insertion_repeat_pat_len,
                          'SD_II_last_dimer': SD_inverted_insertion_last_dimer,
                          'SD_II_dist_between_reps': dist_between_reps,
                          'SD_II_motif_pos' : SD_II_motif_pos,
                          'SD_II_motif_freq_small' : SD_II_motif_freq_small,
                          'SD_II_motif_freq_large' : SD_II_motif_freq_large
                             }
                    self.snap_out_dict = _d
                
                if self.extension:
                      print("ON")
                      rep_pat = II_pattern
                      rep_pat_pos_2 = self.mutant_sequence[
                        (self.indel_position):].index(rep_pat)
                      print(rep_pat_pos_2)
                      #distance between fra pos e Prevc not elongated
                      
                      for i in range (1, rep_pat_pos_2):
                      ###elongation of MH 
                          MH1_1 = self.ref_genome_up [(- MH_lengths -i):] 
                          MH1_1_indel_P1 = MH1_1 + self.INDEL + P1
                          MH1_1_indel_P1revc = self.reverse_complement_converter(seq = MH1_1_indel_P1)
                          #print("MH1_1_indel_P1revc: "+ MH1_1_indel_P1revc)
                          #print("self.INDEL :" + self.INDEL)
                          #print("MH1_1: "+ MH1_1)
                          #print("P1: "+ P1)
                          if MH1_1_indel_P1revc in self.ref_genome_down [ len(P1):]:
                              tmp_pos = self.ref_genome_down[(len(P1)):].index(MH1_1_indel_P1revc)
                              #tmp_pos pos of 2nd template during elongation 
                              #print("2nd template pos: "+str(tmp_pos))
                      
                          else:
                              MH1_2 = self.ref_genome_up[(- MH_lengths -i +1):]
                              #print("elongations finished")
                              #print("MH1_2: "+str(MH1_2))
                              MH1_2revc = self.reverse_complement_converter(seq = MH1_2)
                              MH1_2_indel_P1 = MH1_2 + self.INDEL + P1
                              #print("MH1_2_indel_P1" +MH1_2_indel_P1)
                              MH1_2_indel_P1revc = self.reverse_complement_converter(seq = MH1_2_indel_P1)
                              #print("MH1_2_indel_P1revc: " + MH1_2_indel_P1revc)
                              tmp_pos_2 = self.ref_genome_down[
                                       len(P1):].index(MH1_2_indel_P1revc)
                              #print("tmp_pos_2: " +str(tmp_pos_2))
                              break
                          #if tmp_pos < (len(self.INDEL) + len(P1_2)): break
                          #why i sholud break the loop?
                      
                      #elongation of P
                      for i in range(1, rep_pat_pos_2):
                          P1_1 = self.ref_genome_down[:( + MH_lengths + i)]
                          MH1_2_indel_P1_1 = MH1_2 + self.INDEL + P1_1
                          MH1_2_indel_P1_1revc = self.reverse_complement_converter(seq = MH1_2_indel_P1_1)
                          #print("MH1_2_indel_P1_1revc: " + MH1_2_indel_P1_1revc)

                          if MH1_2_indel_P1_1revc in self.ref_genome_down [ len(P1):]:
                              continue
                          else:
                              P1_2 = P1_1[:-1]
                              P1_2revc = self.reverse_complement_converter(seq = P1_2)
                              #print("P1_2: "+P1_2)
                              MH1_2_indel_P1_2 = MH1_2 + self.INDEL + P1_2
                              II_pattern = self.reverse_complement_converter(seq = MH1_2_indel_P1_2)
                              print("II_pattern : "+ II_pattern)
                              break
                      else:
                          II_pattern=MH1_2_indel_P1revc #se P2 non viene allungato allora il BR rimane con R non allungato
                      

                      rep_pat = II_pattern
                      print ("final pattern found: " + rep_pat)
                      rep_pat_pos = self.mutant_sequence[
                          (self.indel_position+len(P1_2)):].index(rep_pat)
                      #print("rep_pat_pos"+str (rep_pat_pos)) #distance between P1 e P1revc
                      
                      #sequence between the two patterns
                      mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P1_2) + self.indel_length):(self.indel_length +len(P1_2) + rep_pat_pos)]}'
                      #print("snapback - loop :"+mmej_marked_inter_reps_seq)
                      
                       
                      #SD_I_Insertion_mutant_pattern is useful to illustrate the genomic context in the output file
                      SD_I_Insertion_mutant_pattern = self.snap_mutant_pattern_generator(P1=P1_2,P2=P1_2revc,MH1=MH1_2,MH2=MH1_2revc,rep_pat=rep_pat)
                      #print (SD_I_Insertion_mutant_pattern)
                      #defining  SD_II_motif_pos,SD_II_motif_freq_small,SD_II_motif_freq_large 
                      temp_SD_inverted_insertion=get_motifs_freqs(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, large_window=1000,small_window=self.windowsize,motif=rep_pat,indel_type='DEL')
                      #print (temp_SD_inverted_insertion)
          
                      SD_II_motif_pos.append(temp_SD_inverted_insertion[0])
                      SD_II_motif_freq_small.append(temp_SD_inverted_insertion[1])
                      SD_II_motif_freq_large.append(temp_SD_inverted_insertion[2])
            
                      SD_inverted_insertion_P2 = self.reverse_complement_converter(seq = P1_2)
                      SD_inverted_insertion_MH2 = self.reverse_complement_converter(seq = MH1_2)
                      SD_inverted_insertion_repeat_pat = rep_pat
                      dist_between_reps = len(mmej_marked_inter_reps_seq)
                      SD_inverted_insertion_repeat_pat_len=len(SD_inverted_insertion_repeat_pat)
                      SD_inverted_insertion_last_dimer = SD_inverted_insertion_repeat_pat[-2:]
                      
                  
                      
                      # setting output dict
                      _d = {
                            'SD_inverted_insertion': SD_inverted_insertion,
                            'SD_II_mutant_pattern':SD_I_Insertion_mutant_pattern,
                            'SD_II_MHrevc': SD_inverted_insertion_MH2,
                            'SD_II_Prevc': SD_inverted_insertion_P2, 
                            'SD_II_repeat_pat': SD_inverted_insertion_repeat_pat,
                            'SD_II_repeat_pat_len':SD_inverted_insertion_repeat_pat_len,
                            'SD_II_last_dimer': SD_inverted_insertion_last_dimer,
                            'SD_II_dist_between_reps': dist_between_reps,
                            'SD_II_motif_pos' : SD_II_motif_pos,
                            'SD_II_motif_freq_small' : SD_II_motif_freq_small,
                            'SD_II_motif_freq_large' : SD_II_motif_freq_large
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
            'SD_DD_mutant_pattern': the genomic context after the repair mechanism
            'SD_DD_MHrevc': the revc sequence of MH (MH2)
            'SD_DD_Prevc': the revc sequence of P (P2)
            'SD_DD_repeat_pat': the 2nd pattern found downstream the DSB
            'SD_DD_repeat_pat_len': SD_DD_repeat_pat's length
            'SD_DD_last_dimer': SD_DD_repeat_pat's last dimer useful for hidden markov
            'SD_DD_dist_between_reps': distance between the patterns (Y sequence)
            'SD_DD_motif_pos' : position of the pattern from the position of the DSB (pos)
            'SD_DD_motif_freq_small' : frequence of the pattern in small window
            'SD_DD_motif_freq_large' : frequence of the pattern in big window       
            No extension: MH2 is the revc of MH found downstream
                          P2 is the revc of MH found downstream
            Extension: MH1 is MH not extended
                       MH1_1 is MH during the extension
                       MH1_2 is MH with max extension
                       
                       P1 is P not extended
                       P1_1 is P during extension
                       P1_2 is P with max extension

        
        returns:
            SD_direct_deletion: True/False
            SD_deletion_mutant_pattern : the seq after the SD-direct_deletion
            SD_DD_P2': P2 elongated,
            'SD_DD_MH2': MH2 elongated 
            'SD_DD_repeat_pat': MH2 + P2 
            'SD_DD_repeat_pat_len': length of SD_DD_repeat_pat,
            'SD_DD_last_dimer': last dimer of SD_DD_repeat_pat,
            'SD_DD_dist_between_reps': the length og 'random seq'
            'SD_DD_motif_pos' : positions of BR motif 
            'SD_DD_motif_freq_large' : motif frequence 
            
            
        """  
        SD_DD_motif_pos, SD_DD_motif_freq_large, SD_DD_motif_freq_small = list(), list(), list()               
        ###1st step is to bring MH1 and P1 from the ref seq 
        for MH_lengths in self.MH_lengths:
            MH1 = self.ref_genome_up [ - MH_lengths :]  
            P1 = self.ref_genome_down [(self.indel_length):(self.indel_length+MH_lengths)]
            MH1_P1 = MH1 + P1
            #print("MH1_P1: " + MH1_P1)
            #print("self.INDEL :" + self.INDEL)
            #print(len(self.INDEL))
            #print("MH1: "+ MH1)
            #print(len(MH1))
            #print("P1: "+ P1)
            #print(len(P1))
            #print("ex: " + str(self.extension))
            #searching the pattern in self.DSB_down
        
            if MH1_P1 in self.DSB_down[(len(self.INDEL) + len(P1) + 0):]:#1 nt is the necessary for the formation of the loop
                SD_direct = True

            else:
                SD_direct = False
                print ("BR not found")
            
            if SD_direct:
                print("BR found!")
            
                if not self.extension:
                    print("OFF")
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
                
                    #print("MH2: "+MH2)
                    #print("P2: "+P2)
                    #seq between the patterns
                    mmej_marked_inter_reps_seq = f'{self.DSB_down[len(P2) + self.indel_length :(len(P2) + self.indel_length + rep_pat_pos_2)]}'
                    #SD_direct_deletion_mutant_pattern is useful to illustrate the genomic context in the output file
                    SD_direct_deletion_mutant_pattern = self.direct_deletions_mutant_pattern_generator(P2=P2, MH2=MH2,rep_pat=rep_pat)
                    
                    #define SD_DD_motif_pos,SD_DD_motif_freq_small,SD_DD_motif_freq_large
                    temp_SD_direct_deletion=get_motifs_freqs(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, large_window=1000,small_window=self.windowsize,motif=rep_pat,indel_type='DEL')
                    
                    
                    #SD_DD_motif_pos.append(temp_SD_direct_deletion[0])
                    
                    #tutto ciò vale solo per dd in inverted devi mettere il + 1 idiota in up e down
                    up_start = len(self.INDEL) + len (P1) + len(rep_pat)
                    down_start = - len(MH1) - len(rep_pat)
                    print(up_start)
                    print(down_start)
                    print(self.pos_on_chr)
                    print(rep_pat)
                    #SD_DD_motif_pos_2=get_motifs_pos_2(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)
                    
                    #print(SD_DD_motif_pos_2)

                    SD_DD_motif_pos_3=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)
                    print(SD_DD_motif_pos_3)
                    #SD_DD_motif_pos_2.append(temp_SD_direct_deletion[0])
                    
                    #print(temp_SD_direct_deletion)
                    SD_DD_motif_freq_small.append(temp_SD_direct_deletion[1])
                    SD_DD_motif_freq_large.append(temp_SD_direct_deletion[2])

                   #       
                    SD_direct_deletion_P2 = P2
                    SD_direct_deletion_MH2 = MH2
                    SD_direct_deletion_repeat_pat = rep_pat
                    dist_between_reps = len(mmej_marked_inter_reps_seq)
                    SD_direct_deletion_repeat_pat_len=len(SD_direct_deletion_repeat_pat)
                    SD_direct_deletion_last_dimer = SD_direct_deletion_repeat_pat[-2:]
                
                    # setting output dict
                    _d = {
                        'SD_direct_deletion': SD_direct,
                        'SD_DD_mutant_pattern':SD_direct_deletion_mutant_pattern,
                        'SD_DD_MH2': SD_direct_deletion_MH2,
                        'SD_DD_P2': SD_direct_deletion_P2, 
                        'SD_DD_repeat_pat': SD_direct_deletion_repeat_pat,
                        'SD_DD_repeat_pat_len':SD_direct_deletion_repeat_pat_len,
                        'SD_DD_last_dimer': SD_direct_deletion_last_dimer,
                        'SD_DD_dist_between_reps': dist_between_reps,
                        'SD_DD_motif_pos' : SD_DD_motif_pos_3,
                        'SD_DD_motif_freq_small' : SD_DD_motif_freq_small,
                        'SD_DD_motif_freq_large' : SD_DD_motif_freq_large
                        }
                    self.loop_out_dict = _d
                
              
                if self.extension:
                      print("ON")
                      rep_pat = MH1_P1
                      rep_pat_pos_2 = self.mutant_sequence[
                          (self.indel_position + len(self.INDEL) + len (P1)):].index(rep_pat)#distanza fra pos e MH2
                      print("distance from the end of the SD anche the start of the 2° MH in SD_DD: " + str(rep_pat_pos_2))
                      MH1_1 = MH1
                      P1_1 = P1

                      for i in range (0, rep_pat_pos_2):
                      ###elongation of MH 
                          MH1_1 = self.ref_genome_up [(- MH_lengths -i):] 
                          P1_1 = self.ref_genome_down [(self.indel_length):(self.indel_length+MH_lengths+i)]
                          MH1_1_P1 = MH1_1 + P1
                          #print("MH1_1_P1: "+MH1_1_P1)
                          #print("self.INDEL :" + self.INDEL)
                          #print("MH1_1: "+ MH1_1)
                          #print("P1: "+ P1)
                          if MH1_1_P1 in self.ref_genome_down [(len(self.INDEL) + len(P1)):]:
                              tmp_pos = self.ref_genome_down[(len(self.INDEL) + len(P1)):].index(MH1_1_P1)
                              #tmp_pos è la posizione del secondo template
                              #print("tmp_pos: "+str(tmp_pos))

                          else:
                              MH1_2 = self.ref_genome_up[(- MH_lengths -i +1):]
                              #print("elongations finished")
                              #print("MH1_2: "+str(MH1_2))
                              MH1_2_P1 = MH1_2 + P1
                              #print("MH1_2_P1_2: " + MH1_2_P1)
                              tmp_pos = self.ref_genome_down[
                                      (len(self.INDEL) + len(P1)):].index(MH1_2_P1)
                              #print("tmp_pos" +str(tmp_pos))
                              break
                      else:
                          MH1_2 = MH1_1     #se rep_pat_pos_2 è piccolo e il for fa poche iterazione, posso non arrivare al break e quindi posso anche
                          MH1_2_P1 = rep_pat #saltare l'else interno, per questo claude mi ha detto di mettere un else esterno
                                                
                      # elongation of P
                      for i in range(0, rep_pat_pos_2):
                          P1_1 = self.ref_genome_down[(self.indel_length):(self.indel_length+MH_lengths+i)]
                          MH1_2_P1_1 = MH1_2 + P1_1

                          if MH1_2_P1_1 in self.ref_genome_down[(len(self.INDEL) + len(P1_1)):]:
                              continue
                          else:
                              P1_2 = P1_1[:-1]
                              #print("P1_2: "+P1_2)
                              MH1_2_P1_2 = MH1_2 + P1_2
                              #print("MH1_2_P1_2: "+MH1_2_P1_2)
                              break
                      else:
                          MH1_2_P1_2=rep_pat #se P2 non viene allungato allora il BR rimane con R non allungato
                          P1_2 = P1_1 

                      rep_pat = MH1_2_P1_2
                      rep_pat_pos = self.mutant_sequence[
                          (self.indel_position+len(P1_2)):].index(rep_pat)
                      #print("rep_pat_pos"+str (rep_pat_pos)) #distanza fra P1 e MH2
                    
                      #seq from from P1 to MH2
                      mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P1_2) + self.indel_length):(len(P1_2) + rep_pat_pos)]}'
                      #print(mmej_marked_inter_reps_seq)

                      #SD_direct_deletion_mutant_pattern is useful to illustrate the genomic context in the output file
                      SD_direct_deletion_mutant_pattern = self.direct_deletions_mutant_pattern_generator(P2=P1_2, MH2=MH1_2,rep_pat=rep_pat)
                      #define SD_DD_motif_pos,SD_DD_motif_freq_small,SD_DD_motif_freq_large                      
                      temp_SD_direct_deletion=get_motifs_freqs(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, large_window=1000,small_window=self.windowsize,motif=rep_pat,indel_type='DEL')
                    
                    
                      
                      #SD_DD_motif_pos.append(temp_SD_direct_deletion[0])
                      SD_DD_motif_freq_small.append(temp_SD_direct_deletion[1])
                      SD_DD_motif_freq_large.append(temp_SD_direct_deletion[2])

                      
                    #tutto ciò vale solo per dd in inverted devi mettere il + 1 idiota in up e down
                      up_start = len(self.INDEL) + len (P1_2) + len(rep_pat)
                      down_start = - len(MH1_2) - len(rep_pat) 
                      #print(up_start)
                      #print(down_start)
                      #print("pos on chrom: "+self.pos_on_chr)
              
                      #SD_DD_motif_pos=get_motifs_pos_2(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)
                    
                      #print(SD_DD_motif_pos)
                      SD_DD_motif_pos_3=get_motifs_pos_3(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr,motif=rep_pat,windowsize=self.windowsize,down_start=down_start,up_start=up_start)
                      
                      SD_direct_deletion_P2 = P1_2
                      SD_direct_deletion_MH2 = MH1_2
                      SD_direct_deletion_repeat_pat = rep_pat
                      dist_between_reps = len(mmej_marked_inter_reps_seq)
                      SD_direct_deletion_repeat_pat_len=len(SD_direct_deletion_repeat_pat)
                      SD_direct_deletion_last_dimer = SD_direct_deletion_repeat_pat[-2:]
                      
                      
                      
                      # setting output dict
                      _d = {
                          'SD_direct_deletion': SD_direct,
                          'SD_DD_mutant_pattern':SD_direct_deletion_mutant_pattern,
                              'SD_DD_MH2': SD_direct_deletion_MH2,
                              'SD_DD_P2': SD_direct_deletion_P2, 
                              'SD_DD_repeat_pat': SD_direct_deletion_repeat_pat,
                              'SD_DD_repeat_pat_len':SD_direct_deletion_repeat_pat_len,
                              'SD_DD_last_dimer': SD_direct_deletion_last_dimer,
                              'SD_DD_dist_between_reps': dist_between_reps,
                              'SD_DD_motif_pos' : SD_DD_motif_pos_3,
                              'SD_DD_motif_freq_small' : SD_DD_motif_freq_small,
                              'SD_DD_motif_freq_large' : SD_DD_motif_freq_large
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
            SD_DI_mutant_pattern : the seq after the SD-direct_deletion
            SD_DI_P2': P2 elongated,
            'SD_DI_MH2': MH2 elongated 
            'SD_DI_repeat_pat': MH2 + indel + P2 
            'SD_DI_repeat_pat_len': length of SD_DI_repeat_pat,
            'SD_DI_last_dimer': last dimer of SD_DI_repeat_pat,
            'SD_DI_dist_between_reps': the length of 'random seq'
            'SD_DI_motif_pos' : positions of BZR motif 
            'SD_DI_motif_freq_large' : motif frequence 
            
            
        """
        SD_DI_motif_pos, SD_DI_motif_freq_large, SD_DI_motif_freq_small = list(), list(), list()               
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
        
            if DI_pattern in self.ref_genome_down[(len(P1) + 1):]:#1 nt is the necessary for the formation of the loop
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
                
                    #print("MH2: "+MH2)
                    #print("P2: "+P2)

                    # Crearing mmej_marked
                    mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P2) + self.indel_length):(self.indel_length +len(P2) + rep_pat_pos)]}'
                    print("Y seq :"+mmej_marked_inter_reps_seq)
        
                    # set output variables as attributes
                    SD_D_Insertion_mutant_pattern = self.loop_mutant_pattern_generator(P2=P2, MH2=MH2,
                                 rep_pat=rep_pat)
                    #print (SD_D_Insertion_mutant_pattern)
                                 
                    SD_DI_P2 = P2
                    SD_DI_MH2 = MH2
                    SD_DI_repeat_pat = rep_pat
                    SD_DI_between_reps = len(mmej_marked_inter_reps_seq)
                    SD_DI_last_dimer = SD_DI_repeat_pat[-2:]
            
                    temp_SD_direct_insertion=get_motifs_freqs(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, large_window=1000,small_window=self.windowsize,motif=rep_pat,indel_type='DEL')
                    #print (temp_SD_direct_insertion)
            
                    SD_DI_motif_pos.append(temp_SD_direct_insertion[0])
                    SD_DI_motif_freq_small.append(temp_SD_direct_insertion[1])
                    SD_DI_motif_freq_large.append(temp_SD_direct_insertion[2])
            
       
                    #setting output dict
                    _d = {
                        'SD_direct_insertion': SD_direct_insertion, 
                        'SD_DI_mutant_pattern':SD_D_Insertion_mutant_pattern,
                        'SD_DI_P2': SD_DI_P2,'SD_DI_MH2': SD_DI_MH2, 
                        'SD_DI_repeat_pat': SD_DI_repeat_pat,
                        'SD_DI_dist_between_reps': SD_DI_between_reps,
                        'SD_DI_last_dimer': SD_DI_last_dimer,
                        'SD_DI_motif_pos': SD_DI_motif_pos,
                        'SD_DI_motif_freq_small': SD_DI_motif_freq_small,
                        'SD_DI_motif_freq_large': SD_DI_motif_freq_large
                         }
                    self.loop_out_dict = _d 
                
                
                if self.extension:
                      print("ON")
                      rep_pat = DI_pattern
                      rep_pat_pos_2 = self.mutant_sequence[
                        (self.indel_position):].index(rep_pat)
                      #distanza fra pos e MH2
                      

                      for i in range (1, rep_pat_pos_2):
                      ###elongation of MH 
                          MH1_1 = self.ref_genome_up [(- MH_lengths -i):] 
                          #P1_1 = self.ref_genome_down [(self.indel_length):(self.indel_length+MH_lengths+i)]
                          #MH1_1_P1 = MH1_1 + P1
                          MH1_1_indel_P1 = MH1_1 + self.INDEL + P1
                          #print("MH1_1_P1: "+MH1_1_P1)
                          #print("self.INDEL :" + self.INDEL)
                          #print("MH1_1: "+ MH1_1)
                          #print("P1: "+ P1)
                          if MH1_1_indel_P1 in self.ref_genome_down [ len(P1):]:
                              tmp_pos = self.ref_genome_down[(len(P1)):].index(MH1_1_indel_P1)
                              #tmp_pos è la posizione del secondo template una volta allungato
                              #print("tmp_pos: "+str(tmp_pos))

                          else:
                              MH1_2 = self.ref_genome_up[(- MH_lengths -i +1):]
                              #print("elongations finished")
                              #print("MH1_2: "+str(MH1_2))
                              #MH1_2_P1 = MH1_2 + P1
                              MH1_2_indel_P1 = MH1_2 + self.INDEL + P1
                              print("MH1_2_indel_P1_2: " + MH1_2_indel_P1)
                              tmp_pos_2 = self.ref_genome_down[
                                       len(P1):].index(MH1_2_indel_P1)
                              #print("tmp_pos_2: " +str(tmp_pos_2))
                              break
                          #if tmp_pos < (len(self.INDEL) + len(P1_2)): break
                          #why i sholud break the loop?
                        
                      # elongation of P2, confined in extention_space[:*starting position of MH2_2*]
                      for i in range(1, rep_pat_pos_2):
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
                          DI_pattern=MH1_2_P1 #se P2 non viene allungato allora il BR rimane con R non allungato
                
                      # Creating mmej_marked
                      rep_pat = DI_pattern
                      rep_pat_pos_3 = self.mutant_sequence[
                          (self.indel_position+self.indel_length+len(P1_2)):].index(rep_pat)
                      #print("rep_pat_pos_3: "+str (rep_pat_pos_3)) #distanza fra P1 e MH2
                    
                      mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P1_2) + self.indel_length):(len(P1_2) + self.indel_length + rep_pat_pos_3)]}'
                        #seq from P1 to MH2
                      print(mmej_marked_inter_reps_seq)
                
                      SD_D_Insertion_mutant_pattern = self.loop_mutant_pattern_generator(P2=P1_2, MH2=MH1_2,rep_pat=rep_pat)
                      #print(SD_D_Insertion_mutant_pattern)
                      
                      SD_DI_P2 = P1_2
                      SD_DI_MH2 = MH1_2
                      SD_DI_repeat_pat = rep_pat
                      SD_DI_between_reps = len(mmej_marked_inter_reps_seq)
                      SD_DI_last_dimer = SD_DI_repeat_pat[-2:]

                      temp_SD_direct_insertion=get_motifs_freqs(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, large_window=1000,small_window=self.windowsize,motif=rep_pat,indel_type='DEL')
                      #print(temp_SD_direct_insertion)
            
                      SD_DI_motif_pos.append(temp_SD_direct_insertion[0])
                      SD_DI_motif_freq_small.append(temp_SD_direct_insertion[1])
                      SD_DI_motif_freq_large.append(temp_SD_direct_insertion[2])
            
       
                    #setting output dict
                      _d = {
                          'SD_direct_insertion': SD_direct_insertion, 
                          'SD_DI_mutant_pattern':SD_D_Insertion_mutant_pattern,
                          'SD_DI_P2': SD_DI_P2,'SD_DI_MH2': SD_DI_MH2, 
                          'SD_DI_repeat_pat': SD_DI_repeat_pat,
                          'SD_DI_dist_between_reps': SD_DI_between_reps,
                          'SD_DI_last_dimer': SD_DI_last_dimer,
                          'SD_DI_motif_pos': SD_DI_motif_pos,
                          'SD_DI_motif_freq_small': SD_DI_motif_freq_small,
                          'SD_DI_motif_freq_large': SD_DI_motif_freq_large
                           }
                      self.loop_out_dict = _d 

    def sd_direct_substitution(self):
        print("self.INDEL :"+self.INDEL)
        print("self.ANC:"+self.ANC)
        print("self.DER:"+self.DER)
        """
        """

        SD_DS_motif_pos, SD_DS_motif_freq_large, SD_DS_motif_freq_small = list(), list(), list()               
        ###1st step is to bring MH1 and P1 from the ref seq 
        for MH_lengths in self.MH_lengths:
            MH1 = self.ref_genome_up [ - MH_lengths :]  
            P1 = self.ref_genome_down [ len(self.ANC) : len(self.ANC) + MH_lengths ]
            DS_pattern = MH1 + self.INDEL + P1 
            print("DS_pattern: " + DS_pattern)
            #print("MH1_P1: " + MH1_P1)
            print("self.INDEL :" + self.INDEL)
            print("MH1: "+ MH1)
            print("P1: "+ P1)
            #print("ex: " + str(self.extension))
            
            if DS_pattern in self.ref_genome_down [len(self.ANC) + 1:]:
                SD_direct_substitution = True
            else:
                SD_direct_substitution = False 
                print("BZR not found")

            if SD_direct_substitution:
                print("BZR found")
                
                if not self.extension:
                    #MH1 = non allungato, MH1_1 in fase di allungamento, MH1_2 =  MH1 a fine allungamento 
                    print("OFF")
                    rep_pat = DS_pattern
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
                
                    print("MH2: "+MH2)
                    print("P2: "+P2)
                    
                    mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(self.DER))+ len(P2) : (len(self.DER))+ len(P2) + rep_pat_pos]}'
                    print("Y seq :"+mmej_marked_inter_reps_seq)
        
                    
                    # set output variables as attributes
                    SD_D_Substitution_mutant_pattern = self.direct_substitution_mutant_pattern_generator(P2=P2, MH2=MH2,rep_pat=rep_pat)
                    print (SD_D_Substitution_mutant_pattern)
                                 
                    SD_DS_P2 = P2
                    SD_DS_MH2 = MH2
                    SD_DS_repeat_pat = rep_pat
                    SD_DS_repeat_pat_len = len (rep_pat)
                    SD_DS_between_reps = len(mmej_marked_inter_reps_seq)
                    SD_DS_last_dimer = SD_DS_repeat_pat[-2:]
            
                    temp_SD_direct_substitution=get_motifs_freqs(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, large_window=1000,small_window=self.windowsize,motif=rep_pat,indel_type='DEL')
                    print (temp_SD_direct_substitution)
            
                    SD_DS_motif_pos.append(temp_SD_direct_substitution[0])
                    SD_DS_motif_freq_small.append(temp_SD_direct_substitution[1])
                    SD_DS_motif_freq_large.append(temp_SD_direct_substitution[2])
            
       
                    #setting output dict
                    _d = {
                        'SD_direct_substitution': SD_direct_substitution, 
                        'SD_DS_mutant_pattern':SD_D_Substitution_mutant_pattern,
                        'SD_DS_P2': SD_DS_P2,'SD_DS_MH2': SD_DS_MH2, 
                        'SD_DS_repeat_pat': SD_DS_repeat_pat,
                        'SD_DS_repeat_pat_len': SD_DS_repeat_pat_len,
                        'SD_DS_dist_between_reps': SD_DS_between_reps,
                        'SD_DS_last_dimer': SD_DS_last_dimer,
                        'SD_DS_motif_pos': SD_DS_motif_pos,
                        'SD_DS_motif_freq_small': SD_DS_motif_freq_small,
                        'SD_DS_motif_freq_large': SD_DS_motif_freq_large
                         }
                    self.loop_out_dict = _d 
                

                if self.extension:
                    print("ON")
                    rep_pat = DS_pattern
                    rep_pat_pos_2 = self.mutant_sequence[
                            (self.indel_position):].index(rep_pat)
                    #distanza fra pos e MH2
                      

                    for i in range (1, rep_pat_pos_2):
                      ###elongation of MH 
                         MH1_1 = self.ref_genome_up [(- MH_lengths -i):] 
                         MH1_1_indel_P1 = MH1_1 + self.INDEL + P1
                         print("self.INDEL :" + self.INDEL)
                         print("MH1_1: "+ MH1_1)
                         print("P1: "+ P1)
                         if MH1_1_indel_P1 in self.ref_genome_down [ len(self.ANC) + len(P1):]:
                             tmp_pos = self.ref_genome_down[ len(self.ANC) + (len(P1)):].index(MH1_1_indel_P1)
                             #tmp_pos è la posizione del secondo template una volta allungato
                             #print("tmp_pos: "+str(tmp_pos))

                         else:
                             MH1_2 = self.ref_genome_up[(- MH_lengths -i +1):]
                             print("elongations finished")
                             print("MH1_2: "+str(MH1_2))
                             MH1_2_indel_P1 = MH1_2 + self.INDEL + P1
                             print("MH1_2_indel_P1: " + MH1_2_indel_P1)
                             tmp_pos_2 = self.ref_genome_down[
                                     len(self.ANC) + len(P1):].index(MH1_2_indel_P1)
                             #print("tmp_pos_2: " +str(tmp_pos_2))
                             break
                         #if tmp_pos < (len(self.INDEL) + len(P1_2)): break
                         #why i sholud break the loop?
   
                    for i in range(1, rep_pat_pos_2):
                        P1_1 = self.ref_genome_down[
                                len(self.ANC) : len(self.ANC) + MH_lengths+i]
                        print("P1_1:"+P1_1)
                        MH1_2_indel_P1_1 = MH1_2 + self.INDEL + P1_1
                        print(MH1_2_indel_P1_1)

                        if MH1_2_indel_P1_1 in self.ref_genome_down [ len(self.ANC) + len(P1):]:
                            continue
                        else:
                            P1_2 = P1_1[:-1]
                            print("P1_2: "+P1_2)
                            DS_pattern = MH1_2 + self.INDEL + P1_2
                            print("DS_pattern : "+ DS_pattern)
                            break
                    else:
                        
                        DS_pattern=MH1_2_indel_P1 #se P2 non viene allungato allora il BR rimane con R non allungato
                    
                    rep_pat = DS_pattern
                    rep_pat_pos_2 = self.mutant_sequence[
                          (self.indel_position+len(self.INDEL)+len(P1_2)):].index(rep_pat)
                    print("rep_pat_pos_2: "+str (rep_pat_pos_2)) #distanza fra P1 e MH2

                    print("DS_pattern : "+ DS_pattern)
                    mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(self.DER))+ len(P1_2) : (len(self.DER))+ len(P1_2) + rep_pat_pos_2]}'
                    print("Y seq :"+mmej_marked_inter_reps_seq) 
                    
                    # set output variables as attributes
                    SD_D_Substitution_mutant_pattern = self.direct_substitution_mutant_pattern_generator(P2=P1_2, MH2=MH1_2,rep_pat=rep_pat)
                    print (SD_D_Substitution_mutant_pattern)
                                 
                    SD_DS_P2 = P1_2
                    SD_DS_MH2 = MH1_2
                    SD_DS_repeat_pat = rep_pat
                    SD_DS_repeat_pat_len = len (rep_pat)
                    SD_DS_between_reps = len(mmej_marked_inter_reps_seq)
                    SD_DS_last_dimer = SD_DS_repeat_pat[-2:]
            
                    temp_SD_direct_substitution=get_motifs_freqs(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, large_window=1000,small_window=self.windowsize,motif=rep_pat,indel_type='DEL')
                    print (temp_SD_direct_substitution)
            
                    SD_DS_motif_pos.append(temp_SD_direct_substitution[0])
                    SD_DS_motif_freq_small.append(temp_SD_direct_substitution[1])
                    SD_DS_motif_freq_large.append(temp_SD_direct_substitution[2])
            
       
                    #setting output dict
                    _d = {
                        'SD_direct_substitution': SD_direct_substitution, 
                        'SD_DS_mutant_pattern':SD_D_Substitution_mutant_pattern,
                        'SD_DS_P2': SD_DS_P2,'SD_DS_MH2': SD_DS_MH2, 
                        'SD_DS_repeat_pat': SD_DS_repeat_pat,
                        'SD_DS_repeat_pat_len': SD_DS_repeat_pat_len,
                        'SD_DS_dist_between_reps': SD_DS_between_reps,
                        'SD_DS_last_dimer': SD_DS_last_dimer,
                        'SD_DS_motif_pos': SD_DS_motif_pos,
                        'SD_DS_motif_freq_small': SD_DS_motif_freq_small,
                        'SD_DS_motif_freq_large': SD_DS_motif_freq_large
                         }
                    self.loop_out_dict = _d 

    def sd_inverted_substitution(self):
        print("self.ANC:"+self.ANC)
        print("self.DER:"+self.DER)
        """
        """
        
        SD_IS_motif_pos, SD_IS_motif_freq_large, SD_IS_motif_freq_small = list(), list(), list()               
        ###1st step is to bring MH1 and P1 from the ref seq 
        for MH_lengths in self.MH_lengths:
            MH1 = self.ref_genome_up [ - MH_lengths :]  
            P1 = self.ref_genome_down [ len(self.ANC) : len(self.ANC) + MH_lengths ]
            #create the IS_pattern and search it
            IS_pattern = self.reverse_complement_converter(seq = MH1 + self.INDEL + P1)
            print("IS_pattern: " + IS_pattern)
            print("self.INDEL :" + self.INDEL)
            print("MH1: "+ MH1)
            print("P1: "+ P1)
            #print("ex: " + str(self.extension))
            #search the pattern after 1° P
            if IS_pattern in self.ref_genome_down [len(self.ANC) + MH_lengths:]:
                SD_inverted_substitution = True
            else:
                SD_inverted_substitution = False 
                print("RZB not found")

            if SD_inverted_substitution:
                print("RZB found")
               
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
                
                    print("MHrevc: " + MHrevc)
                    print("Prevc: " + Prevc)

                    mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(self.DER))+ len(P1) : (len(self.DER))+ len(P1) + rep_pat_pos]}'
                    print("Y seq :"+mmej_marked_inter_reps_seq)
        
                    
                    # set output variables as attributes
                    SD_I_Substitution_mutant_pattern = self.inverted_substitution_mutant_pattern_generator(P1=P1, MH1=MH1,rep_pat=rep_pat)
                    print (SD_I_Substitution_mutant_pattern)
                                 
                    SD_IS_Prevc = Prevc
                    SD_IS_MHrevc = MHrevc
                    SD_IS_repeat_pat = rep_pat
                    SD_IS_repeat_pat_len = len (rep_pat)
                    SD_IS_between_reps = len(mmej_marked_inter_reps_seq)
                    SD_IS_last_dimer = SD_IS_repeat_pat[-2:]
            
                    temp_SD_inverted_substitution=get_motifs_freqs(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, large_window=1000,small_window=self.windowsize,motif=rep_pat,indel_type='DEL')
                    print (temp_SD_inverted_substitution)
            
                    SD_IS_motif_pos.append(temp_SD_inverted_substitution[0])
                    SD_IS_motif_freq_small.append(temp_SD_inverted_substitution[1])
                    SD_IS_motif_freq_large.append(temp_SD_inverted_substitution[2])
            
       
                    #setting output dict
                    _d = {
                        'SD_inverted_substitution': SD_inverted_substitution, 
                        'SD_IS_mutant_pattern':SD_I_Substitution_mutant_pattern,
                        'SD_IS_Prevc': SD_IS_Prevc,'SD_DS_MHrevc': SD_IS_MHrevc, 
                        'SD_IS_repeat_pat': SD_IS_repeat_pat,
                        'SD_IS_repeat_pat_len': SD_IS_repeat_pat_len,
                        'SD_IS_dist_between_reps': SD_IS_between_reps,
                        'SD_IS_last_dimer': SD_IS_last_dimer,
                        'SD_IS_motif_pos': SD_IS_motif_pos,
                        'SD_IS_motif_freq_small': SD_IS_motif_freq_small,
                        'SD_IS_motif_freq_large': SD_IS_motif_freq_large
                         }
                    self.snap_out_dict = _d 

                if self.extension:
                    print("ON")
                    rep_pat = IS_pattern
                    rep_pat_pos_2 = self.mutant_sequence[
                            (self.indel_position):].index(rep_pat)
                    #distance from pos to Prevc
                      

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
                             #print("tmp_pos_2: " +str(tmp_pos_2))
                             break
                    for i in range(1, rep_pat_pos_2):
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

                    rep_pat = IS_pattern
                    rep_pat_pos_2 = self.mutant_sequence[
                          (self.indel_position+len(self.INDEL)+len(P1_2)):].index(rep_pat)
                    print("rep_pat_pos_2: "+str (rep_pat_pos_2)) #distanza fra P1 e MH2

                    print("IS_pattern : "+ IS_pattern)
                    
                    mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(self.DER))+ len(P1_2) : (len(self.DER))+ len(P1_2) + rep_pat_pos_2]}'
                    print("Y seq :"+mmej_marked_inter_reps_seq)
        
                    
                    # set output variables as attributes
                    SD_I_Substitution_mutant_pattern = self.inverted_substitution_mutant_pattern_generator(P1=P1_2, MH1=MH1_2,rep_pat=rep_pat)
                    print (SD_I_Substitution_mutant_pattern)
                                 
                    SD_IS_Prevc = self.reverse_complement_converter(seq = P1_2) 
                    SD_IS_MHrevc = self.reverse_complement_converter(seq = MH1_2)
                    SD_IS_repeat_pat = rep_pat
                    SD_IS_repeat_pat_len = len (rep_pat)
                    SD_IS_between_reps = len(mmej_marked_inter_reps_seq)
                    SD_IS_last_dimer = SD_IS_repeat_pat[-2:]
            
                    temp_SD_inverted_substitution=get_motifs_freqs(ref=self.refFA, CHR=self.chrom, POS=self.pos_on_chr, large_window=1000,small_window=self.windowsize,motif=rep_pat,indel_type='DEL')
                    print (temp_SD_inverted_substitution)
            
                    SD_IS_motif_pos.append(temp_SD_inverted_substitution[0])
                    SD_IS_motif_freq_small.append(temp_SD_inverted_substitution[1])
                    SD_IS_motif_freq_large.append(temp_SD_inverted_substitution[2])
            
       
                    #setting output dict
                    _d = {
                        'SD_inverted_substitution': SD_inverted_substitution, 
                        'SD_IS_mutant_pattern':SD_I_Substitution_mutant_pattern,
                        'SD_IS_Prevc': SD_IS_Prevc,'SD_IS_MHrevc': SD_IS_MHrevc, 
                        'SD_IS_repeat_pat': SD_IS_repeat_pat,
                        'SD_IS_repeat_pat_len': SD_IS_repeat_pat_len,
                        'SD_IS_dist_between_reps': SD_IS_between_reps,
                        'SD_IS_last_dimer': SD_IS_last_dimer,
                        'SD_IS_motif_pos': SD_IS_motif_pos,
                        'SD_IS_motif_freq_small': SD_IS_motif_freq_small,
                        'SD_IS_motif_freq_large': SD_IS_motif_freq_large
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
            'del_mmej', 'del_mmejl','del_mmej_cand', 'del_mmej_marked' ,'del_mmej_marked_on_ref',
            'del_last_dimer','del_mmej_cand_len','del_mmej_motif_pos','del_mmej_freq_small',
            'del_mmej_freq_large',
            
            # SD_inverted_insertion
            'SD_inverted_insertion','SD_II_mutant_pattern','SD_II_MHrevc','SD_II_Prevc',
            'SD_II_repeat_pat',#'SD_II_repeat_pat_len', 
            'SD_II_last_dimer','SD_II_dist_between_reps','SD_II_motif_pos','SD_II_motif_freq_small',
            'SD_II_motif_freq_large',

            #SD-direct_deletion
            'SD_direct_deletion', 'SD_DD_mutant_pattern','SD_DD_MH2','SD_DD_P2', 
            'SD_DD_repeat_pat','SD_DD_repeat_pat_len','SD_DD_last_dimer',
            'SD_DD_dist_between_reps','SD_DD_motif_pos','SD_DD_motif_freq_small',
            'SD_DD_motif_freq_large',
            
            # SD_direct_insertion
            'SD_direct_insertion','SD_DI_mutant_pattern','SD_DI_MH2','SD_DI_P2',
            'SD_DI_repeat_pat','SD_DI_dist_between_reps','SD_DI_last_dimer',
            'SD_DI_motif_pos','SD_DI_motif_freq_small','SD_DI_motif_freq_large', 
            
            # SD_inverted_deletion
            'SD_inverted_deletion','SD_ID_mutant_pattern','SD_ID_MHrevc','SD_ID_Prevc',
            'SD_ID_repeat_pat','SD_ID_repeat_pat_len','SD_ID_last_dimer','SD_ID_dist_between_reps',
            'SD_ID_motif_pos','SD_ID_motif_freq_small','SD_ID_motif_freq_large',
            #SD_direct_substitution
            'SD_direct_substitution','SD_DS_mutant_pattern','SD_DS_P2','SD_DS_MH2',
            'SD_DS_repeat_pat','SD_DS_repeat_pat_len','SD_DS_dist_between_reps','SD_DS_last_dimer',
            'SD_DS_motif_pos','SD_DS_motif_freq_small','SD_DS_motif_freq_large',
            #SD_inverted_substitution
             'SD_inverted_substitution','SD_IS_mutant_pattern','SD_IS_Prevc','SD_IS_MHrevc',
             'SD_IS_repeat_pat','SD_IS_repeat_pat_len','SD_IS_dist_between_reps','SD_IS_last_dimer',
             'SD_IS_motif_pos','SD_IS_motif_freq_small','SD_IS_motif_freq_large'                         
            # polymerase slippage
            #'pol_slip' , 'pol_slip_submotif' , 'pol_slippage_repeatsIndel',
            #'pol_slippage_repeatsDownstream', 'pos_slip_pos'
            ]

        _Dtypes={
            'del_mmej':str, #'del_mmej_cand':object, #'del_mmej_marked':object ,'del_mmej_marked_on_ref':object,
#            'del_last_dimer':object,'del_mmej_cand_len':object,
            
            #SD_inverted_insertion
            'SD_inverted_insertion':str,'SD_II_mutant_pattern':str,
            'SD_II_MHrevc':str, 'SD_II_Prevc':str,
            'SD_II_repeat_pat':str,'SD_II_last_dimer':str,'SD_II_dist_between_reps':str,
            'SD_II_motif_pos':str,'SD_II_motif_freq_small':str,'SD_II_motif_freq_large':str,
            
            #SD-direct_deletion
            'SD_direct_deletion':str,'SD_DD_mutant_pattern':str,
            'SD_DD_MH2':str,'SD_DD_P2':str,
            'SD_DD_repeat_pat':str,'SD_DD_repeat_pat_len':str,'SD_DD_last_dimer':str,
            'SD_DD_dist_between_reps':str,
            'SD_DD_motif_pos':str, 'SD_DD_motif_freq_small':str, 'SD_DD_motif_freq_large':str,  

            # SD_direct_insertion
            'SD_direct_insertion':str,'SD_DI_mutant_pattern':str,
            'SD_DI_MH2':str,'SD_DI_P2':str,
            'SD_DI_repeat_pat':str,'SD_DI_dist_between_reps':str, 'SD_DI_last_dimer':str,
            'SD_DI_motif_pos':str, 'SD_DI_motif_freq_small':str, 'SD_DI_motif_freq_large':str,
            
            # SD_inverted_deletion
            'SD_inverted_deletion':str,'SD_ID_mutant_pattern':str,
            'SD_ID_MHrevc':str,'SD_ID_Prevc':str,
            'SD_ID_repeat_pat':str,'SD_ID_repeat_pat_len':str,'SD_ID_last_dimer':str,
            'SD_ID_dist_between_reps':str, 'SD_ID_motif_pos':str, 'SD_ID_motif_freq_small':str,
            'SD_ID_motif_freq_large':str,
            #SD_direct_substitution
            'SD_direct_substitution':str,'SD_DS_mutant_pattern':str,
            'SD_DS_MH2':str,'SD_DS_P2':str, 
            'SD_DS_repeat_pat':str, 'SD_DS_repeat_pat_len':str,'SD_DS_last_dimer':str,
            'SD_DS_motif_pos':str, 'SD_DS_motif_freq_small':str, 'SD_DS_motif_freq_large':str,
            #SD_inverted_substitution
            'SD_inverted_substitution':str,'SD_IS_mutant_pattern':str,
            'SD_IS_Prevc':str,'SD_IS_MHrevc':str,
            'SD_IS_repeat_pat':str,'SD_IS_repeat_pat_len':str,'SD_IS_dist_between_reps':str,'SD_IS_last_dimer':str,
            'SD_IS_motif_pos':str,'SD_IS_motif_freq_small':str,'SD_IS_motif_freq_large':str
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

        # setting the right values to match corresponding columns
        ex_df.loc[0,[key for key in d.keys()]] = np.array([val for val in d.values()], dtype=object)
        return ex_df
        
    # The function below are functions the generate a representation
    # of the different patterns per repair mechanism

    # Deletion MMEJ
    def del_mutant_mmej_marked(self, mmej_cand: str):
        """
        Create string that represent the mutated genome after MMEJ. 
        The pattern of an MMEJ that results in a deletion is:
        upstream seq -> MMEJ -> DSB -> DSB to MMEJ -> downstream seq
        Args:
            mmej_cand (str): the MMEJ pattern
        Returns:
            The pattern as a string: upstream_seq (10 bp) * MH * downstream_seq (10bp)
        """
        if self.flip:
            return f'{self.DSB_down[0:10][::-1]}*[{mmej_cand[::-1]}]*{self.DSB_up[(len(self.DSB_up)-len(mmej_cand)-10):(len(self.DSB_up)-len(mmej_cand))][::-1]}'
        return f'{self.DSB_up[(len(self.DSB_up)-len(mmej_cand)-10):(len(self.DSB_up)-len(mmej_cand))]}*[{mmej_cand}]*{self.DSB_down[0:10]}'

    def del_mmej_marked_on_reference(self, mmej_cand: str,
                                       mmej_cand_position: int ):
        """
        Create string that represent the reference genome before MMEJ:
        upstream seq -> MMEJ -> DSB -> DSB to MMEJ -> MMEJ -> downstream seq 
        Args:
            mmej_cand (str): thr MMEJ pattern
            mmej_cand_position (int): the position of the MMEJ seq 
                relative to the DSB
        Returns:
            The pattern as a string: upstream_seq (10 bp) * MH *|* downstream_seq (10bp)
        """
        before_DSB = f'{self.ref_genome_up[(len(self.ref_genome_up)-len(mmej_cand) -10):(len(self.ref_genome_up)-len(mmej_cand))]}*[{mmej_cand}]*|'
        DSB_to_MMEJ_seq = f'{self.ref_genome_down[0:(mmej_cand_position-len(mmej_cand)-1)]}'
        # print(mmej_cand_position, mmej_cand) 
        MMEJ_till_end = self.ref_genome_down[(mmej_cand_position-1):(mmej_cand_position+10)]
        if self.flip:
            before_DSB = f'*[{mmej_cand[::-1]}]*{self.ref_genome_up[(len(self.ref_genome_up)-len(mmej_cand) -10):(len(self.ref_genome_up)-len(mmej_cand))][::-1]}'
            return f'{MMEJ_till_end[::-1]}|*[{mmej_cand[::-1]}]*{DSB_to_MMEJ_seq[::-1]}{before_DSB}'
        return f'{before_DSB}{DSB_to_MMEJ_seq}*[{mmej_cand}]*{MMEJ_till_end}'


    # Insertion (SD-Snap back) MMEJ
    def snap_mutant_pattern_generator(self, P1:str, P2: str, MH1:str, MH2: str,
                                 rep_pat:str):
        """
        Create the pattern that should match the reference sequence in
        cases where SD-Snap back is possible, the pattern is:
        [MH2->INS->P2] -> random seq -> invert_and_complement([MH2->INS->P2])
        Args:
            P1, P2, MH1, MH2 (str): parts of the pattern as shown in fig.1A from:
            Khodaverdian, V. Y. et al. Secondary structure forming sequences drive SD-MMEJ 
            repair of DNA double-strand breaks. Nucleic Acids Res. 45, 12848–12861 (2017).
            rep_pat (str): the following pattern: [MH2->INS->P2]
        Returns:
            The pattern itself in format 10bpupstream*MH2*INS P2 randomseq*P1 INVINS MH1 10bpdownstream
        """
        
        inv_ins = self.reverse_complement_converter(seq = self.INDEL)
        rep_pat_pos = self.mutant_sequence[(self.indel_position+len(P2)):].index(rep_pat)
        mmej_marked_rep_up = f'*MH1[{MH1}]|*INS[{self.INDEL}]P1[{P1}]'
        mmej_marked_up = f'{self.mutant_sequence[(self.indel_position-len(MH2)-10):(self.indel_position-len(MH2))]}{mmej_marked_rep_up}'
        mmej_marked_rep_down = f'*P2[{P2}][{inv_ins}]MH2[{MH2}]{self.DSB_down[(rep_pat_pos+len(rep_pat)+len(P2)):(rep_pat_pos+len(rep_pat)+len(P2)+10)]}'
        mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P2) + self.indel_length):(len(P2) + rep_pat_pos)]}'

        return f'{mmej_marked_up}{mmej_marked_inter_reps_seq}{mmej_marked_rep_down}'

    # Insertion (SD-Loop out) MMEJ
    def loop_mutant_pattern_generator(self, P2: str, MH2: str,
                                 rep_pat:str):
        """
        Create the pattern that should match the reference sequence in
        cases where SD-Snap back is possible, the pattern is:
        [MH2->INS->P2] -> random seq -> [MH2->INS->P2]
        Args:
            P2, MH2 (str):parts of the pattern as shown in fig.1A from:
            Khodaverdian, V. Y. et al. Secondary structure forming sequences drive SD-MMEJ 
            repair of DNA double-strand breaks. Nucleic Acids Res. 45, 12848–12861 (2017).
            rep_pat (str): the following pattern: [MH2->INS->P2]
        Returns:
            The mutant genomic context.
        """
        
        rep_pat_pos = self.mutant_sequence[(self.indel_position+len(P2)):].index(rep_pat)
        mmej_marked_rep_up = f'*MH1[{MH2}]|*INS[{self.INDEL}]P1[{P2}]'
        mmej_marked_up = f'{self.mutant_sequence[(self.indel_position-len(MH2)-10):(self.indel_position-len(MH2))]}{mmej_marked_rep_up}'
        mmej_marked_rep_down = f'*MH2[{MH2}][{self.INDEL}]P2[{P2}]{self.DSB_down[(rep_pat_pos+len(rep_pat)+len(P2)):(rep_pat_pos+len(rep_pat)+len(P2)+10)]}'
        mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P2) + self.indel_length):(len(P2) + rep_pat_pos)]}'

        return f'{mmej_marked_up}{mmej_marked_inter_reps_seq}{mmej_marked_rep_down}'
    
    # fabio:for sd_direct_deletions
    def direct_deletions_mutant_pattern_generator(self, P2: str, MH2: str,
                                 rep_pat:str):
        """
        Create the MUTATED pattern that should occurs when a SD_direct_deletion acts, match the reference sequence in
        cases where SD-Snap back is possible, the pattern is:
        [MH1-DELETED-P2] -> random seq -> [MH2->P2]
        Args:
            P2, MH2 (str):parts of the pattern as shown in fabio's drawings 
            Returns:
            The mutant genomic context.
        """
 
        rep_pat_pos_2 = self.mutant_sequence[
                    (self.indel_position + len(self.INDEL) + len (P2)):].index(rep_pat)
        
        rep_pat_pos = self.mutant_sequence[(self.indel_position+len(P2)):].index(rep_pat)
        mmej_marked_rep_up = f'*MH1[{MH2}]|[DELETION]|P1[{P2}]'
        mmej_marked_up = f'{self.mutant_sequence[(self.indel_position-len(MH2)-10):(self.indel_position-len(MH2))]}{mmej_marked_rep_up}'
        mmej_marked_rep_down = f'*MH2[{MH2}]|P2[{P2}]{self.DSB_down[(rep_pat_pos+len(rep_pat)+len(P2)):(rep_pat_pos+len(rep_pat)+len(P2)+10)]}'
        #mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P2) + self.indel_length):(len(P2) + rep_pat_pos)]}'

        mmej_marked_inter_reps_seq_2 = f'{self.DSB_down[len(P2) + self.indel_length:(len(P2) + self.indel_length + rep_pat_pos_2)]}'
       
        return f'{mmej_marked_up}{mmej_marked_inter_reps_seq_2}{mmej_marked_rep_down}'
    def inverted_deletions_mutant_pattern_generator(self, P1: str, MH1: str,
                                 rep_pat:str):
        """
        Create the MUTATED pattern that should occurs when a SD_inverted_deletion acts, match the reference sequence in
        cases where SD-Snap back is possible, the pattern is:
        [MH1-DELETED-P1] -> random seq -> [MH2->P2] (MH2->P2 are the reverse complement of MH1 and P1)
        Args:
            P2, MH2 (str):parts of the pattern as shown in fabio's drawings 
            rep_pat: the pattern founded = [MH2->P2]
            Returns:
            The mutant genomic context.
        """
 
        P2revc = self.reverse_complement_converter(seq = P1)
        MH2revc = self.reverse_complement_converter(seq = MH1)
        rep_pat_pos_2 = self.mutant_sequence[
                    (self.indel_position + len(self.INDEL) + len (P1)):].index(rep_pat)
        
        rep_pat_pos = self.mutant_sequence[(self.indel_position+len(P1)):].index(rep_pat)
        mmej_marked_rep_up = f'*MH1[{MH1}]|[DELETION]|P1[{P1}]'
        mmej_marked_up = f'{self.mutant_sequence[(self.indel_position-len(MH1)-10):(self.indel_position-len(MH1))]}{mmej_marked_rep_up}'
        mmej_marked_rep_down = f'*P2revc[{P2revc}]|MH2revc[{MH2revc}]{self.DSB_down[(rep_pat_pos+len(rep_pat)+len(P1)):(rep_pat_pos+len(rep_pat)+len(P1)+10)]}'
        #mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P2) + self.indel_length):(len(P2) + rep_pat_pos)]}'

        mmej_marked_inter_reps_seq_2 = f'{self.DSB_down[len(P1) + self.indel_length:(len(P1) + self.indel_length + rep_pat_pos_2)]}'
       
        return f'{mmej_marked_up}{mmej_marked_inter_reps_seq_2}{mmej_marked_rep_down}'

    def direct_substitution_mutant_pattern_generator(self, P2: str, MH2: str, rep_pat:str):
    
        rep_pat_pos = self.mutant_sequence[(self.indel_position+len(P2)):].index(rep_pat)
        mmej_marked_rep_up = f'*MH1[{MH2}]|*Z_SEQ[{self.INDEL}]P1[{P2}]'
        mmej_marked_up = f'{self.mutant_sequence[(self.indel_position-len(MH2)-10):(self.indel_position-len(MH2))]}{mmej_marked_rep_up}'
        mmej_marked_rep_down = f'*MH2[{MH2}]*Z_SEQ [{self.INDEL}]P2[{P2}]{self.DSB_down[(rep_pat_pos+len(rep_pat)+len(P2)):(rep_pat_pos+len(rep_pat)+len(P2)+10)]}'
        mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(self.DER))+ len(P2) : (len(self.DER))+ len(P2) + rep_pat_pos]}'


        return f'{mmej_marked_up}{mmej_marked_inter_reps_seq}{mmej_marked_rep_down}'

    def inverted_substitution_mutant_pattern_generator(self, P1: str, MH1: str,
            rep_pat:str):

        """
        Create the MUTATED pattern that should occurs when a SD_inverted_substitution acts, match the reference sequence in
        cases where SD-Snap back is possible, the pattern is:
        [MH1-Z-REVC-P1] -> random seq -> [MH2->Z_SEQ->P2] (MH2->P2 are the reverse complement of MH1 and P1)
        Args:
            P2, MH2 (str):parts of the pattern as shown in fabio's drawings
            rep_pat: the pattern founded = [MH2->Zrevc->P2]
            Returns:
            The mutant genomic context.
        """

        P2revc = self.reverse_complement_converter(seq = P1)
        MH2revc = self.reverse_complement_converter(seq = MH1)
        rep_pat_pos_2 = self.mutant_sequence[
                (self.indel_position + len(self.INDEL) + len (P1)):].index(rep_pat)
        Z_seq = self.reverse_complement_converter(seq = self.INDEL)
        
        rep_pat_pos = self.mutant_sequence[(self.indel_position+len(P1)):].index(rep_pat)
        mmej_marked_rep_up = f'*MH1[{MH1}]|*Z_REVC[{self.INDEL}]P1[{P1}]'
        mmej_marked_up = f'{self.mutant_sequence[(self.indel_position-len(MH1)-10):(self.indel_position-len(MH1))]}{mmej_marked_rep_up}'
        mmej_marked_rep_down = f'*MH2revc[{MH2revc}]Z_SEQ[{Z_seq}]P2revc[{P2revc}]{self.DSB_down[(rep_pat_pos+len(rep_pat)+len(P1)):(rep_pat_pos+len(rep_pat)+len(P1)+10)]}'
        mmej_marked_inter_reps_seq_2 = f'{self.DSB_down[len(P1) + self.indel_length:(len(P1) + self.indel_length + rep_pat_pos_2)]}'

        return f'{mmej_marked_up}{mmej_marked_inter_reps_seq_2}{mmej_marked_rep_down}'
