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

class emMEJrealignment:
    def __init__(self, 
                ANC: str, DER: str,
                indel_type: str, flip: bool, include_context: bool,
                windowsize: int, refFA: str, MH_len_early_stop: int, 
                chrom: str, pos_on_chr:int):
        self.indel_position = windowsize #fabri: OK, this is weird but try to understand that.
        self.indel_type = indel_type 
        self.flip = flip        
        self.include_context = include_context
        if MH_len_early_stop != 0:
            self.MH_early_stop = int(MH_len_early_stop)
        else: self.MH_early_stop = "None"
        if (self.indel_type == 'INS'): self.indel_length = len(DER) - 1
        if (self.indel_type == 'DEL'): self.indel_length = len(ANC) - 1
        
        # making all VCF format set to the right standart
        if (self.indel_type == 'INS'):
            # cropping indel sequence to the right size. fabri: is this for substitutions? Then why before length-1?
            self.INDEL =  DER[len(ANC):]
            # getting a -30->indel->+30 reference sequence using get_ref_context, defined in src/pysam_getfasta.py. Returns sequence.
            self.ref_seq = get_ref_context(refFA=refFA, chrom=chrom, 
                            indel_pos=pos_on_chr,
                            context_window_size=windowsize,
                            indel_seq=self.INDEL)
            self.ref_seq = self.ref_seq.upper()
            # generating the mutant context
            self.mutant_sequence = mutant_context_generator(
                        reference_contex=self.ref_seq,
                        alt=DER,
                        ref=ANC,
                        window_size=(windowsize))

            self.ANC = ANC
            self.INDEL =  DER[len(ANC):] #fabri: duplicated line?
            # Flipping sequences (looking upward, i.e. in 5' direction from indel, without complement)
            if self.flip:
                self.ref_seq = self.ref_seq[::-1]
                self.mutant_sequence = self.mutant_sequence[::-1]
                self.ANC = self.ANC[::-1]
                self.INDEL = self.INDEL[::-1]
                self.indel_position = self.indel_position + self.indel_length +1
            
        if (self.indel_type == 'DEL'):
            windowsize = windowsize + self.indel_length
            self.indel_position = windowsize
            # cropping indel sequence to the right size
            self.ANC = ANC[len(DER):]
            self.INDEL = DER[len(DER):]
            # self.indel_length = len(self.ANC)
            
            # Take care of the cases of fliped sequences
            if self.flip:
                self.ANC = self.ANC[::-1]
                self.indel_position = self.indel_position -self.indel_length +1
            # getting a -30->indel->+30 refenrence context
            self.ref_seq = get_ref_context(
                        refFA=refFA, chrom=chrom,
                        indel_pos=pos_on_chr,
                        context_window_size=(windowsize),
                         indel_seq=self.INDEL)
            self.ref_seq = self.ref_seq.upper()
            
            if self.flip: self.ref_seq = self.ref_seq[::-1]
            # generating the mutant context #fabri: the mutant context is the alternative?
            self.mutant_sequence = mutant_context_generator(
                        reference_contex=self.ref_seq,
                        alt=self.INDEL, 
                        ref=self.ANC, 
                        window_size=windowsize) 
                           
        ### Inputs Dtype chacking ###    
        # assert check_argument_types() 
        # type_checker(self.__dict__, self.__init__.__annotations__.values(),mode='class') 

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

    def seq_indel_spliter(self):
        """
        seq_indel_spliter will get self as an input and will return
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
        if self.MH_early_stop != "None":
            mmej_motifs = mmej_motifs[mmej_motifs.str.len() == self.MH_early_stop]
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

    def microhomology_detection(self):
        self.seq_indel_spliter() # get seqs upstream and downstream

        # initializing general attributes in order to create a complete
        # attribute set for all indels
        self.microhomology = False
        self.snap_out_dict = {}
        self.loop_out_dict = {}
        self.pol_slip_dict = {}
        self.del_out_dict = {}

        """
        Since the DSB in snap-back and loop-out ist ocur in between the 
        two MH, there is no need to check if suck MH exist on both sides 
        of the breake (as with trans-MMEJ and deletions).
        Therefore we call the snap-back and loop-out functions right away if
        an indel is an insertion
        """
        if self.indel_type == 'INS': 
            self.sd_snap_back_MMEJ()
            self.sd_loop_out_MMEJ()
            self.pol_slip()

        if self.indel_type == 'DEL': 
            self.deletion_MMEJ()


    def deletion_MMEJ(self):
        """
        Finds the longest subsequence (a motif) that exist both before and after
        the DSB, and ends at the last base of the deletion sequence
        """
        # set default values for mmej_cand
        del_mmej = False
        del_mmej_cand, del_mmej_cand_len ,del_mmej_marked= np.nan, np.nan, np.nan
        del_last_dimer, mmej_cand = np.nan ,np.nan
        distance_to_indel_position, del_mmej_marked_on_ref = np.nan, np.nan
        # 1st condition:
        self.first_MH_finder() # finding all possible MH motifs
        self.first_MH_indicator()

        # 2nd condition
        # reverse integrating over all motifs and look whether it exists before the break. If not breaks and longest is considered as the potential MH
        for mh_motif in self.mmej_motifs[::-1]:
            if mh_motif == self.ref_genome_up[-len(mh_motif):]:
                del_mmej = True
                mmej_cand = del_mmej_cand = mh_motif
                mmej_cand_position = self.indel_length + 1 #fabri: why is this within the loop? should be above for efficiency
                #also, ugly and misleading name for variable. This is better as simply indel_length+1, or indel_length_plus1
                self.distance_to_indel_position = self.indel_length - len(mmej_cand)
                break
            else: continue
        
        # get sequences with marked components of the mechanism pattern 
        if del_mmej:
            del_mmej_marked = self.del_mutant_mmej_marked(mmej_cand)
            del_mmej_marked_on_ref = self.del_mmej_marked_on_reference(mmej_cand,
                                       mmej_cand_position)
        #define last dimer, because used later for 2nd order Markov Chain                
            if len(mmej_cand) >= 2: del_last_dimer = mmej_cand[-2:]
            else: del_last_dimer = self.ref_genome_up[-2] + mmej_cand
            del_last_dimer = del_last_dimer.upper()
            del_mmej_cand_len = int(len(del_mmej_cand))

        else:
            distance_to_indel_position = np.nan
            del_mmej_marked_on_ref = np.nan
            del_mmej = False
            del_mmej_cand = mmej_cand
            self.distance_to_indel_position = distance_to_indel_position
        
        _d = {
            'del_mmej': del_mmej,'del_mmej_cand': del_mmej_cand,
            'del_mmej_marked': del_mmej_marked ,'del_mmej_marked_on_ref':del_mmej_marked_on_ref,
            'del_last_dimer':del_last_dimer, 'del_mmej_cand_len':del_mmej_cand_len
            }
        
        self.del_out_dict = _d

            
    def pol_slip(self):
        """
        A function that looks for the pattern of a polymerase slippage, i.e. a submotif preceding the indel that is repeated tp form the indel

        """
        submotif_befor_DSB="";submotif_after_DSB=""; #added to avoid bug in Guy's version
        for i in range(1,len(self.INDEL)):
            submotif_befor_DSB = self.DSB_up[-i:] #DSB_up (see up in glossary): context right before beginning of indel (-i bases).
            submotif_after_DSB = self.INDEL[-i:]
            if submotif_befor_DSB == submotif_after_DSB:
                continue
            else:
                submotif_befor_DSB = submotif_befor_DSB[1:]
                submotif_after_DSB = submotif_after_DSB[1:] #1: because the first position/index(0) does not fulfil the condition. So we take the motif after 0 (which is 1:).
                break
        if len(submotif_after_DSB) > 0:
            submotif_pos = np.array([m.end() for m in #submotif_pos is the end positions of the submotif in the insertion sequence (e.g. TA if insertion is TATATA, so indices are 2,4,6 (1based))
                reg.finditer(submotif_befor_DSB, self.INDEL, overlapped=False)])
            #print(submotif_pos)
            if submotif_befor_DSB*len(submotif_pos) == self.INDEL: #checking that by repeating a submotif I get indel right length and indel
                pol_slip = True
                pol_slip_submotif = submotif_befor_DSB
                pol_slippage_times = len(submotif_pos)

                _d = { #pol_slip is flag whether is slippage or not
                    'pol_slip': pol_slip,'pol_slip_submotif': pol_slip_submotif, #pol_slip_motif is the submotif
                    'pol_slippage_times': pol_slippage_times , 
                    'pol_slippage_last_dimer':submotif_befor_DSB[-2:], 'pol_slip_motif_len':len(pol_slip_submotif) #last dimer before insertion - reported here because we use it to start the chain for markov in MMEJ calcuations. Here useless.
                    } #length of submotif
                
                self.pol_slip_dict = _d
        #print(self.pol_slip_dict)

    
    def sd_snap_back_MMEJ(self):
        """
        sd_loop_out_MMEJ will detect SD-Snap back MMEJ by looking for the following pattern:
        [MH2->INS->P2] -> random seq (with length>0) -> invert_and_complement([MH2->INS->P2])
        This pattern is based on the one represented in fig.1B from:
        Khodaverdian, V. Y. et al. Secondary structure forming sequences drive SD-MMEJ 
        repair of DNA double-strand breaks. Nucleic Acids Res. 45, 12848–12861 (2017).

        Args: 
            self attributes of the EMMEJrealignment object
        Returns:
            snap_mmej_marked (str): the whole pattern of SD-snap-back MMEJ 
                marked on the mutant context
            snap_P1 (str): P1 as represented in Fig1B.
            snap_P2 (str): P2 as represented in Fig1B.
            snap_mh1 (str): MH1 as represented in Fig1B.
            snap_mh2 (str): MH2 as represented in Fig1B.
            snap_inv_comp_repeat_pat (str):The 1st repeat as represented in Fig1B
                (the one serves as a template for thr elongation of the loop).
            snap_repeat_pat (str): The 2nd repeat as represented in Fig1B
                (the one that get sythesized by elongation of the loop).
        """

        #1st step is to construct an initial seed pattern, 
        #which is the minimal pattern to define a snap-back,
        #called minimal_pattern, as follows:
        #1bp -> inv_comp(DER) -> 1bp
        
        mh2 = self.DSB_up[(-1):]
        mh1 = self.reverse_complement_converter(seq = mh2)
        P2 = self.DSB_down[(self.indel_length):(self.indel_length+1)] 
        P1 = self.reverse_complement_converter(seq = P2)
        
        minimal_pattern = f'{mh2}{self.INDEL}{P2}'
        inv_comp_minimal_pattern = self.reverse_complement_converter(seq=minimal_pattern)
        
        # making sure that inv_comp_minimal_pattern is in self.DSB_down
        # before start to find the complite pattern (MH1,MH2,P1 and P2)
        if inv_comp_minimal_pattern in self.DSB_down[len(inv_comp_minimal_pattern):]:
            SD_snap_back = True
        else: SD_snap_back = False
       
        # 2nd step is to find the MHs (5' side of pattern, i.e. MH2 and MH1),
        # which should also be invert-complementary to one another).
        # Search range here to 20bp (longer
        # MH seq would be called SSA which is another mechanism)
        MH_max_range = 20
        if SD_snap_back :
            for i in range(1,MH_max_range):
                mh2 = self.DSB_up[(-i):]
                mh1 = self.reverse_complement_converter(seq = mh2)
                if f'{inv_comp_minimal_pattern[:-1]}{mh1}' in self.DSB_down[len(minimal_pattern): ]: continue
                # making sure that if the only pattern that exist is minimal_pattern,
                # we won't shorten mh1
                elif i==1: break
                else: 
                    # taking the last mh2 and mh1 that matched
                    mh2 = mh2[1:]
                    mh1 = mh1[:-1]
                    break  

        # 3rd step is to find P2 and P1 (should also be invert-complementary 
        # to one another) only if (fabri??)mh1 is found downstream
        # to the insertion end point
        
        if SD_snap_back :
            # start at the first nucleotide that comes right after the insertion,
            # calculate the invert-complementary version and see if it exist downstream
            # to the insertion end point + P2. if true, elongate by 1bp

            for i in range(1,len(self.DSB_down)):
                P2 = self.DSB_down[self.indel_length:(self.indel_length + i)]
                inv_comp_rep_pat = f'{mh2}{self.INDEL}{P2}'
                rep_pat = self.reverse_complement_converter(seq = inv_comp_rep_pat)
                if rep_pat in self.DSB_down[
                    (self.indel_length + len(P2)): ]:
                    continue
                
                # making sure that if the only pattern that exist is minimal_pattern,
                # we won't shorten P2
                elif i==1: break
                else:
                    # taking the longest P2 and P1
                    P2 = P2[:-1]
                    P1 = self.reverse_complement_converter(seq = P2)
                    # taking the longest inv_comp_rep_pat and rep_pat
                    inv_comp_rep_pat = inv_comp_rep_pat[:-1]
                    rep_pat = rep_pat[1:]
                    break
            rep_pat_pos = self.mutant_sequence[(self.indel_position+len(P2)):].index(rep_pat)
            mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P2) + self.indel_length):(len(P2) + rep_pat_pos)]}'
            # set different variables as attributes
            snap_mmej_marked = self.snap_mutant_pattern_generator(P1=P1, P2=P2, MH1=mh1, MH2=mh2,
                                 rep_pat=rep_pat)

            snap_P1, snap_P2 = P1, P2
            snap_mh1, snap_mh2 = mh1, mh2
            snap_repeat_pat, snap_inv_comp_repeat_pat = rep_pat, inv_comp_rep_pat
            dist_between_reps = len(mmej_marked_inter_reps_seq)
            snap_last_dime = inv_comp_rep_pat[-2:]

            _d = {'SD_snap_back': SD_snap_back, 'snap_mmej_marked':snap_mmej_marked, 
                    'snap_P1': snap_P1, 'snap_P2':snap_P2, 'snap_mh1':snap_mh1, 
                'snap_mh2': snap_mh2, 'snap_repeat_pat':snap_repeat_pat,
                'snap_inv_comp_repeat_pat': snap_inv_comp_repeat_pat,
                'snap_last_dimer':snap_last_dime, 'snap_dist_between_reps': dist_between_reps}
            
            
            self.snap_out_dict = _d


    def sd_loop_out_MMEJ(self):
        """
        ####### CORRECTION ###########
        sd_loop_out_MMEJ will detect sd loop out MMEJ by looking for the following pattern:
        [MH2->INS->P2] -> random seq -> [MH2->INS->P2]
        This pattern is based on the one represented in fig.1A from:
        Khodaverdian, V. Y. et al. Secondary structure forming sequences drive SD-MMEJ 
        repair of DNA double-strand breaks. Nucleic Acids Res. 45, 12848–12861 (2017).
        Procedure works by first expanding MH2, then P2. This does not guarantee finding all, nor longest.
        Not a very good approach.

        Args: 
            self attributes of the EMMEJrealignment object

        returns:
            loop_mmej_marked (str): the whole pattern of SD-loop-out MMEJ 
                marked on the mutant context
            loop_P2 (str): P2 as represented in Fig1A.
            loop_mh2 (str):mh2 (the actual syntetic MMEJ) as represented in Fig1A.
            loop_repeat_pat (str): The repeat that contains the insertion (1st rep
                in Fig1A).
        """  

        ###1st step is to construct a small minimal_pattern: 1bp -> DER -> 1bp
        #which is the smallest pattern that can be defined as loop-out 
        MH1 = self.DSB_up[(-1):]
        P2 = self.DSB_down[(self.indel_length):(self.indel_length+1)]
        minimal_pattern = MH1 + self.INDEL + P2

        # making sure that minimal_pattern is in self.DSB_down
        if minimal_pattern in self.DSB_down[(len(self.INDEL) + len(P2) + 1):]:
            SD_loop_out = True
        else: SD_loop_out = False
        
        if SD_loop_out:
            #NB: SD Loop Out demands the formation of a stem-loop structure, 
            #which might not form if the two MH->ins->P patterns are adjacent. 
            #We must ensure that there is at least 1bp between the two. 
            #fabri: is this done, and is it appropriate?
            
            # position of the MH->ins->P 2nd pattern
            minimal_pattern_pos = self.DSB_down[
                (len(self.INDEL) + len(P2)):].index(minimal_pattern)
            # extention_space is the longest sequence that mh2 can be found in. 
            #Note that in snapback this is not necessary because patterns have opposite orientation
            extention_space = self.DSB_down[
                (len(self.INDEL) + len(P2)):
                (len(self.INDEL) + len(P2) +minimal_pattern_pos)]
            
            # Enlarge minimal_pattern by increasing mh2, then finding position of minimal_pattern (tmp_pos) downstream (within DSB_down)
            if len(extention_space) > 1:
                # elongation of mh2, confined in extention_space
                for i in range(1,len(extention_space)):
                    mh2 = self.DSB_up[(-i):]
                    minimal_pattern = mh2 + self.INDEL + P2

                    if minimal_pattern in self.DSB_down[(len(self.INDEL) + len(P2)):]:
                        tmp_pos = self.DSB_down[(len(self.INDEL) + len(P2)):].index(minimal_pattern)
                    else:
                        mh2 = self.DSB_up[(-i+1):]
                        minimal_pattern = mh2 + self.INDEL + P2
                        tmp_pos = self.DSB_down[
                            (len(self.INDEL) + len(P2)):].index(minimal_pattern)
                        break
                    if tmp_pos < (len(self.INDEL) + len(P2)): break

                # elongation of P2, confined in extention_space[:*starting position of mh2*]
                for i in range(1,len(extention_space[:(tmp_pos+1)])):
                    P2 = self.DSB_down[(self.indel_length):(self.indel_length+i)]
                    minimal_pattern = mh2 + self.INDEL + P2

                    if minimal_pattern in self.DSB_down[(len(self.INDEL) + len(P2)):]:
                        continue
                    else:
                        P2 = P2[:-1]
                        minimal_pattern = mh2 + self.INDEL + P2
                        break
            else: mh2 = MH1
            # Crearing mmej_marked
            rep_pat = minimal_pattern
            rep_pat_pos = self.mutant_sequence[
                (self.indel_position+len(P2)):].index(rep_pat)
            mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P2) + self.indel_length):(len(P2) + rep_pat_pos)]}'
            
            # set output variables as attributes
            loop_mmej_marked = self.loop_mutant_pattern_generator(P2=P2, MH2=mh2,
                                 rep_pat=rep_pat)
            loop_P2 = P2
            loop_mh2 = mh2
            loop_repeat_pat = rep_pat
            dist_between_reps = len(mmej_marked_inter_reps_seq)
            loop_last_dimer = loop_repeat_pat[-2:]
            # setting output dict
            _d = {
                'SD_loop_out': SD_loop_out, 'loop_mmej_marked':loop_mmej_marked,
                    'loop_P2': loop_P2,'loop_mh2': loop_mh2, 
                    'loop_repeat_pat': loop_repeat_pat,
                    'loop_dist_between_reps': dist_between_reps,
                    'loop_last_dimer': loop_last_dimer
                 }
            self.loop_out_dict = _d 


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
            'del_mmej', 'del_mmej_cand', 'del_mmej_marked' ,'del_mmej_marked_on_ref',
            'del_last_dimer','del_mmej_cand_len',
            # snap-back
            'SD_snap_back', 'snap_mmej_marked', 'snap_P1', 'snap_P2','snap_mh1', 
            'snap_mh2', 'snap_repeat_pat','snap_inv_comp_repeat_pat', 'snap_last_dimer',
            'snap_dist_between_reps',
            # loop-out
            'SD_loop_out','loop_mmej_marked', 'loop_P2','loop_mh2','loop_repeat_pat',
            'loop_last_dimer', 'loop_dist_between_reps',
            # polymerase slippage
            'pol_slip', 'pol_slip_submotif', 'pol_slippage_times',
            'pol_slippage_last_dimer', 'pol_slip_motif_len']
        
        if self.include_context: _colnames = _colnames + ['ref_genome_context', 'mutant_sequence']

        ex_df = pd.DataFrame(columns=_colnames)
        ex_df.loc[0,:] = np.nan
        # merging dicts from all paths
        d = {**self.del_out_dict, 
            **self.snap_out_dict, **self.loop_out_dict, **self.pol_slip_dict}
        if self.include_context: d = {**d, **self.context}

        # setting the right values to match corresponding columns
        ex_df.loc[0,[key for key in d.keys()]] = [val for val in d.values()]
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
        mmej_marked_rep_up = f'*MH2[{MH2}]|*INS[{self.INDEL}]P2[{P2}]'
        mmej_marked_up = f'{self.mutant_sequence[(self.indel_position-len(MH2)-10):(self.indel_position-len(MH2))]}{mmej_marked_rep_up}'
        mmej_marked_rep_down = f'*P1[{P1}][{inv_ins}]MH1[{MH1}]{self.DSB_down[(rep_pat_pos+len(rep_pat)+len(P2)):(rep_pat_pos+len(rep_pat)+len(P2)+10)]}'
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
            The pattern itself.
        """
        
        rep_pat_pos = self.mutant_sequence[(self.indel_position+len(P2)):].index(rep_pat)
        mmej_marked_rep_up = f'*MH2[{MH2}]|*INS[{self.INDEL}]P2[{P2}]'
        mmej_marked_up = f'{self.mutant_sequence[(self.indel_position-len(MH2)-10):(self.indel_position-len(MH2))]}{mmej_marked_rep_up}'
        mmej_marked_rep_down = f'*MH1[{MH2}][{self.INDEL}]P1[{P2}]{self.DSB_down[(rep_pat_pos+len(rep_pat)+len(P2)):(rep_pat_pos+len(rep_pat)+len(P2)+10)]}'
        mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P2) + self.indel_length):(len(P2) + rep_pat_pos)]}'

        return f'{mmej_marked_up}{mmej_marked_inter_reps_seq}{mmej_marked_rep_down}'

