""" Code Glossary:
    up -> upstream
    down -> downstream
    ref -> referencr sequence
    individual -> referencr sequence with the indel embeded inside
    cand -> candidate
    dsb -> double strand break
    TS -> tamplate_switch (relevant to insertions)
    inv_comp -> invrs and complementary sequence
"""
# importing libraries
import sys
import re
# import typing

import pandas as pd
import numpy as np
from pysam import FastaFile
from helper import type_checker

sys.path.append('/home/labs/alevy/guyta/guy_master_project/src')
from pysam_getfasta import *

class MicroHomology:
    def __init__(self, 
                ancestral_sequence: str, indel_sequence: str,
                indel_type: str, flip: bool, include_context: bool,
                windowsize: int, refFA: str, MH_len_early_stop: int, 
                chr: str, pos_on_chr:int):
        self.indel_position = windowsize
        self.indel_type = indel_type 
        self.flip = flip        
        self.include_context = include_context
        if MH_len_early_stop != "None":
            self.MH_early_stop = int(MH_len_early_stop)
        else: self.MH_early_stop = "None"
        if (self.indel_type == 'INS'): self.indel_length = len(indel_sequence) - 1
        if (self.indel_type == 'DEL'): self.indel_length = len(ancestral_sequence) - 1
        
        # making all VCF format set to the right standart
        if (self.indel_type == 'INS'):
            # cropping indel sequence to the right size
            self.indel_sequence =  indel_sequence[len(ancestral_sequence):]
            # getting a -30->indel->+30 refenrence context
            self.ref_seq = get_ref_context(refFA=refFA, chr=chr, 
                            indel_pos=pos_on_chr,
                            context_window_size=windowsize,
                            indel_seq=self.indel_sequence)
            self.ref_seq = self.ref_seq.upper()
            # generating the accession context
            self.accession_sequence = accession_context_generator(
                        reference_contex=self.ref_seq,
                        alt=indel_sequence,
                        ref=ancestral_sequence,
                        window_size=(windowsize))

            self.ancestral_sequence = ancestral_sequence
            self.indel_sequence =  indel_sequence[len(ancestral_sequence):]
            # Take care of the cases of fliped sequences
            if self.flip:
                self.ref_seq = self.ref_seq[::-1]
                self.accession_sequence = self.accession_sequence[::-1]
                self.ancestral_sequence = self.ancestral_sequence[::-1]
                self.indel_sequence = self.indel_sequence[::-1]
                self.indel_position = self.indel_position + self.indel_length +1
            
        if (self.indel_type == 'DEL'):
            windowsize = windowsize + self.indel_length
            self.indel_position = windowsize
            # cropping indel sequence to the right size
            self.ancestral_sequence = ancestral_sequence[len(indel_sequence):]
            self.indel_sequence = indel_sequence[len(indel_sequence):]
            # self.indel_length = len(self.ancestral_sequence)
            
            # Take care of the cases of fliped sequences
            if self.flip:
                self.ancestral_sequence = self.ancestral_sequence[::-1]
                self.indel_position = self.indel_position -self.indel_length +1
            # getting a -30->indel->+30 refenrence context
            self.ref_seq = get_ref_context(
                        refFA=refFA, chr=chr,
                        indel_pos=pos_on_chr,
                        context_window_size=(windowsize),
                         indel_seq=self.indel_sequence)
            self.ref_seq = self.ref_seq.upper()
            
            if self.flip: self.ref_seq = self.ref_seq[::-1]
            # generating the accession context
            self.accession_sequence = accession_context_generator(
                        reference_contex=self.ref_seq,
                        alt=self.indel_sequence, 
                        ref=self.ancestral_sequence, 
                        window_size=windowsize) 
                           
        ### Inputs Dtype chacking ###    
        # assert check_argument_types() 
        # type_checker(self.__dict__, self.__init__.__annotations__.values(),mode='class') 

        if self.include_context:
            self.context = {
                'ref_genome_context':self.ref_seq, 
                'accession_sequence': self.accession_sequence}

        self.tab = str.maketrans("ACTG", "TGAC")
        self.microhomology = self.microhomology_detection()
        self.ex_data = self.export_data()


    def reverse_complement_converter(self, seq: str):
        """
        A fuction that calculate reverse complements
        """
        return seq.translate(self.tab)[::-1]

    def seq_indel_spliter(self):
        """
        seq_indel_spliter will get self as an input and will return
        the upstream and downstream sequences of both the reference and
        the individual as an output.
        """       
        self.DSB_up = self.accession_sequence[:(self.indel_position)] 
        self.DSB_down = self.accession_sequence[(self.indel_position):] 
        self.ref_genome_up = self.ref_seq[:(self.indel_position)]
        self.ref_genome_down = self.ref_seq[(self.indel_position):]


    def first_MH_finder(self): # change the name because not really a finder
        """ ## For deletion MMEJ: create a list of all possible MH motifs
        Mutate:
            self.mmej_motifs (list): a list of all possible MH motifs    
        """
        mmej_motifs = pd.Series([str(self.ancestral_sequence)[i:] for i in 
                              range((self.indel_length-1), -1, -1)])
        if self.MH_early_stop != "None":
            mmej_motifs = mmej_motifs[mmej_motifs.str.len() == self.MH_early_stop]
        self.mmej_motifs = mmej_motifs.to_list()
        
    def first_MH_indicator(self): # change the name because not really a finder
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
    it will take both reference and the accession sequences
    as an input and will give True/False flage as an output,
    in cases of deletion mmej it will also give mmej_marked_on_ref.
    There are two conditions that must exist in order to detect
    microhomology in deletions indels:
    1. the microhomology sequence must exist both rigth before the
    indel and downstream (when looking at the reference sequence).
    2. the distance between the two microhomilogy sequences is equal to
    the deletion's length
    optional:
    3. the sequence around the 2nd microhomology sequence in the reference
    wont exist in the individual's sequence
    In case of insertions, the pathway must go trhough the deletion step
    but with an alternative 2nd condition:
    2alt. The function will look for the following pattern -> repeat
    - microhomology - repeat.
    3. The length of the repeat must be smaller then the indel_len.
    """

    def microhomology_detection(self):
        self.seq_indel_spliter() # splitter

        # initializing general attributes in order to create a complete
        # attribute set for all indels, TS means tamplate_switch
        self.microhomology = False
        self.snap_out_dict = {}
        self.loop_out_dict = {}
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

        if self.indel_type == 'DEL': 
            self.deletion_MMEJ()


    def deletion_MMEJ(self):
        """
        Finds the longest subsequence (a motif) that exist both before and after
        the DSB, and ends at the last base of the deletion sequence
        """
        # set defult values for mmej_cand
        del_mmej = False
        del_mmej_cand, del_mmej_cand_len ,del_mmej_marked= np.nan, np.nan, np.nan
        del_last_dimer, mmej_cand = np.nan ,np.nan
        distance_to_indel_position, del_mmej_marked_on_ref = np.nan, np.nan
        # 1st condition:
        self.first_MH_finder() # finding all possible MH motifs
        self.first_MH_indicator()

        # 2nd condition
        # reverse interating over all motifs and look whether it exists
        # before the break or not
        for mh_motif in self.mmej_motifs[::-1]:
            if mh_motif == self.ref_genome_up[-len(mh_motif):]:
                del_mmej = True
                mmej_cand = del_mmej_cand = mh_motif
                mmej_cand_position = self.indel_length + 1
                self.distance_to_indel_position = self.indel_length - len(mmej_cand)
                break
            else: continue
        
        # get sequences with marked components of the mechanism pattern 
        if del_mmej:
            del_mmej_marked = self.del_accsession_mmej_marked(mmej_cand)
            del_mmej_marked_on_ref = self.del_mmej_marked_on_reference(mmej_cand,
                                       mmej_cand_position)
                
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

            
    def sd_snap_back_MMEJ(self):
        """
        sd_loop_out_MMEJ will detect SD-Snap back MEMJ by looking for the following pattern:
        [MH2->INS->P2] -> random seq (with length>0) -> invert_and_complement([MH2->INS->P2])
        This pattern is based on the one represented in fig.1B from:
        Khodaverdian, V. Y. et al. Secondary structure forming sequences drive SD-MMEJ 
        repair of DNA double-strand breaks. Nucleic Acids Res. 45, 12848–12861 (2017).

        Args: 
            self attributes of the MicroHomology object
        Returns:
            snap_mmej_marked (str): the whole pattern of SD-snap-back MMEJ 
                marked on the accession context
            snap_P1 (str): P1 as represented in Fig1B.
            snap_P2 (str): P2 as represented in Fig1B.
            snap_mh1 (str): MH1 as represented in Fig1B.
            snap_mh2 (str): MH2 as represented in Fig1B.
            snap_inv_comp_repeat_pat (str):The 1st repeat as represented in Fig1B
                (the one serves as a template for thr elongation of the loop).
            snap_repeat_pat (str): The 2nd repeat as represented in Fig1B
                (the one that get sythesized by elongation of the loop).
        """

        """
        1st step is to construct a small pattern as follow:
        1bp -> inv_comp(indel_sequence) -> 1bp
        we'll call this pattern: fondamental pattern (fond_pat) since this is
        the smallest pattern that can be defined as snap-back
        """

        mh2 = self.DSB_up[(-1):]
        mh1 = self.reverse_complement_converter(seq = mh2)
        P2 = self.DSB_down[(self.indel_length):(self.indel_length+1)] 
        P1 = self.reverse_complement_converter(seq = P2)
        
        fond_pat = f'{mh2}{self.indel_sequence}{P2}'
        inv_comp_fond_pat = self.reverse_complement_converter(seq=fond_pat)
        
        # making sure that inv_comp_fond_pat is in self.DSB_down
        # before start to find the complite pattern (MH1,MH2,P1 and P2)
        if inv_comp_fond_pat in self.DSB_down[len(inv_comp_fond_pat):]:
            SD_snap_back = True
        else: SD_snap_back = False
       
        # 2nd step is to find MH2 and MH1 (should also be invert-complementary 
        # to one another), conditioning the search range here to 20bp (longer
        # MH seq would be called SSA which is another mechanism)
        MH_max_range = 20
        if SD_snap_back :
            for i in range(1,MH_max_range):
                mh2 = self.DSB_up[(-i):]
                mh1 = self.reverse_complement_converter(seq = mh2)
                if f'{inv_comp_fond_pat[:-1]}{mh1}' in self.DSB_down[len(fond_pat): ]: continue
                # making sure that if the only pattern that exist is fond_pat,
                # we won't shorten mh1
                elif i==1: break
                else: 
                    # taking the last mh2 and mh1 that matched
                    mh2 = mh2[1:]
                    mh1 = mh1[:-1]
                    break  

        # 3rd step is to find P2 and P1 (should also be invert-complementary 
        # to one another) only if mh1 is found downstream
        # to the insertion end point
        
        if SD_snap_back :
            # start at the first nucleotide that comes right after the insertion,
            # calculate the invert-complementary version and see if it exist downstream
            # to the insertion end point, if true, elongate by 1bp

            for i in range(1,len(self.DSB_down)):
                P2 = self.DSB_down[
                    self.indel_length:(self.indel_length + i)]
                inv_comp_rep_pat = f'{mh2}{self.indel_sequence}{P2}'
                rep_pat = self.reverse_complement_converter(
                    seq = inv_comp_rep_pat)
                if rep_pat in self.DSB_down[
                    (self.indel_length + len(P2)): ]:
                    continue
                
                # making sure that if the only pattern that exist is fond_pat,
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
            rep_pat_pos = self.accession_sequence[(self.indel_position+len(P2)):].index(rep_pat)
            mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P2) + self.indel_length):(len(P2) + rep_pat_pos)]}'
            # set different variables as attributes
            snap_mmej_marked = self.snap_accession_pattern_generator(P1=P1, P2=P2, MH1=mh1, MH2=mh2,
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

        Args: 
            self attributes of the MicroHomology object

        returns:
            loop_mmej_marked (str): the whole pattern of SD-loop-out MMEJ 
                marked on the accession context
            loop_P2 (str): P2 as represented in Fig1A.
            loop_mh2 (str):mh2 (the actual syntetic MMEJ) as represented in Fig1A.
            loop_repeat_pat (str): The repeat that contains the insertion (1st rep
                in Fig1A).
        """  
        """
        1st step is to construct a small pattern as follow:
        1bp -> indel_sequence -> 1bp
        we'll call this pattern: fondamental pattern (fond_pat) since this is
        the smallest pattern that can be defined as snap-back
        """ 
        MH1 = self.DSB_up[(-1):]
        P2 = self.DSB_down[(self.indel_length):(self.indel_length+1)]
        fond_pat = MH1 + self.indel_sequence + P2

        # making sure that fond_pat is in self.DSB_down
        # before start to find the complite pattern (MH1,MH2,P1 and P2)
        if fond_pat in self.DSB_down[(len(self.indel_sequence) + len(P2) + 1):]:
            SD_loop_out = True
        else: SD_loop_out = False
        
        if SD_loop_out:
            """
            SD Loop Out demends the formation of a stem-loop structure.
            stem-loop structure cant form if the two MH->ins->P patterns
            are adjacent to one another.
            We must ensure that there is at least 1bp between the two.
            """
            # position of the MH->ins->P 2nd pattern
            fond_pat_pos = self.DSB_down[
                (len(self.indel_sequence) + len(P2)):].index(fond_pat)
            # extention_space is the longest sequence that mh2 can be found in 
            extention_space = self.DSB_down[
                (len(self.indel_sequence) + len(P2)):
                (len(self.indel_sequence) + len(P2) +fond_pat_pos)]
            
            # in cases were the extention_space is longer then 1,
            # enlarge fond_pat to longest possible one
            if len(extention_space) > 1:
                # elongation of mh2, confined in extention_space
                for i in range(1,len(extention_space)):
                    mh2 = self.DSB_up[(-i):]
                    fond_pat = mh2 + self.indel_sequence + P2

                    if fond_pat in self.DSB_down[
                        (len(self.indel_sequence) + len(P2)):]:
                        tmp_pos = self.DSB_down[
                            (len(self.indel_sequence) + len(P2)):].index(fond_pat)
                    else:
                        mh2 = self.DSB_up[(-i+1):]
                        fond_pat = mh2 + self.indel_sequence + P2
                        tmp_pos = self.DSB_down[
                            (len(self.indel_sequence) + len(P2)):].index(fond_pat)
                        break
                    if tmp_pos < (len(self.indel_sequence) + len(P2)): break

                # elongation of P2, confined in extention_space[:*starting position of mh2*]
                for i in range(1,len(extention_space[:(tmp_pos+1)])):
                    P2 = self.DSB_down[(self.indel_length):(self.indel_length+i)]
                    fond_pat = mh2 + self.indel_sequence + P2

                    if fond_pat in self.DSB_down[(len(self.indel_sequence) + len(P2)):]:
                        continue
                    else:
                        P2 = P2[:-1]
                        fond_pat = mh2 + self.indel_sequence + P2
                        break
            else: mh2 = MH1
            # Crearing mmej_marked
            rep_pat = fond_pat
            rep_pat_pos = self.accession_sequence[
                (self.indel_position+len(P2)):].index(rep_pat)
            mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P2) + self.indel_length):(len(P2) + rep_pat_pos)]}'
            
            # set output variables as attributes
            loop_mmej_marked = self.loop_accession_pattern_generator(P2=P2, MH2=mh2,
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
        Encapsulate all the relevant values into a dataframe object with standart
        format for all MMEJ pathways.
        Args:
            self MicroHomology inputs.
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
            'loop_last_dimer', 'loop_dist_between_reps']
        
        if self.include_context: _colnames = _colnames + ['ref_genome_context', 'accession_sequence']

        ex_df = pd.DataFrame(columns=_colnames)
        ex_df.loc[0,:] = np.nan
        # merging dicts from all paths
        d = {**self.del_out_dict, 
            **self.snap_out_dict, **self.loop_out_dict}
        if self.include_context: d = {**d, **self.context}

        # setting the right values to thiere coresponding columns
        ex_df.loc[0,[key for key in d.keys()]] = [val for val in d.values()]
        return ex_df


    # The function below are functions the generate a representation
    # of the different patterns per repair mechanism

    # Deletion MMEJ
    def del_accsession_mmej_marked(self, mmej_cand: str):
        """
        Create the pattern of a MMEJ that result in a deletion
        upstream seq -> MMEJ -> DSB -> DBS to MMEJ -> downstream seq
        Args:
            mmej_cand (str): thr MMEJ pattern
        Returns:
            The pattern itself.
        """
        if self.flip:
            return f'{self.DSB_down[0:10][::-1]}*[{mmej_cand[::-1]}]*{self.DSB_up[(len(self.DSB_up)-len(mmej_cand)-10):(len(self.DSB_up)-len(mmej_cand))][::-1]}'
        return f'{self.DSB_up[(len(self.DSB_up)-len(mmej_cand)-10):(len(self.DSB_up)-len(mmej_cand))]}*[{mmej_cand}]*{self.DSB_down[0:10]}'

    def del_mmej_marked_on_reference(self, mmej_cand: str,
                                       mmej_cand_position: int ):
        """
        Create the pattern of a MMEJ that result in a deletion
        on the reference genome:
        upstream seq -> MMEJ -> DSB -> DBS to MMEJ -> MMEJ -> downstream seq 
        Args:
            mmej_cand (str): thr MMEJ pattern
            mmej_cand_position (int): the position of the MMEJ seq 
                relative to the DSB
        Returns:
            The pattern itself.
        """
        before_DBS = f'{self.ref_genome_up[(len(self.ref_genome_up)-len(mmej_cand) -10):(len(self.ref_genome_up)-len(mmej_cand))]}*[{mmej_cand}]*|'
        DSB_to_MMEJ_seq = f'{self.ref_genome_down[0:(mmej_cand_position-len(mmej_cand)-1)]}'
        # print(mmej_cand_position, mmej_cand) 
        MMEJ_till_end = self.ref_genome_down[(mmej_cand_position-1):(mmej_cand_position+10)]
        if self.flip:
            before_DBS = f'*[{mmej_cand[::-1]}]*{self.ref_genome_up[(len(self.ref_genome_up)-len(mmej_cand) -10):(len(self.ref_genome_up)-len(mmej_cand))][::-1]}'
            return f'{MMEJ_till_end[::-1]}|*[{mmej_cand[::-1]}]*{DSB_to_MMEJ_seq[::-1]}{before_DBS}'
        return f'{before_DBS}{DSB_to_MMEJ_seq}*[{mmej_cand}]*{MMEJ_till_end}'


    # Insertion (SD-Snap back) MMEJ
    def snap_accession_pattern_generator(self, P1:str, P2: str, MH1:str, MH2: str,
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
            The pattern itself.
        """
        
        inv_ins = self.reverse_complement_converter(seq = self.indel_sequence)
        rep_pat_pos = self.accession_sequence[(self.indel_position+len(P2)):].index(rep_pat)
        mmej_marked_rep_up = f'*MH2[{MH2}]|*INS[{self.indel_sequence}]P2[{P2}]'
        mmej_marked_up = f'{self.accession_sequence[(self.indel_position-len(MH2)-10):(self.indel_position-len(MH2))]}{mmej_marked_rep_up}'
        mmej_marked_rep_down = f'*P1[{P1}][{inv_ins}]MH1[{MH1}]{self.DSB_down[(rep_pat_pos+len(rep_pat)+len(P2)):(rep_pat_pos+len(rep_pat)+len(P2)+10)]}'
        mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P2) + self.indel_length):(len(P2) + rep_pat_pos)]}'

        return f'{mmej_marked_up}{mmej_marked_inter_reps_seq}{mmej_marked_rep_down}'

    # Insertion (SD-Loop out) MMEJ
    def loop_accession_pattern_generator(self, P2: str, MH2: str,
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
        
        rep_pat_pos = self.accession_sequence[(self.indel_position+len(P2)):].index(rep_pat)
        mmej_marked_rep_up = f'*MH2[{MH2}]|*INS[{self.indel_sequence}]P2[{P2}]'
        mmej_marked_up = f'{self.accession_sequence[(self.indel_position-len(MH2)-10):(self.indel_position-len(MH2))]}{mmej_marked_rep_up}'
        mmej_marked_rep_down = f'*MH1[{MH2}][{self.indel_sequence}]P1[{P2}]{self.DSB_down[(rep_pat_pos+len(rep_pat)+len(P2)):(rep_pat_pos+len(rep_pat)+len(P2)+10)]}'
        mmej_marked_inter_reps_seq = f'{self.DSB_down[(len(P2) + self.indel_length):(len(P2) + rep_pat_pos)]}'

        return f'{mmej_marked_up}{mmej_marked_inter_reps_seq}{mmej_marked_rep_down}'

