# 29.06.2022

import pandas as pd
from pysam import FastaFile

def get_ref_context(refFA: FastaFile,  chr: str, 
                     indel_seq: str, indel_pos: int,
                     context_window_size: int) -> str:
    """
    get fasta context from a reference genome fasta file that was previosly 
    indexed using tabix
    Args:
        refFA (FastaFile): a FastaFile object (from pysam)
        chr (str): chromosome name as in the reference file 
        indel_seq (str): indel sequence
        indel_pos (int): indel position
        context_window_size (int): number of bases +/- the indel
    Returns:
        (srt): the reference sequence +/- window size 
    """
    return refFA.fetch(chr, (indel_pos - context_window_size),
                             (indel_pos + len(indel_seq) + context_window_size+1))


def accession_context_generator(reference_contex : str,
                                alt : str, ref : str, window_size: int):
    """
    Reconstruct the accession sequence 
    Args:
        reference_contex (str): the reference genome around the indel 
        alt (str): the indel ALT sequence
        ref (str): the indel REF sequence
        window_size (int): number of bases +/- the indel
    Returns:
        accession_context (str): the accession sequence 
    """
    # matching window_size to be 0 indexed
    window_size = window_size-1
    # # Chacking that tha reference context match the REF
    # reference_contex_REF = reference_contex[
    #     window_size:(len(ref)+window_size)]

    # assert (reference_contex_REF == ref), f'REF dose not match the reference context, expected: {ref}, got: {reference_contex_REF}, Indel position (0 based): {window_size},alt: {alt}, reference_contex: {reference_contex}'
    # Generating accession context
    per_indel_seq = reference_contex[:(window_size)]
    post_indel_seq = reference_contex[(window_size + len(ref)):]
    
    # Chacking that accsession context was is correct
    accession_context = f'{per_indel_seq}{alt}{post_indel_seq}'
    # accession_context_check = accession_context[
    #     window_size:(len(alt)+window_size)]
    # assert (accession_context_check == alt), f'accession_context_check dose not match ALT, expected: {alt}, got: {accession_context_check}'
    
    return accession_context
