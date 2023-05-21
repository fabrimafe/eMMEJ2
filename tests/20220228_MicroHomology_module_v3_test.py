# 28.02.2022
# importing the relevant libraries/modules:
import datetime
import timeit
import pandas as pd
import numpy as np
import sys
print('Libreries imported')
sys.path.append('/home/labs/alevy/guyta/guy_master_project/src')
from MicroHomology_module_v3 import MicroHomology
from MMEJ_2nd_order_MM_v2 import motif_probabily_calc
from data_exploration_util import *
pd.options.display.max_colwidth = 2400
print('_________________________________________________________________________________________________________')
print('Modules imported')
mem_reporter()

"""
This script is a test for the MicroHomology_module_v3 module
to make sure it runs properly and return the same expected 
results as I develop it. 
"""

def accession_context_generator(reference_contex : str,
                                alt : str, ref : str, window_size: int):
    """
    Generating the accession context by removing the ref sequence
    from the reference genome context and adding the alt sequence
    Args:
        reference_contex (str): the reference genome sequence
        ref (str): the anccestral allele
        alt (str): the derived allele
        window_size (int): the window size that was used in the
            vcftools getfasta process (1 based)
    Returns:
        accession_context (str): the mutated accession context
    """
    # matching window_size to be 0 indexed
    window_size = window_size-1
    # # Chacking that tha reference context match the REF
    reference_contex_REF = reference_contex[
        window_size:(len(ref)+window_size)]

    assert (reference_contex_REF == ref), f'REF dose not match the reference context, expected: {ref}, got: {reference_contex_REF}'
    # Generating accession context   
    pre_indel_seq = reference_contex[:(window_size)]
    # alt = alt[len(ref):]
    if len(alt)>len(ref):
        alt = alt[len(ref):]
        post_indel_seq = reference_contex[(window_size+1):] # +len(alt)+len(ref)
    else:
        post_indel_seq = reference_contex[(window_size+len(ref)):]
    # Chacking that accsession context was is correct
    accession_context = f'{pre_indel_seq}{alt}{post_indel_seq}'
    accession_context_check = accession_context[
        window_size:(len(alt)+window_size)]
    assert (accession_context_check == alt), f'accession_context_check dose not match ALT, expected: {alt}, got: {accession_context_check}'
    return accession_context

# Creating a fake dataset with fake examples
fake_df = pd.DataFrame(columns = ['ALT','REF', 'indel_type', 'indel_length', 'indel_pos', 
                                  'Read', 'reference_context_seq_2K', 'reference_context_seq_120'])

# Example 1: Deletion with high probability to be MMEJ
derived = 'C'
window_size = 120
indel_position = window_size -1 # converting window size to 0 base

fake_ref_120 = 'C'*115 + 'AGGT' + 'CGGCCAGGT' + 'C'*111
ancestrel = 'CGGCCAGGT'
indel_length = len(ancestrel) - len(derived)
fake_read = accession_context_generator(
    reference_contex = fake_ref_120, ref=ancestrel,alt=derived, window_size=120)

fake_ref_2K = 'C'*996 + 'AGGT' + 'CGGCC' + 'AGGT' + 'C'*991
fake_df.loc[0, ['ALT','REF', 'indel_type', 'indel_length','indel_pos',
               'Read', 'reference_context_seq_2K', 'reference_context_seq_120']] =  derived, ancestrel, 'DEL',indel_length, window_size, fake_read, fake_ref_2K, fake_ref_120


# Example 2: Deletion with low probability to be MMEJ
# fake_read = 'C'*20+'AGGT'*20+'C'*15 + 'AGGT' + 'C'*121

fake_ref_2K = 'C'*881 + 'C'*20+'AGGT'*20+'C'*15 + 'AGGT' + 'AGGT'*1 + 'C'*20 + 'AGGT'*10 + 'C'*936
fake_ref_120 = 'C'*20+'AGGT'*20+'C'*15 + 'AGGT' + 'CAGGT'*1 + 'C'*20 + 'AGGT'*10 + 'C'*57

ancestrel = 'CAGGT'
derived = 'C'
indel_length = len(ancestrel) - len(derived)
window_size = 120
indel_position = window_size -1 # converting window size to 0 base
fake_read = accession_context_generator(
    reference_contex = fake_ref_120, ref=ancestrel,alt=derived, window_size=120)
fake_df.loc[1, ['ALT','REF', 'indel_type', 'indel_length','indel_pos',
               'Read', 'reference_context_seq_2K', 'reference_context_seq_120']] =  derived, ancestrel, 'DEL',indel_length, window_size, fake_read, fake_ref_2K, fake_ref_120


# Example 3: Insertion (Trans) with High probability to be MMEJ
# working version ###################
pre_indel = 'C'*115
MH1 = 'AGGT'
rep = 'CCAG'
MH2 = 'ACG'
post_pattern = 'C'*110
fake_read = f'{pre_indel}{MH1}{rep}{MH2}{rep}{post_pattern}'
fake_ref_2K = 'C'*996 + f'{MH1}{MH1}{rep}{MH2}' + 'C'*989
fake_ref_120 = f'{pre_indel}{MH1}{MH1}{rep}{MH2}{post_pattern}'
fake_df.loc[2, ['ALT','REF', 'indel_type', 'indel_length','indel_pos',
               'Read', 'reference_context_seq_2K', 'reference_context_seq_120']] =  'CCCCCAGGT','', 'INS',len('CCCCCAGGT'), 119, fake_read, fake_ref_2K, fake_ref_120
# working version ###################

pre_indel = 'C'*115
MH1 = 'AGGT'
rep = 'CCAG'
MH2 = 'ACG'
post_pattern = 'C'*110
fake_read = f'{pre_indel}{MH1}{rep}{MH2}{rep}{post_pattern}'
fake_ref_2K = 'C'*996 + f'{MH1}{MH1}{rep}{MH2}' + 'C'*989
fake_ref_120 = f'{pre_indel}{MH1}{rep}{MH2}{post_pattern}'
# print(f'fake_read: {fake_read}')
# print(f'fake_ref_120: {fake_ref_120}')

derived = 'CCAGACG'
ancestrel = 'C'
indel_length = len(derived) - len(ancestrel)
window_size = 120
indel_position = window_size -1 # converting window size to 0 base
fake_read = accession_context_generator(
    reference_contex = fake_ref_120, ref=ancestrel,alt=derived, window_size=window_size)
# print(f'fake_read: {fake_read}')
fake_df.loc[2, ['ALT','REF', 'indel_type', 'indel_length','indel_pos',
               'Read', 'reference_context_seq_2K', 'reference_context_seq_120']] =  derived, ancestrel, 'INS',indel_length, 119, fake_read, fake_ref_2K, fake_ref_120

# Example 4: Insertion (Trans) with low probability to be MMEJ
pre_indel = 'C'*115
MH1 = 'AGGT'
rep = 'CCAG'
MH2 = 'ACG'
post_pattern = 'C'*110
fake_read = f'{pre_indel}{MH1}{rep}{MH2}{rep}{post_pattern}'
fake_ref_2K = ('CCAGACGCCAG' * 180) + 'C'*20
fake_ref_120 = f'{pre_indel}{MH1}{MH1}{rep}{MH2}{post_pattern}'
fake_df.loc[3, ['ALT','REF', 'indel_type', 'indel_length','indel_pos',
               'Read', 'reference_context_seq_2K', 'reference_context_seq_120']] =  'CCCCCAGGT','', 'INS',len('CCCCCAGGT'), 119, fake_read, fake_ref_2K, fake_ref_120


# Example 5: Insertion (Snap back) with High probability
# fake_read = 'AAAAAAAA'*10+'C'*35 + 'AGGT' + 'TGGA' +'CCCCCCCCCC'+'AAAAAAAA'+ 'GGGGGGGGGG'+ 'TCCA'+'ACCT' + 'AAAAAAAA'*9 + 'C'*9
# fake_ref_2K = ('A'*75)+('AAAAAAAA'*10+'C'*35)*8+'AGGT'+('CCCCCCCCCC'+'AAAAAAAA'+ 'GGGGGGGGGG'+ 'TCCA'+'ACCT' + 'AAAAAAAA'*9 + 'C'*9)*8+('C'*65)
# fake_ref_120 = 'AAAAAAAA'*10+'C'*35 + 'AGGT' +'CCCCCCCCCC'+'AAAAAAAA'+ 'GGGGGGGGGG'+ 'TCCA'+'ACCT' + 'AAAAAAAA'*9 + 'C'*9
# # fake_read_example1 = 'AAAAAAAA'*10+'C'*35 + ' MH2[AGGT]' + ' INS[TGGA]' +' P2[CCCCCCCCCC]'+'AAAAAAAA'+ ' P1[GGGGGGGGGG]'+ '[TCCA]'+' MH1[ACCT]' + 'AAAAAAAA'*9 + 'C'*9
# # print('fake_read_example1:', fake_read_example1,'\n')
# fake_df.loc[4, ['ALT','REF', 'indel_type', 'indel_length','indel_pos',
#                'Read', 'reference_context_seq_2K', 'reference_context_seq_120']] =  'TGGA','', 'INS',len('TGGA'), 119, fake_read, fake_ref_2K, fake_ref_120


fake_ref_2K = ('A'*75)+('AAAAAAAA'*10+'C'*35)*8+'AGGT'+('CCCCCCCCCC'+'AAAAAAAA'+ 'GGGGGGGGGG'+ 'TCCA'+'ACCT' + 'AAAAAAAA'*9 + 'C'*9)*8+('C'*65)
fake_ref_120 = 'AAAAAAAA'*10+'C'*35 + 'AGGTT' +'CCCCCCCCCC'+'AAAAAAAA'+ 'GGGGGGGGGG'+ 'TCCA'+'ACCT' + 'AAAAAAAA'*9 + 'C'*9
# print(f'fake_r120: {fake_ref_120}')
derived = 'TTGGA'
ancestrel = 'T'
indel_length = len(derived) - len(ancestrel)
window_size = 120
indel_position = window_size -1 # converting window size to 0 base
fake_read = accession_context_generator(
    reference_contex = fake_ref_120, ref=ancestrel,alt=derived, window_size=window_size)
fake_df.loc[4, ['ALT','REF', 'indel_type', 'indel_length','indel_pos',
               'Read', 'reference_context_seq_2K', 'reference_context_seq_120']] =  derived, ancestrel, 'INS',indel_length, 119, fake_read, fake_ref_2K, fake_ref_120


# Example 6: Insertion (Loop-out) with High probability
fake_read = 'AAAAAAAA'*10+'C'*35 + 'AGGT' + 'TGGA' +'CCCCCCCCCC'+'AAAAAAAA'+ 'AGGT' +'TGGA'+ 'CCCCCCCCCC'+ 'TAAAAAAA'*9 + 'C'*9
fake_ref_2K = ('A'*75)+('AAAAAAAA'*10+'C'*35)*8+'AGGT' + ('CCCCCCCCCC'+'AAAAAAAA'+ 'AGGT'+'TGGA' + 'CCCCCCCCCC'+ 'AAAAAAAA'*9 + 'C'*9)*8+('C'*65)
fake_ref_120 = 'AAAAAAAA'*10+'C'*35 + 'AGGT' +'CCCCCCCCCC'+'AAAAAAAA'+'AGGT'+'TGGA' + 'CCCCCCCCCC'+ 'TAAAAAAA'*9 + 'C'*9
# fake_read_example2 = 'AAAAAAAA'*10+'C'*35 + ' MH2[AGGT]' + ' INS[TGGA]' +' P2[CCCCCCCCCC]'+'AAAAAAAA'+ ' MH1[AGGT]' +'TGGA'+ ' P1[CCCCCCCCCC]'+ 'TAAAAAAA'*9 + 'C'*9
# print('fake_read_example2:', fake_read_example2,'\n')
fake_df.loc[5, ['ALT','REF', 'indel_type', 'indel_length','indel_pos',
               'Read', 'reference_context_seq_2K', 'reference_context_seq_120']] =  'TGGA','', 'INS',len('TGGA'), 119, fake_read, fake_ref_2K, fake_ref_120

df = fake_df

assert (df.shape == (6, 8)), f'Wrong df.shape, Expected : (6, 8), Got: {df.shape}'

# print(df.shape)

"""
Detecting microhomology:
The data must have the following columns:
# accession_context -> A sequence with the indel and no gaps (type: str)
# derived_allele -> The indel itself (i.e. the sequence that got inserted/deleted)
# ancesstral_allele -> The reference sequence at the indel position 
    (in a case of insertion the value should be: "", in case of deletion 
    the value should be: the deleted sequence)
# reference_context_seq_120 -> the reference genome at a window of +-120bp from the indel position
    (i.e. -120bp ->INDEL->+ 120bp)
# indel_pos -> the indel position in the read
# indel_length -> the length of the indel
# indel_type -> the type of the indel (INS\DEL)
"""
microhomology_object = df.apply( lambda x: MicroHomology(x['Read'],
                                                            x['ALT'],
                                        x['REF'],
                                        x['reference_context_seq_120'],
                                        x['indel_pos'],
                                        # x['indel_length'],
                                        x['indel_type']), axis = 1)

print('microhomology_object assinged')

df2 = None
for i in microhomology_object.index:
    if (i == 0):
        df2 = microhomology_object[i].ex_data
    else:
        df2 = pd.concat([df2, microhomology_object[i].ex_data])

print('microhomology_object.ex_data extracted')

df2.reset_index(inplace = True)
df = df.join(df2, how = 'right')
df.drop(columns = ['index'], inplace = True)
print('Data and microhomology_objects joined')

# creating empty columns for the probabilities
df['del_mmej_prob'] = np.nan
df['trans_mmej_prob'] = np.nan
df['SD_snap_back_prob'] = np.nan
df['SD_loop_out_prob'] = np.nan
# df['SD_snap_back_prob'] = np.nan
# reassining 'None' values to be False when needed for ferther filtrarion
df['del_mmej'].replace(to_replace = [None], value = False, inplace = True)
df['SD_trans'].replace(to_replace = [None], value = False, inplace = True)
df['SD_snap_back'].replace(to_replace = [None], value = False, inplace = True)
df['SD_loop_out'].replace(to_replace = [None], value = False, inplace = True)
"""
The probability to find mmej will be calculate seperatly since each scenario has a different
motif. Other then the motif itself, the rest is similar.
"""
print('Markov model start')
mem_reporter()
# creating a df with only deletions
df_with2K_del_mmej = df[((df['del_mmej'] == True) & 
           (df['indel_type'] == 'DEL'))]

print('DEL mmej probability (Markov) starts')
df.loc[
    ((df['del_mmej'] == True) & (df['indel_type'] == 'DEL')),
    'del_mmej_prob'] = df_with2K_del_mmej.apply(lambda x: motif_probabily_calc(
                    _sequance = x['reference_context_seq_2K'] , 
                                 _motif = x['del_mmej_cand'],
                                _indel_len = int(x['indel_length']),
                    _memory_dimer = x['del_last_dimer']) , axis = 1 )
print('DEL mmej probability (Markov) ends')

# creating a df with only insertions
df_with2K_ins_mmej = df[(df['SD_trans'] == True)].copy()
df_with2K_ins_mmej.loc[:,'trans_mmej_cand_len'] = df_with2K_ins_mmej['trans_mmej_cand_len'].astype('int')
print('INS mmej probability (Markov) starts')
df.loc[
    (df['SD_trans'] == True),
    'trans_mmej_prob'] = df_with2K_ins_mmej.apply(lambda x: motif_probabily_calc(
                    _sequance = x['reference_context_seq_2K'] , 
                                 _motif = x['trans_TS_repeat'],
                                _indel_len = x['trans_mmej_cand_len'],
                    _memory_dimer = x['trans_last_dimer']), axis = 1 )
print('INS mmej probability (Markov) ends')

# creating a df with only mmej snap back
df_with2K_mmej_snap_back = df[df['SD_snap_back'] == True].copy()

print('snap back mmej probability (Markov) starts')
df.loc[
    (df['SD_snap_back'] == True),
    'SD_snap_back_prob'] = df_with2K_mmej_snap_back.apply(
    lambda x: motif_probabily_calc(_sequance = x['reference_context_seq_2K'] , 
                                 _motif = x['snap_repeat_pat'],
                                _indel_len = x['snap_dist_between_reps'],
                    _memory_dimer = x['snap_last_dimer']), axis = 1 )

print('snap back mmej probability (Markov) ends')

# creating a df with only mmej Loop out
df_with2K_mmej_loop_out = df[df['SD_loop_out'] == True].copy()

print('Loop out mmej probability (Markov) starts')
df.loc[
    (df['SD_loop_out'] == True),
    'SD_loop_out_prob'] = df_with2K_mmej_loop_out.apply(
    lambda x: motif_probabily_calc(_sequance = x['reference_context_seq_2K'] , 
                                 _motif = x['loop_repeat_pat'],
                                _indel_len = x['loop_dist_between_reps'],
                    _memory_dimer = x['loop_last_dimer']), axis = 1 )

print('Loop out mmej probability (Markov) ends')
print('Markov model ends')
mem_reporter()

"""
Sice the Markov model calculates the probability to find 2 motifs in a
defined range and a given context, this probability is the probability
that something isn't MMEJ, therefore if one wants to have the probability
of something to be a MMEJ, he should do 1-probability
"""
df['del_mmej_prob'] = 1 - df['del_mmej_prob']
df['trans_mmej_prob'] = 1 - df['trans_mmej_prob']
df['SD_snap_back_prob'] = 1 - df['SD_snap_back_prob']
df['SD_loop_out_prob'] = 1 - df['SD_loop_out_prob']

# pd.set_option("display.max_columns")
print(df)
"""
Testing results:
"""
# setting expected mmej_cand results:
mem_reporter()
print('\n##### TEST RESULTS #####')
# mmej
# deletions
expected_del_mmej = [True, True, False, False, False, False]
 
assert (list(df['del_mmej']) == expected_del_mmej), f"""Wrong del_mmej results.
 Expected {expected_del_mmej}, Got: {list(df['del_mmej'])}"""
print('#del_mmej Passed')

# trans
expected_SD_trans = [False, False, True, True, False, False]
 
# assert (list(df['SD_trans']) == expected_SD_trans), f"""Wrong SD_trans results.
#  Expected {expected_SD_trans}, Got: {list(df['SD_trans'])}"""
# print('#SD_trans Passed')

# snap-back
expected_snap = [False, False, False, False, True, False]
 
assert (list(df['SD_snap_back']) == expected_snap), f"""Wrong SD_snap_back results.
 Expected {expected_snap}, Got: {list(df['SD_snap_back'])}"""
print('#SD_snap_back Passed')

# loop- out
expected_loop_mmej = [False, False, False, False, False, True]
 
assert (list(df['SD_loop_out']) == expected_loop_mmej), f"""Wrong SD_loop_out results.
 Expected {expected_loop_mmej}, Got: {list(df['SD_loop_out'])}"""
print('#SD_loop_out Passed')

# mmej_cand_del
expected_del_mmej_res = ['CCAGGT', 'AGGT', np.nan, np.nan, np.nan, np.nan]
assert (list(df['del_mmej_cand']) == expected_del_mmej_res), f"""Wrong del_mmej_cand results. Expected {expected_del_mmej_res},
                                                        Got: {list(df['del_mmej_cand'])}"""

# # trans_reps_pat
# expected_trans_reps_pat = [np.nan, np.nan, 'CCAGACGCCAG', 'CCAGACGCCAG', np.nan, np.nan]
# assert (list(df['trans_reps_pat']) == expected_trans_reps_pat), f"""Wrong trans_reps_pat results. Expected {expected_trans_reps_pat},
#                                                         Got: {list(df['trans_reps_pat'])}"""

# mmej_cand_snap
expected_snap_repeat_pat = [np.nan, np.nan, np.nan, np.nan, 'GGGGGGGGGGTCCAACCT', np.nan]
assert (list(df['snap_repeat_pat']) == expected_snap_repeat_pat), f"""Wrong snap_repeat_pat results. Expected {expected_snap_repeat_pat},
                                                        Got: {list(df['snap_repeat_pat'])}"""

# mmej_cand_loop
expected_loop_repeat_pat = [np.nan, np.nan, np.nan, np.nan, np.nan, 'AGGTTGGACCCCCCCCCC']
assert (list(df['loop_repeat_pat']) == expected_loop_repeat_pat), f"""Wrong loop_repeat_pat results. Expected {expected_loop_repeat_pat},
                                                        Got: {list(df['loop_repeat_pat'])}"""

print('#del_mmej_cand, trans_reps_pat, snap_repeat_pat, loop_repeat_pat Passed')


# del_last_dimer
expected_del_last_dimer = ['GT', 'GT', np.nan, np.nan, np.nan, np.nan]
assert (list(df['del_last_dimer']) == expected_del_last_dimer), f"""Wrong del_last_dimer results.
 Expected {expected_del_last_dimer}, Got: {list(df['del_last_dimer'])}"""

# # trans_last_dimer
# expected_trans_last_dimer = [np.nan, np.nan, 'AG', 'AG', np.nan, np.nan]
# assert (list(df['trans_last_dimer']) == expected_trans_last_dimer), f"""Wrong trans_last_dimer results.
#  Expected {expected_trans_last_dimer}, Got: {list(df['trans_last_dimer'])}"""

# snap_last_dimer
expected_snap_last_dimer = [np.nan, np.nan, np.nan, np.nan, 'CC', np.nan]
assert (list(df['snap_last_dimer']) == expected_snap_last_dimer), f"""Wrong snap_last_dimer results.
 Expected {expected_snap_last_dimer}, Got: {list(df['snap_last_dimer'])}"""

# loop_last_dimer
expected_loop_last_dimer = [np.nan, np.nan, np.nan, np.nan, np.nan, 'CC']
assert (list(df['loop_last_dimer']) == expected_loop_last_dimer), f"""Wrong loop_last_dimer results.
 Expected {expected_loop_last_dimer}, Got: {list(df['loop_last_dimer'])}"""


print('last_dimer Passed')

# del_mmej_marked_on_ref
expected_mmej_marked_on_ref = ['CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*[CCAGGT]*|CGG*[CCAGGT]*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC', 
    'GGTAGGTAGGTAGGTCCCCCCCCCCCCCCC*[AGGT]*|C*[AGGT]*CCCCCCCCCCCCCCCCCCCAGGTAGGTAGG', np.nan, np.nan, np.nan, np.nan]
assert (list(df['del_mmej_marked_on_ref']) == expected_mmej_marked_on_ref), f"""Wrong del_mmej_marked_on_ref results.
 Expected {expected_mmej_marked_on_ref}, Got: {list(df['del_mmej_marked_on_ref'])}"""
print('#del_mmej_marked_on_ref Passed')

# del_mmej_marked
expected_del_mmej_marked = ['CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*[CCAGGT]*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC', 
    'GGTAGGTAGGTAGGTCCCCCCCCCCCCCCC*[AGGT]*CCCCCCCCCCCCCCCCCCCCCAGGTAGGTA', np.nan, np.nan, np.nan, np.nan]
assert (list(df['del_mmej_marked']) == expected_del_mmej_marked), f"""Wrong del_mmej_marked results.
 Expected {expected_del_mmej_marked}, Got: {list(df['del_mmej_marked'])}"""
print('#del_mmej_marked Passed')

# # trans_mmej_marked
# expected_trans_mmej_marked = [np.nan, np.nan, 'CCCCCCCCCCCCCCCCCCCCCCCCCCAGGT|*rep1[CCAG]*ACG*rep2[CCAG]*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC', 
#             'CCCCCCCCCCCCCCCCCCCCCCCCCCAGGT|*rep1[CCAG]*ACG*rep2[CCAG]*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC', np.nan, np.nan]
# assert (list(df['trans_mmej_marked']) == expected_trans_mmej_marked), f"""Wrong trans_mmej_marked results.
#  Expected {expected_trans_mmej_marked}, Got: {list(df['trans_mmej_marked'])}"""
# print('#trans_mmej_marked Passed')

# snap_mmej_marked
expected_snap_mmej_marked = [np.nan, np.nan, np.nan, np.nan, 
    'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*MH2[AGGT]|*INS[TGGA]P2[CCCCCCCCCC]AAAAAAAA*P1[GGGGGGGGGG][TCCA]MH1[ACCT]AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
     np.nan]
assert (list(df['snap_mmej_marked']) == expected_snap_mmej_marked), f"""Wrong snap_mmej_marked results.
 Expected {expected_snap_mmej_marked}, Got: {list(df['snap_mmej_marked'])}"""
print('#snap_mmej_marked Passed')

# loop_mmej_marked
expected_loop_mmej_marked = [np.nan, np.nan, np.nan, np.nan, np.nan,
 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*MH2[AGGT]|*INS[TGGA]P2[CCCCCCCCCC]AAAAAAAA*MH1[AGGT][TGGA]P1[CCCCCCCCCC]TAAAAAAATAAAAAAATAAAAAAATAAAAA']
assert (list(df['loop_mmej_marked']) == expected_loop_mmej_marked), f"""Wrong loop_mmej_marked results.
 Expected {expected_loop_mmej_marked}, Got: {list(df['loop_mmej_marked'])}"""
print('#loop_mmej_marked Passed')

print('######### ALL TESTS PASSED ########')