# 06.02.2022
# importing libraries
import regex as re
import itertools

import pandas as pd
import numpy as np

# from VcfProcess_module import VcfProcess
# from MicroHomology_module_v3 import MicroHomology
from helper import type_checker
import typing

def get_motifs_pos(context_seq: str, motif: str) -> list:
    """
    Find all the positions of which the motif is accure (including
    overlaps), then calculate the distance berween each location
    to the next one.
    Args:
        context_seq (str): the context sequence of the indel
        motif (str): the motif sequence
    Returns:
        motifs_dist_ls (list): a list with the distances
    """
    motifs_inx = np.array([m.start() for m in 
            re.finditer(motif, context_seq, overlapped=True)])
    motifs_end_inx = motifs_inx + len(motif)
    motifs_dist_ls = []
    
    for i in range(len(motifs_inx)):        
        next_pos_ls = motifs_inx[motifs_inx>=motifs_end_inx[i]]
        if len(next_pos_ls)>0: 
            next_pos = next_pos_ls[0]
            motifs_dist_ls.append(next_pos-motifs_end_inx[i])
        else: break
    return np.array(motifs_dist_ls)

def get_q(context_seq: str, motif: str,
            window_size: int, obs_pos: int) -> float: 
    """
    Args:
        contex_seq (str): the context sequence
        motif (str): the motif sequence
        window_size (int): last position to look for the motif (we set 50)
        obs_pos (int): position of the observed motif
    Returns:
        q (float): 
    """
    context_seq = context_seq[round(len(context_seq)/2):round(len(context_seq)/2+window_size)]
    # Generate a fake vector that indicate the positions at which the motif was found
    # 1.
    mmej_pos = np.zeros(window_size)
    # 2.
    mmej_motif_pos = np.array([m.start() for m in 
            re.finditer(motif, context_seq, overlapped=True)])
    mmej_motif_pos = mmej_motif_pos[mmej_motif_pos<window_size]
    mmej_pos[mmej_motif_pos] = 1
    # 3. Normalize
    occurence_n = len(mmej_pos[mmej_pos == 1])
    mmej_pos = mmej_pos/occurence_n
    
    # 4. and 5. 
    NHEJ_dist = np.array([1/i for i in range(1,window_size+1)])/window_size
    
    q = mmej_pos[obs_pos]/ (mmej_pos[obs_pos] + NHEJ_dist[obs_pos])
    return q


def get_nmer_combinations(n_mer: int):
        """
        This function will take the length of the kmer as an input
        and will generate a list of all possible combinations in k
        length, composed from the four nucleotides.
        Args:
            n_mer (int): length of the kmer
        """
        nuc = ['A','T','C','G']
        # creating an empty list
        n_mer_list = []
        # creating all possible kmers
        for output in itertools.product(nuc, repeat= n_mer):
            n_mer_list.append(''.join(output))
        n_mer_list = sorted(n_mer_list, key=lambda x: x[-1])
        return n_mer_list

def get_kmers_counts(seq: str, k_size: int):
    """
    Takes a sequence and count the appearance of all
    possible kmers (of size=k_size+1) in this sequence
    Args:
        seq (str): a DNA sequence
        k_size (int): length of the kmer
    Returns:
        pd.Series with the kmers and thier counts
    """
    sequence = seq.upper()
    # k_size = int(k_size)
    # memmory allocating : vector in length of the sequence, for kmers
    kmers = np.zeros(len(sequence), dtype = 'object')
    # spliting the sequence into kmers
    for i in range(0, (len(sequence) - 1)):
        if i>(len(sequence)-(k_size + 1)): break
        # generating a kmer that is twice the length of the desiered kmer
        else: kmers[i] = sequence[i:(i+(k_size + 1))]
    # calculating the freq of each twice long kmer in the sequence
    return (pd.Series(kmers).value_counts().drop(labels=[0]))

def get_equilibrium_state(sequence: str, n_mer=2 ):
    """
    Calculates the dimer frequency in the context
    which is the equivalent to the equilibrium state
    """
    
    dimers = get_nmer_combinations(n_mer=n_mer)
    dimers_freq = pd.DataFrame(index=dimers)
    dimers_freq['counts'] = dimers_freq.index.map(lambda d: np.array([m.end() for m in 
                re.finditer(d, sequence, overlapped=True)]).shape[0]).T
    dimers_freq['dimers_freq'] = (dimers_freq['counts']/(len(sequence)-1))+10**(-6)
    return dimers_freq['dimers_freq']

def kmers_counts_vector2dataframe(kmers_counts_vec: np.array, 
                k_size: int ) -> pd.DataFrame:
    """
    Assigning kmers counts from a vector into dataframe
    Args:
        kmers_counts_vec (pd.Series): a vector with counts per kmer
            in the sequence
        k_size (int): kmer length
    Returns:
        A pd.DataFrame with the kmers counts
    """
    # creating an empty kmers_freq_table
    kmers_freq_table = pd.DataFrame(index= get_nmer_combinations(k_size),
                    columns= get_nmer_combinations(k_size-1))
    # filling the empty kmers_freq_table
    for i in range(0,len(kmers_counts_vec)):
        value = kmers_counts_vec[i]
        if 'N' in kmers_counts_vec.index[i]:
            pass
        else:
            row = kmers_counts_vec.index[i][:k_size]
            col = kmers_counts_vec.index[i][k_size:]
            kmers_freq_table.loc[row,col] = value
    return kmers_freq_table

def norm_TM(TM: pd.DataFrame) -> pd.DataFrame:
    """
    Normalizing the transition matrix.
    Each sum(row)==1
    Args:
        TM (pd.DataFrame): transition matrix
    """
    TM.fillna(0, inplace = True)
    # adding 10**(-6) to the whole matrix
    TM = TM + 10**(-6)
    # normalizing by column (each column sums to 1)
    return TM.apply(lambda x: x/x.sum(), axis=1)


def TM_chack(norm_TM: pd.DataFrame):
    """
    Basicaly the probability of moving from any kmer (row, say AG)
    to any nucleotide should sum up to 1 (that is, sum(AG->A, AG->T, AG->C, AG->G))
    Args:
        norm_TM (pd.DataFrame): a normalized transition matrix
    """
    norm_TM_row_sum = norm_TM.sum(axis=1).round(decimals=1)
    norm_TM_row_sum = (norm_TM_row_sum != 1.00)      
    assert any(norm_TM_row_sum) != 1, f"""Transition matrix normalization problem
    The following row does not sums to 1.0
    {norm_TM_row_sum}
    Transition Matrix:
    {norm_TM_row_sum}"""
    
def get_TM(seq: str, k_size: int) -> pd.DataFrame:
    """
    Generates a transition matrix based on a sequence
    Args:
        seq (str): the DNA sequence on whitch the transition
            matrix will be calculated.
    Returns:
        TM (pd.DataFrame): the transition matrix
    """
    kmers_counts_vec = get_kmers_counts(seq=seq, k_size=k_size)
    kmers_freq_table = kmers_counts_vector2dataframe(
                kmers_counts_vec=kmers_counts_vec, k_size=k_size)
    TM = norm_TM(TM=kmers_freq_table)
    TM_chack(TM)

    return TM


# def get_motif_dimers(motif: str) -> list: # , memory_dimer: str
#     """
#     split the motif's sequence into dimers with propogating step=1bp
#     Args:
#         motif (str): the motif's sequence
#     Returns:
#         motif_dimers (list): a list of dimers from the motif
#     """
#     motif_dimers = []
#     # adding dimers from the motif's sequence, progressing in steps
#     # of 1 base per iteration
#     [motif_dimers.append(motif[i : i+2]) for i in range(0,len(motif)-1)]
#     return motif_dimers


# def get_motif_prob(motif: str,
#                     motifs_d: int, 
#                     context_seq: str,
#                     TM: pd.DataFrame) -> float:
#     """
#     Calculating the probability to find a motif given a
#     transition matrix (that is based on the reference genome)
#     Args:
#         motif (str): the motif sequence
#         motifs_d (int): distance between two motifs
#         TM (pd.DataFrame): a transition matrix (dimer -> monomer)
#     Returns:
#         motif_prob (float): the probability to find the motif
#     """
#     motif_prob = 1
#     if (motifs_d==0): 
#         motif_prob = len(re.findall((motif*2), 
#                     context_seq, overlapped=True))/len(context_seq) 
#     elif len(motif) <= 2: pass
#     elif len(motif) > 2:
#         # iterating through both Bn (bases in the motif) and (Bn-1,Bn-2) 
#         # (from motif_dimers)
#         for b in range(len(motif)-2):
#             # calculating the product of P(Bn|Bn-1,Bn-2)
#             # print(f'P({motif[b+2:b+3]}|{motif[b:b+2]})')
#             # print(TM.loc[motif[b:b+2],motif[b+2:b+3]])
#             motif_prob = motif_prob * TM.loc[motif[b:b+2],motif[b+2:b+3]]
#     # print(f'motif prob P(M) = {motif_prob}')
#     return motif_prob

# def get_motif_first_dimer_state_prob(state_df: pd.DataFrame, 
#                             motif_first_dimer: str, motif: str) -> list:
#     """
#     Getting the probability to have the memory_dimer at each state
#     Args:
#         state_df (pd.DataFrame): a dataframe with all states prob
#             distribution 
#         motif (str): the motif sequence
#         motif_1st_dimer (str): the motif's 1st dimer
#     Returns:
#         a list of probability to have the memory_dimer at each state
#     """
#     # getting P(motif's 1st kmer at S(n)) from the states matrix
#     if len(motif) > 2: 
#         return [state_df.loc[motif_first_dimer,i] for i in state_df.columns]
#     elif len(motif) == 2: 
#         return [state_df.loc[motif,i] for i in state_df.columns]
#     elif len(motif) == 1: 
#         return [state_df.loc[state_df.index.str.endswith(motif),i].sum() for i in state_df.columns]


# def get_motif_prob_per_state(first_dimer_state_prob: list, 
#                         motif_prob: float) -> pd.Series:
#     """
#     Calculating: 
#     P(having the memory_dimer at each state) * P(motif seq)
#     Args:
#         first_dimer_state_prob (list):a list of probability to have 
#                         the motif's 1st dimer at each state
#         motif_prob (float): the probability to find the motif
#     Returns:
#         1-(motif prob per state) (np.array)
#     """      
#     return np.absolute(1 - (np.array(first_dimer_state_prob)* motif_prob))


# def get_motif_prob_in_range(motif_prob_per_state: list):
#     """
#     Calculating the probability to find the motif at the range(S(0),S(n))
#     Formula:
#     let Px(n) be the probability to find the motif in a specific state 
#     (location) such that:
#         Px(n) =  P(motif's last kmer at S(n))*P(motif_prob)
#     as : 1 - [ ∏ (1-Px(n))], (n=0 : n=last state number)
#     Args:
#         motif_prob_per_state (np.array): an array with prob of
#             not finding the motif at each state
#     Return:
#         The prob to find the motif at any distance from the 1st motif
#         starting at 0 up until the last state. (float)
#     """
#     return np.prod(motif_prob_per_state)

    
def get_markov_states(transition_matrix: pd.DataFrame, 
                            first_kmer: str,
                            number_of_states: int,
                            eq_state=None, eq=False) -> pd.DataFrame:
        """
        get_markov_states will take a fasta context, motif and a
        number of states. The function will construct a 2nd order markov chain
        and produce a dataframe with n states of the chain.
        """
        # seting the memory of the model
        memory = 2
        first_kmer = first_kmer.upper()
        # seting the first kmer to start with
        first_kmer = first_kmer[-memory:]
        M = transition_matrix
        if eq==True:
            S0=eq_state
        else:
            # Creating the first state vector (shape 16*1, all 0)
            S0 = pd.Series(np.zeros([M.shape[0]]),
                        index = M.index)
            # assigning the first state to 1
            S0[first_kmer] = 1
        # creating an empty dataframe for the states 
        states_df = pd.DataFrame(S0)
        # creating current state vector (shape 16*1)
        Sn = pd.Series(np.zeros([M.shape[0]]),
                      index = M.index)

        # This is the main loop of the function
        for i in range(1,number_of_states):
            """
            The idea here is that in each iteration, we look at the one that comes before
            (which we'll call state t-1), and will calculate the current state (t)
            """
            if i == 1 : 
                # for the firts iteration, take the memory base to be the last base
                # of the first kmer
                memory_base = first_kmer[-1] 
                # Calculating the first transition
                S4 = S0 @ M
                # adding the memory base to the next state
                S4.index = [memory_base + S4.index[b] for b in range(0,len(S4.index))]
                # assigning the next state to be the current state
                Sn[S4.index] = S4

            else:
                bases = ['A','T','C','G']
                # Crearing a temporary vector that will store the current state (t) (shape 16*1)
                # since we are using the last state for all iterations through the bases
                Sni = pd.Series(np.zeros([M.shape[0]]),
                              index = M.index)
                # iterating through the four bases
                for base in bases:
                    # setting the memory base
                    memory_base = base
                    # Creaing a temporary vector that will contatin only the relevant enteries from state t-1
                    # which are all the enteries that ends with the base that is currently iterating (shape 16*1)
                    Sn2 = pd.Series(np.zeros([M.shape[0]]),
                              index = M.index)
                    # Taking all the relevant enteries form state t-1(ends with the current base) (shape 4*1)
                    Sn_base = Sn[Sn.index.str.endswith(base)] 
                    # assign the relevant enteries to the temporary vector Sn2
                    Sn2[Sn_base[Sn_base.index.str.endswith(base)].index] = Sn_base
                    # multiply the temporary vector Sn2 and the transition matrix
                    S4 = Sn2 @ M
                    S4.index = [memory_base + S4.index[b] for b in range(0,len(S4.index))]
                    # Now, S4 contains the transitions from the current base to all other bases
                    # storing the next state of transition from the current base to all the other,
                    # in the next state storing temporary vector.
                    Sni[S4[S4.index.str.startswith(base)].index] = S4[S4.index.str.startswith(base)]
                # assigning the next state to be the current state
                Sn = Sni
            # change the state's name to its iteration
            Sn.rename(i, inplace = True)
            # add the new state to states_df
            states_df = pd.concat([states_df, Sn], axis = 1)
        # states_df = states_df + (1*10**-6)
        return states_df


def get_pMHpos(TM: pd.DataFrame, state_df: pd.DataFrame,
                motif: str, motifs_d:int):
    # print('==================')
    # print(motifs_d, motif, motif[:2])
    # probability for the 1st dimer
    if len(motif) == 1:
        mot_starting_prob = state_df.loc[state_df.index.str.endswith(motif),(motifs_d-1)].sum()
    elif len(motif) == 2:
        mot_starting_prob = state_df.loc[motif[:2], motifs_d] #* state_df.loc[state_df.index.str.endswith(motif),(motifs_d-1)].sum()
        # print(motif[:2], (motifs_d),mot_starting_prob)
    elif len(motif) > 2:
        mot_starting_prob = state_df.loc[motif[:2], motifs_d]
        for b in range(len(motif)-2):
            # calculating the product of P(Bn|Bn-1,Bn-2)
            # print(f'P({motif[b+2:b+3]}|{motif[b:b+2]})')
            # print(TM.loc[motif[b:b+2],motif[b+2:b+3]])
            mot_starting_prob = mot_starting_prob * TM.loc[motif[b:b+2],motif[b+2:b+3]]
    
    return mot_starting_prob

def get_LK(motif_prob_per_state: pd.Series) -> float:
    """
    Calculating liklyhood of finding a motif in a specific distance 
    from the 1st motif.
    Formula:
    let Px(n) be the probability to find the motif in a specific state 
    (location) such that:
        Px(n) =  P(motif's last kmer at S(n))*P(motif_prob) 
        =>
        Liklyhood = ∏ [(1-Px(0))*...*(1-Px(n-1))]*(Px(n))
    Args:
        motif_prob_per_state (np.array): motif prob per state
    Returns:
        The likelyhood of finding a motif in a given location 
        and context (float)
    """
    Nminus1_last_states_prod = 1 - np.array(motif_prob_per_state)[-1]
    return (np.prod(motif_prob_per_state[:-1]) * Nminus1_last_states_prod)

def get_p_value(motif_prob_per_state: np.array, 
                motifs_d: int, early_stop: float) -> float:
    """
    Calculating p-value of finding a motif in a specific distance 
    from the 1st motif.
    Formula:
    Let Ln be the likelihood to find the 2nd motif in a given distance
    from the 1st motif, then:
    p_val = ΣLn, (1<n<distance between the two motifs)
    Args:
        motif_prob_per_state (np.array): motif prob per state
        motifs_d (int):  the distance between the two motifs
        early_stop (float): a value that if passed, the loop
            will break
    Returns:
        the p-value of finding a motif in a specific distance 
        from the 1st motif.
    """
    p_val = 0
    for i in range((motifs_d-1), 1, -1):
        p_val = p_val + get_LK(
            motif_prob_per_state=motif_prob_per_state[:i])
        # if p_val > early_stop: break
    return p_val


def main_markovian_process(sequence: str , motif: str,
                    motifs_d: int, memory_dimer: str,
                    early_stop: float) -> float: 
    # print('----------------------------------------------')
    # print(indel_seq, motif, memory_dimer, motifs_d)
    
    """
    A function that generates all information needed to 
    calculate the likelihood and p-value of a given motif
    at a given state and context 
    Args:
        sequence (str):  
        motif (str): the motif's sequence
        motifs_d (int):  the distance between the two motifs 
            (used to call is indel_len in the previous version)
        memory_dimer (str): the last dimer that comes befor the DSB
    Returns:
        LK (float): the likelihood to find the motif in a specific
            state (corespond to location)
        p_value (float): the p-value of find the motif in a specific
            state (corespond to location)
    """
    ### Inputs Dtype chacking ###
    # type_checker(motif_probabily_calc.__annotations__, locals(), mode='func')

    memory = 2
    motifs_d = motifs_d+memory
    M = get_TM(sequence, memory)
    # calling get_markov_states to create a df with all
    # the transition probability in the given state range
    state_df = get_markov_states(M, first_kmer = memory_dimer,
                                number_of_states = (motifs_d+1))

    pMHpos = get_pMHpos(TM=M, state_df=state_df,
                motif=motif, motifs_d=motifs_d)


    pMHpos_before_obs_pos = np.absolute(1- np.array([get_pMHpos(TM=M, state_df=state_df,
                motif=motif, motifs_d=i) for i in range(2, motifs_d+1)]))

    LK = get_LK(motif_prob_per_state=pMHpos_before_obs_pos)
    p_value = get_p_value(motif_prob_per_state=pMHpos_before_obs_pos,
                     motifs_d=motifs_d, early_stop=early_stop)
    # Calculating the freq_motif_eq
    eq_state = get_equilibrium_state(sequence=sequence, n_mer=2)
    eq_state_df = get_markov_states(M, first_kmer=memory_dimer,
                                number_of_states=3, 
                                eq_state=eq_state,
                                eq=True)
    eq_state_df.rename(columns={'dimers_freq':0}, inplace=True)
    freq_motif_eq = get_pMHpos(TM=M, state_df=eq_state_df,
                motif=motif, motifs_d=2)
    
    return np.array((LK,p_value, pMHpos, freq_motif_eq))  

