# 06.02.2022
# importing libraries
import re
import itertools

import pandas as pd
import numpy as np

# from VcfProcess_module import VcfProcess
# from MicroHomology_module_v3 import MicroHomology
from helper import type_checker
import typing



def nucleotides_combination_generator(n_mer: int):
        """
        This function will take the length of the kmer as an input
        and will generate a list of all possible combinations in k
        length, composed from the four nucleotides.
        """
        nuc = ['A','T','C','G']
        # creating an empty list
        n_mer_list = []
        # creating all possible kmers
        for output in itertools.product(nuc, repeat= n_mer):
            n_mer_list.append(''.join(output))
        n_mer_list = sorted(n_mer_list, key=lambda x: x[-1])
        return n_mer_list
    

def transition_matrix_generator(_sequance, k_size):
        """
        This function will take the sequance and kmer length as an input.
        By calling nucleotides_combination_generator() it will generate all
        possible kmers with length==k. Then it will calculate the freq of
        each kmer and produce a dataframe containing all the possible kmers
        and thier freq in a given sequence.
        Example of principle:
        desiered kmer len = 2 , 'AT' followed by 'GG'-> calculate freq of 'ATGG'
        -> the freq is equal to the probabilty of having 'AT' followed by 'GG'
        -> assingning it to the 'AT' column and 'GG' row.
        """
        # making sure that arguments are in the right format
        sequance = _sequance.upper()
        k_size = int(k_size)
        # memmory allocating : vector in length of the sequance, for kmers
        kmers = np.zeros(len(sequance), dtype = 'object')
        # spliting the sequance into kmers
        for i in range(0, (len(sequance) - 1)):
            if i>(len(sequance)-(k_size + 1)): break
            # generating a kmer that is twice the length of the desiered kmer
            else: kmers[i] = sequance[i:(i+(k_size + 1))]
        # calculating the freq of each twice long kmer in the sequence
        kmers_freq = (pd.Series(kmers).value_counts().drop(labels=[0]))#/
#                     len(kmers)).sort_index(key=lambda x : x.str.lower())
#         print(kmers_freq)
       
        # memmory allocating: DataFrame with all possible kmers as rows and columns
        kmers_freq_table = pd.DataFrame(index= nucleotides_combination_generator(k_size),
                    columns= nucleotides_combination_generator(k_size-1))
        
        # assinging the twice-long-kmer freq to the desiered kmer.
        for i in range(0,len(kmers_freq)):
            value = kmers_freq[i]
            if 'N' in kmers_freq.index[i]:
                pass
            else:
                row = kmers_freq.index[i][:k_size]
                col = kmers_freq.index[i][k_size:]
                kmers_freq_table.loc[row,col] = value

#         print(kmers_freq_table)
#         kmers_freq_table.fillna((10**(-6)), inplace = True)
        kmers_freq_table.fillna(0, inplace = True)
#         print(kmers_freq_table)
        # adding 1 to the whole matrix
        kmers_freq_table = kmers_freq_table + 10**(-6)
#         print(kmers_freq_table)
        # normalizing by column (each column sums to 1)
        norm_kmers_freq_table = kmers_freq_table.apply(lambda x: x/x.sum(), axis=1)
        
#         print(norm_kmers_freq_table)
        
        #### This section is to check that the normalization of norm_kmers_freq_table
        # was done rigth. Basicaly the probability of moving from any kmer (row, say AG)
        # to any nucleotide should sum up to 1 (that is, sum(AG->A, AG->T, AG->C, AG->G))
        # if this condition does not hold, the following assert jump:
        ## To check if this assert works, activate the following 1 line: ##
        # norm_kmers_freq_table.iloc[0,0] = 10**6

        norm_kmers_freq_table_row_sum = norm_kmers_freq_table.sum(axis=1).round(decimals=1)
        norm_kmers_freq_table_row_sum = (norm_kmers_freq_table_row_sum != 1.00)      
        assert any(norm_kmers_freq_table_row_sum) != 1, f"""Transition matrix normalization problem
        The following row does not sums to 1.0
        {norm_kmers_freq_table_row_sum}
        Transition Matrix:
        {norm_kmers_freq_table}"""
        #####
        
#         print(f'Summing each state, all rows should sum to 1 \n{norm_kmers_freq_table.apply(sum, axis = 1)}')
        return norm_kmers_freq_table
  
    
# %%timeit
def markov_states_calculator(transition_matrix, _sequance: str, 
                            _first_kmer: str,
                            _number_of_states: int):
        """
        markov_states_calculator will take a fasta context, motif and a
        number of states. The function will construct a 2nd order markov chain
        and produce a dataframe with n states of the chain.
        """
        # seting the memory of the model
        memory = 2
        first_kmer = _first_kmer.upper()
        fasta_context_seq = _sequance.upper()
        # seting the first kmer to start with
        first_kmer = first_kmer[-memory:]
    #     print(first_kmer)
        # Creating the transition matrix -> M, by calling transition_matrix_generator (shape 16*4)
    #     M = transition_matrix_generator(fasta_context_seq, memory)
        M = transition_matrix
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
        for i in range(1,_number_of_states):
            """
            The idea here is that in each iteration, we look at the one that comes before
            (which we'll call state t-1), and will calculate the current state (t)
            """
        #     print(f'####### iteration num:{i} #######')
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
    #             print(S4,S4.sum(), memory_base,'\n',
    #                  S0)
        #         print(Sn, Sn.sum())

            else:
                bases = ['A','T','C','G']
                # Crearing a temporary vector that will stor the current state (t) (shape 16*1)
                # since we are using the last state for all iterations through the bases
                Sni = pd.Series(np.zeros([M.shape[0]]),
                              index = M.index)
                # iterating through the four bases
                for base in bases:
                    # setting the memory base
                    memory_base = base
        #             print('##', memory_base)
                    # Creaing a temporary vector that will contatin only the relevant enteries from state t-1
                    # which are all the enteries that ends with the base that is currently iterating (shape 16*1)
                    Sn2 = pd.Series(np.zeros([M.shape[0]]),
                              index = M.index)
                    # Taking all the relevant enteries form state t-1(ends with the current base) (shape 4*1)
                    Sn_base = Sn[Sn.index.str.endswith(base)] 
                    # assign the relevant enteries to the temporary vector Sn2
                    Sn2[Sn_base[Sn_base.index.str.endswith(base)].index] = Sn_base
        #             print('Sn_base',Sn_base, Sn_base[Sn_base.index.str.endswith(base)].index[0])
        #             print('Sn2- the relevant part of the last state',Sn2)
                    # multiply the temporary vector Sn2 and the transition matrix
                    S4 = Sn2 @ M
                    S4.index = [memory_base + S4.index[b] for b in range(0,len(S4.index))]
                    # Now, S4 contains the transitions from the current base to all other bases
        #             print('#### S4- the next state',S4, S4.sum())
                    # storing the next state of transition from the current base to all the other,
                    # in the next state storing temporary vector.
                    Sni[S4[S4.index.str.startswith(base)].index] = S4[S4.index.str.startswith(base)]
                # assigning the next state to be the current state
                Sn = Sni
            # change the state's name to its iteration
            Sn.rename(i, inplace = True)
            # add the new state to states_df
            states_df = pd.concat([states_df, Sn], axis = 1)
        states_df = states_df + (1*10**-6) # 06.01.2022, Ask Fabrizio if this is done in the right way
        return states_df


def motif_probabily_calc(_sequance: str , 
                                 _motif: str,
                                _indel_len: int,
                        _memory_dimer: str):
        """
        This function will call markov_states_calculator in order to calculate
        all states probabilities, then it'll caculate the probability to find
        a given motif within a given range.
        """
        ### Inputs Dtype chacking ###
        # type_checker(motif_probabily_calc.__annotations__, locals(), mode='func')

        
        memory = 2
        M = transition_matrix_generator(_sequance, memory)
#         print(M)
        # if _memory_dimer == None: _memory_dimer = _motif
        # _motif = _motif.upper()
        # _memory_dimer = _memory_dimer.upper()

        # calling markov_states_calculator to create a df with all
        # the transition probability in the given state range
        state_df = markov_states_calculator(M,_sequance = _sequance , 
                                     _first_kmer = _memory_dimer,
                                    _number_of_states = _indel_len)

        """
        This if statment is to validate that the calculation of
        state_df was done rigth. 
        Basically what we chack is that the sum of all state n 
        probabilities that ends with certain nucleotide is equal 
        to the sum of the probabilities that starts with the very
        same nucleotide at state n+1. 
        """
        if _indel_len > 1:
            nuc = _memory_dimer[-1] # plug here base which you want to check
            n = _indel_len -2  # plug here a state that you want to compare with the following one
#             stateN = state_df.iloc[:,n]
            # Summing up the probabilites at state n that ends with nuc
            stateN_sum = round(
                state_df.iloc[:,n][
                    state_df.iloc[:,n].index.str.endswith(nuc)].sum(), 6)
#             stateN_plu1 = state_df.iloc[:,(n+1)]
            # Summing up the probabilites at state n + 1 that starts with nuc
            stateN_plu1_sum = round(
                state_df.iloc[:,(n+1)][
                    state_df.iloc[:,(n+1)].index.str.startswith(nuc)].sum(), 6)

            assert stateN_sum == stateN_plu1_sum, f"""
            State matrix probem (cheked nucleotide: {nuc}): sum of state {n} :
            (Sum = {stateN_plu1_sum})
            does not match sum of state {n+1} :
            (Sum = {stateN_plu1_sum})
            {state_df}"""
           
        # print(state_df)
        """
        This section is to calculate the probability to find the the mmej motif
        given the sequence, unrelated to the current state (analog to the position
        in the range in which we are calculating the probabily for)
        This is done by calculating the product of P(N|N-1,N-2) from the 
        transition matrix.
        """
        # Creating a list with the first two (N-1,N-2) dimers since they are in a
        # negative position to the mmej motif, they need to calc seperatly
        motif_dimers = [_motif[-2:], 
                   (_motif[-1:] + _motif[:1])]
#         print(motif_dimers)
        # adding the rest (N-1,N-2) dimers, this list is a list of (the N-1,N-2 in P(N|N-1,N-2)).
        [motif_dimers.append(_motif[i : i+2]) for i in range(0,len(_motif)-1)]
#         print(motif_dimers)
        # for motifs with length < memmory (len(mmej_motif) < 2), one can't set
        # the last motif dimer as the memmory dimer (the N-1,N-2 in P(N|N-1,N-2)).
        # therefore one should provid it manualy whene calling the function.
        # this is _memory_dimer.
        if len(_motif) == 1 : motif_dimers = [_memory_dimer]
#         print(motif_dimers)
        # setting initial probability to 1
        motif_prob = 1

        # iterating through both N (bases in the motif) and (N-1,N-2) 
        # (from motif_dimers)
        for b,dimer in zip(range(0,len(_motif)),motif_dimers):
            # calculating the product of P(N|N-1,N-2)
#             print(dimer,_motif)
#             print(M.loc[dimer,_motif[b:b+1]])
            motif_prob = motif_prob * M.loc[dimer,_motif[b:b+1]]
            # This is to check manualy that calculation was done right
        
        #### if one wants to look at the calculations, activate the next two lines:
        #     print(f'P({_motif[b:b+1]} | {dimer}) = {M.loc[dimer,_motif[b:b+1]]}, motif probability = {motif_prob}')
        # print(f'Transition matrix:\n{M} \nState df:{state_df}') 
        ####


        """
        Calculating the probability to find the motif at the range(S0,SN)
        as : 1 - [product of [1-(P(motif's last kmer at SN) * P(motif_prob))]]
        """
        # getting P(motif's last kmer at SN) from the states matrix
        state_prob = [state_df.loc[motif_dimers[0],i] for i in state_df.columns]

        # calculating 1-(P(motif's last kmer at SN) * P(motif_prob))
        state_motif_prob = np.array([(1-(state_prob[s] * motif_prob)) for s in range(0,len(state_prob))])

        # calculating the total probability to find the motif at some state 
        # in a range(S0,SN)
        total_motif_prob = 1 - np.prod(state_motif_prob)
#         print(state_motif_prob,'\n',state_prob)
#         print(total_motif_prob)
        
        return total_motif_prob




# """
# This part is to check that calculations from markov_states_calculator done rigth
# """
# # %%timeit
# # for simplicty, I'll call the function for only 5 states
# # output = markov_states_calculator(transition_matrix_generator(fasta_context_seq_2K, 2),
# #                                   _sequance = fasta_context_seq_2K , 
# #                                  _first_kmer = mmej_cand,
# #                                 _number_of_states = 8)
# # print(f'The function output: \n{output}')
# # print(f'Summing each state, all columns should sum to 1 \n{output.sum()}')

# # ## A check that calculations was done rigth: ##
# # nuc = 'G' # plug here base which you want to check
# # n = 5 # plug here a state that you want to compare with the following one
# # stateN = output.iloc[:,n]
# # print(f'Base: {nuc}, State: {n}, Sum: {stateN[stateN.index.str.endswith(nuc)].sum()}')
# # stateN_plu1 = output.iloc[:,(n+1)]
# # print(f'Base: {nuc}, State: {n+1}, Sum: {stateN_plu1[stateN_plu1.index.str.startswith(nuc)].sum()}')





# %%timeit
"""
Example for mmej motif that is longer then the memory of the model.
In cases like this, the memory dimer is the last dimer of the motif,
So in this example the motif is AGGT, therefore the probability to 
find the motif occur again in distance 0bp from the first occurance,
(i.e. AGGT|AGGT) is the P(AGGT|GT) or:
Ps0(GT)*P(A|GT)*P(G|TA)*P(G|AG)*P(T|GG)
"""
# motif_probabily_calc(_sequance = fasta_context_seq_2K , 
#                                  _motif = mmej_cand,
#                                 _indel_len = indel_length,
#                     _memory_dimer = None)



"""
Example for mmej motif that is shorter then the memory of the model.
In cases like this, the memory dimer is the last dimer that comes rigth
before the indel position, and one should manualy provid it when calling
the function.
So for example if the motif is T, and the last dimer before the indel is AT,
therefore the probability to find the motif occur again in distance 0bp from the
first occurance, (i.e. T|AT) is the P(T|AT) or:
Ps0(AT)*P(T|AT).

the probability to find the motif occur again within a distance of 3bp from the
first occurance, (i.e. T|AT) is by the law of total probabilities
Ps0(AT)*P(T|AT) * Ps1(AT)*P(T|AT) * Ps2(AT)*P(T|AT).
"""
# motif_probabily_calc(_sequance = fasta_context_seq_2K , 
#                                  _motif = 'T',
#                                 _indel_len = 3,
#                     _memory_dimer = 'AT')


