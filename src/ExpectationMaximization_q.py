"""
Module implementing an indel-length based Expectation Maximization algorithm to estimate the proportion
of different indel-inducing molecular mechanisms
"""
import regex as re
import pandas as pd
import numpy as np
import ast

#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import seaborn as sns


from scipy.signal import savgol_filter

pd.set_option('display.max_columns', 200)

class EMq:
    def __init__(self, data: pd.DataFrame, 
                initial_theta: float,
                convergence_threshold: float,
                window_size: int, indel_length_distribution_type: str, MH_lengths: list) -> float: #, 
                #MM_lk: float) -> float:
        """
        This class encapsulate all steps that is needed in order
        to perform an Expectation Maximization algorithm
        Args:
            data (pd.DataFrame): the data to operate on
            mechanism_prob (str): a column name that corespond to the probability
                (of any kind) that represents an indel
            initial_theta (float): initial proportion of the mechanism
            convergence_threshold (float): a cutoff to break the main loop when 
                riched a good enougth convergence
        Return:
            The value that the EM gor converged to (float)
        """
        self.df = data
        self.initial_theta = initial_theta
        self.convergence_threshold = convergence_threshold
        self.indel_length_distribution_type = indel_length_distribution_type
        self.window_size = window_size
        #self.MM_lk = MM_lk
        self.log = pd.DataFrame(columns=['MMEJ_theta', 'NHEJ_theta', 'minus_log_likelihood'])
        #self.indel_length_dist_log = pd.DataFrame(columns=['MMEJ_indel_len_dist', 'NHEJ_indel_len_dist'])
        # exeption handling
        assert (self.df['indel_len'].max() < window_size), f"Window size must be at least the window size that was used for RMdetector.pf, current window size={window_size}"
        self.EM_main_loop(theta_a=self.initial_theta)

    def plot_indel_len_dist(self,theta_a: float):
        plt.plot(self.df['r_nMMEJ'], self.df['r_nNHEJ'], 'ko', alpha=0.2)
        plt.plot(theta_a, 0.02, 'ro')
        plt.plot(0.02, (1-theta_a), 'go')
        plt.xlim([0,1])
        plt.ylim([0,1])
        plt.title(f'Iteration number = {self.iter}')
        plt.xlabel('MMEJ')
        plt.ylabel('NHEJ')
        plt.show()


    def update_log(self, theta_a: float) -> None:
        LL = self.df['r_nMMEJ']*theta_a + self.df['r_nNHEJ']*(1-theta_a)
        LL = np.log(LL)
        LL = LL.sum()
        LL = LL * (-1)
        self.log.loc[self.iter, 'MMEJ_theta'] = theta_a
        self.log.loc[self.iter, 'NHEJ_theta'] = 1 - theta_a
        self.log.loc[self.iter, 'minus_log_likelihood'] = LL
        if hasattr(self, 'indel_length_dist_log'):
            self.indel_length_dist_log=pd.concat([self.indel_length_dist_log,pd.DataFrame(["MMEJ"]+self.indel_len_dist_mmej.tolist()).T],ignore_index=True)
            self.indel_length_dist_log=pd.concat([self.indel_length_dist_log,pd.DataFrame(["NHEJ"]+self.indel_len_dist_nhej.tolist()).T],ignore_index=True)
        else:
            self.indel_length_dist_log=pd.DataFrame(["MMEJ"]+self.indel_len_dist_mmej.tolist()).T
            self.indel_length_dist_log=pd.concat([self.indel_length_dist_log,pd.DataFrame(["NHEJ"]+self.indel_len_dist_nhej.tolist()).T],ignore_index=True)

    def get_motif_pos(self,context_seq:str, 
                    motif: str) -> float:
        """
        Args:
            contex_seq (str): the context sequence
            motif (str): the motif sequence
            window_size (int): last position to look for the motif (we set 50)
            
        Returns:
            mmej_pos (np.array): a normalized vector with all position of motifs
        """
        context_reg = context_seq[(round(len(context_seq)/2)):round(len(context_seq)/2+self.window_size)]
        #context_reg = context_seq[(round(len(context_seq)/2)):round((len(context_seq)/2)+self.window_size)]
        if type(motif) == str:
            mmej_pos = np.zeros(self.window_size)
            mmej_motif_pos = np.array([m.end() for m in 
                re.finditer(motif, context_reg, overlapped=True)])
                
            mmej_motif_pos = mmej_motif_pos[mmej_motif_pos<self.window_size]
            # if the motif isnt in context_reg, return a 0 vector 
            if len(mmej_motif_pos) == 0:
                return np.zeros([self.window_size], dtype='int')
            mmej_pos[mmej_motif_pos] = 1
            indel_pos = mmej_pos.nonzero()[0] #returns indices of nonzero elements for dimensions numpy array (first with [0])
            return np.array(indel_pos)
        else:
            return np.ones([self.window_size], dtype='int')
    
    def get_conditional_prob_mechanisms(self, motif: str,obs_len: int, indel_pos: list) -> list: #L: float
        """
        L here is the likelihood that we used to weigh the q
        """
#        if type(motif) == str:
        if len(motif)>0:
            obs_len_freq = np.take(a=self.indel_len_dist_mmej, indices=obs_len)
            indel_len_sum = np.take(a=self.indel_len_dist_mmej,indices=indel_pos).sum()
            p_MMEJ = (obs_len_freq/indel_len_sum) #*L #
        else:
            p_MMEJ = 0
        p_NHEJ = np.take(a = self.indel_len_dist_nhej, indices=obs_len)/self.indel_len_dist_nhej.sum()
        #p_MMEJ = p_MMEJ/(p_MMEJ+p_NHEJ); pNHEJ=1-pMMEJ
        return [p_MMEJ, p_NHEJ]

    def get_indel_length_dist(self, mechanism :str) -> float:
        """
        ADD DOCS HERE
        """
        indel_len_dist = np.array([(self.df.loc[(self.df['indel_len'] == i), mechanism].sum() /
                    self.df.loc[:, mechanism].sum()) for i in range(1,(self.window_size+1))], dtype="float")
        # adding a bias for regularization
        indel_len_dist = indel_len_dist + (self.dist_bias)
        #indel_len_dist = np.array([(i/indel_len_dist.sum()) for i in indel_len_dist])
        return indel_len_dist

    def E_step(self, theta_a: float)-> pd.Series:
        """
        Expectation step
        """
        # update P(indel l | MMEJ) and P(indel l | NHEJ)
        self.df.loc[: ,['r_nMMEJ', 'r_nNHEJ']] = np.array(
            self.df.apply(lambda x: self.get_conditional_prob_mechanisms(motif=x['del_mmej_cand'][0], 
                    indel_pos=x['motif_position'], obs_len=(x['indel_len']-1)), #, L=x[self.MM_lk]), 
                    result_type='expand', axis=1))
        #print(self.df['r_nMMEJ'])
        #print(self.df['r_nNHEJ'])

        self.update_log(theta_a=theta_a)
        # Realigments normalizing
        # calculating realignment wigths (realignment_w)
        self.df['realignment_w'] = ((self.df['r_nMMEJ']*theta_a)+(self.df['r_nNHEJ'])*(1-theta_a))
        self.df['realignment_w'] = self.df.groupby('variant_id',sort=False).apply(
                lambda x: x['realignment_w']/x['realignment_w'].sum()).reset_index()['realignment_w']
        # Multiply by theta_a and normalization
        self.df['r_nMMEJ']=(self.df['r_nMMEJ']*theta_a)/((self.df['r_nMMEJ']*theta_a)+(self.df['r_nNHEJ'])*(1-theta_a))
        self.df['r_nNHEJ'] = 1 - self.df['r_nMMEJ']
        
        # Multiplying realignment wigths and Likelihoods
        self.df['r_nMMEJ'] = self.df['r_nMMEJ'] * self.df['realignment_w']
        self.df['r_nNHEJ'] = self.df['r_nNHEJ'] * self.df['realignment_w']
        
    def M_step(self) -> float:
        """
        A function that performs the Maximization step when maximizing by indel-length.
        """
        # maximize indel length distributions
        self.indel_len_dist_mmej = self.get_indel_length_dist(mechanism='r_nMMEJ')
        self.indel_len_dist_nhej = self.get_indel_length_dist(mechanism='r_nNHEJ')
        if self.indel_length_distribution_type == "uniform":
                self.indel_len_dist_mmej.fill(1/len(self.indel_len_dist_mmej))
                self.indel_len_dist_nhej.fill(1/len(self.indel_len_dist_nhej))
                self.indel_len_dist_mmej[0]=0
                self.indel_len_dist_nhej[0]=0
        elif self.indel_length_distribution_type == "savitzky_golay":
                #https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-for-a-dataset
                nhej_sg=savgol_filter(self.indel_len_dist_nhej, 50, 1) 
                mmej_sg=savgol_filter(self.indel_len_dist_mmej, 50, 1) 
                self.indel_len_dist_nhej=nhej_sg
                self.indel_len_dist_mmej=mmej_sg
                self.indel_len_dist_mmej[0]=0
                self.indel_len_dist_nhej[0]=0
        self.indel_len_dist_nhej=self.indel_len_dist_nhej/np.sum(self.indel_len_dist_nhej)
        self.indel_len_dist_mmej=self.indel_len_dist_mmej/np.sum(self.indel_len_dist_mmej)

        # maximize the proportion of MMEJ 
        theta_a_t = self.df.loc[:,'r_nMMEJ'].sum()/self.n_variants
        return theta_a_t

    def EM_main_loop(self, theta_a: float) -> None:
        """
        The core loop of the EM
        Args:
            theta_a (float): proportion of the mechanism
        Mutate:
            theta_a (float): the value that the EM was
                converged to (represents the converged proportion 
                    of the mechanism).
            posterior_decoding (pd.Series): the posterior decoding 
                using theta_a       
        """
        self.n_variants = len(np.unique(np.array(self.df['variant_id'])))
        self.dist_bias = 1*10**-4

        #---- Initialize the distributions of indel lengths for MMEJ and NHEJ  ------------------------
        mmej_init_dist = np.array([(len(self.df.loc[self.df['indel_len'] == i, 'variant_id'].unique())/
                            self.n_variants) for i in range(1,(self.window_size+1))])
        mmej_init_dist = mmej_init_dist + (self.dist_bias)
        self.indel_len_dist_mmej = np.array([(i/mmej_init_dist.sum()) for i in mmej_init_dist])
        
        nhej_init_dist = np.array([(len(self.df.loc[((self.df['indel_len'] == i) & 
                                                (self.df['del_mmej_cand'].isna() == True)), 'variant_id'].unique())
                                        /(len(self.df.loc[(self.df['del_mmej_cand'].isna() == True), 'variant_id'].unique())+1)) for i in range(1,(self.window_size+1))])
        nhej_init_dist = nhej_init_dist + (self.dist_bias)
        self.indel_len_dist_nhej = np.array([(i/nhej_init_dist.sum()) for i in nhej_init_dist])
        #print(sum(self.indel_len_dist_mmej)) #fabri: they look a bit weird. first element is very small. they sum to 1.
        #print(sum(self.indel_len_dist_nhej))

        # In cases were there are 0 NHEJ, the distribution takes 0.9 of the values of the MMEJ
        if all(nhej_init_dist) == 0:
            nhej_init_dist = np.array(mmej_init_dist) * 0.9
            nhej_init_dist = np.array([i/nhej_init_dist.sum() for i in nhej_init_dist])

        
        #---- Initialize the vectors of motif positions for each indel  ------------------------
        #self.df.loc[:, 'motif_position'] = self.df.apply(lambda x: self.get_motif_pos(motif=x['del_mmej_cand'], 
        #            context_seq=x[f'ref_context_seq_{self.window_size}bp']), axis=1)

        print("fetch motif positions")
        imhl=0
        self.df.loc[:,'motif_position']=self.df['del_mmej_motif_pos'].apply(lambda x: x[0])
        self.df['motif_position'] = self.df['motif_position'].apply(lambda x: '-9999' if x=='' else x ) 
        self.df.loc[:,'motif_position']=self.df['motif_position'].str.split(",") #.values.tolist()
        self.df['motif_position'] = self.df['motif_position'].apply(lambda x: np.array(x, dtype=np.int))
        #print(self.df.loc[:,'motif_position'])
        
        self.df['motif_position'] = self.df['motif_position'].apply(lambda x: [num for num in x if num >= 0])
        
        print("Start EM loop")
        not_converged=True
        self.iter = 0
        while not_converged:
            # if (self.iter%3) == 0:
            #     self.plot_indel_len_dist(theta_a=theta_a)
            self.E_step(theta_a=theta_a) # E-step           
            theta_a_t = self.M_step() # M-step
            distance_from_convergence=(theta_a_t-theta_a) ** 2
            theta_a=theta_a_t
            self.iter +=1
            if distance_from_convergence < self.convergence_threshold:
                not_converged=False
        
        self.del_MMEJ_post_decoding = self.df.loc[:, 'r_nMMEJ']
        self.del_NHEJ_post_decoding = self.df.loc[:, 'r_nNHEJ']
        self.theta_a = theta_a


def sortdesc_list_zerofirst(lst):
    """
    sort a list in descending order but keeping 0 first
    """
    zeros = [x for x in lst if x == 0]
    non_zeros = [x for x in lst if x != 0]
    sorted_non_zeros = sorted(non_zeros, reverse=True)
    sorted_list = zeros + sorted_non_zeros
    return sorted_list
