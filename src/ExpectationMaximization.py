"""
This module is an impementation of the EM algorithm
"""
import sys

import pandas as pd
import numpy as np
pd.set_option('display.max_columns', 200)

class ExpectationMaximization:
    def __init__(self, data: pd.DataFrame, 
                mechanism_prob: str,  
                initial_theta: float,
                convergence_threshold: float) -> float:
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
        self.mechanism_prob = mechanism_prob
        self.initial_theta = initial_theta
        self.convergence_threshold = convergence_threshold
        self.EM_main_loop(theta_a=self.initial_theta)
        
    def E_step(self, theta_a: float)-> pd.Series:
        """
        A function that performes the E-step as follow:
        (P(being a repair mechanism) * mechanism's theta) / 
        [(P(being a repair mechanism) * mechanism's theta) + 
        ((1-P(being a repair mechanism)) * (1-mechanism's theta))] 
        Args:
            theta_a (float): proportion of the mechanism
        Returns:
            A vector with the E-step calculated per element (per indel)
        """
        self.df['Estep']=self.df[self.mechanism_prob]*theta_a+(1-self.df[self.mechanism_prob])*(1-theta_a)
        Estep=self.df.groupby('variant_id',sort=False).apply(lambda x: x['Estep']/sum(x['Estep']))
        return Estep.reset_index()['Estep']
    
    def create_variant_id(self, chrom, original_pos):
        """Define name (id) for a set of potential realignments of the same original indel
        """
        return(chrom+"_"+str(original_pos).split(".")[0])

    def M_step(self) -> float:
        """
        A function that performes the M-step: maximaizing the theta_a and the weights of the realignments
        Rerurns:
            The Avg value of the E-step vector (float)
        """
        n_variants = len(np.unique(np.array(self.df['variant_id'])))
        theta_a_t=sum(self.df['Estep']*self.df[self.mechanism_prob])/n_variants
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
        self.df['variant_id'] = self.df.apply(lambda x: self.create_variant_id(x['CHR'], x['original_pos']), axis=1)
        not_converged=True

        while not_converged:
            self.df['Estep'] = self.E_step(theta_a=theta_a) # E-step
            theta_a_t = self.M_step() # M-step
            distance_from_convergence=(theta_a_t-theta_a) ** 2
            theta_a=theta_a_t
            if distance_from_convergence < self.convergence_threshold:
                not_converged=False
        
        self.posterior_decoding = self.E_step(theta_a=theta_a)
        self.theta_a = theta_a
        print(f'### Old EM Final theta = {self.theta_a}, ({self.mechanism_prob})')

