

# import libreries
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from helper import type_checker


class ROCCurve:
    def __init__(self, data: pd.DataFrame, 
                mechanism_prob: str,  
                mechanism_name: str):
        self.df = data
        self.mechanism_prob = mechanism_prob
        self.mechanism_name = mechanism_name
        
        
    
    def confution_matrix_generation(self, cutoff) -> pd.Series:
        """
        Takes a dataframe with probabilities of indels to be del MMEJ,
        then uses a cutoff to True/False positive/negative
        Args: 
            cutoff (float): A cutoff to split T/F 
        Mutate:
            self.mechanism_counts (pd.Series): A series containing:
                      T_pos, F_pos, T_neg, F_neg
        """
        T_pos = len(self.df.loc[((self.df[self.mechanism_prob] > cutoff) & 
                            (self.df['ground_true'] == self.mechanism_name)), :])
        F_pos = len(self.df.loc[((self.df[self.mechanism_prob] >= cutoff) & 
                            (self.df['ground_true'] != self.mechanism_name)), :])
        T_neg = len(self.df.loc[((self.df[self.mechanism_prob] < cutoff) & 
                            (self.df['ground_true'] != self.mechanism_name)), :])
        F_neg = len(self.df.loc[((self.df[self.mechanism_prob] <= cutoff) & 
                            (self.df['ground_true'] == self.mechanism_name)), :])

        mechanism_counts = pd.Series(index=['T_pos', 'F_pos', 'T_neg', 'F_neg'],
                    data=[T_pos, F_pos, T_neg, F_neg])
        self.mechanism_counts = mechanism_counts
        
    def get_T_pos_rate(self) -> float:
        """
        This function calculates the true positive rate (sensitivity)
        of a given data.
        T positive = (T positive)/(T positive + F negative)
        """
        return (self.mechanism_counts['T_pos']/(self.mechanism_counts['T_pos']
                                + self.mechanism_counts['F_neg']))

    def get_F_pos_rate(self) -> float:
        """
        This function calculates the false positive rate (1-sensitivity)
        of a given data.
        F positive = (F positive)/(F positive + T negative)
        """
        return (self.mechanism_counts['F_pos']/(self.mechanism_counts['F_pos']
                                + self.mechanism_counts['T_neg']))

    def get_precition(self) -> float:
        """
        This function calculates the precition of a given data
        precition = (T positive)/(T positive + F positive)
        """
        return (self.mechanism_counts['T_pos']/(self.mechanism_counts['T_pos']
                                + self.mechanism_counts['F_pos']))
    
    def ROC_cuerve_data_generation(self, lower_lim=0, upper_lim=1, 
                                            step=0.1) -> float:
        """
        This function calculates the data points for the ROC curve using 
        cutoffs from a given range and step.
        Args:
            lower_lim (float): lowest cutoff
            upper_lim (float): highest cutoff
        Mutate:
            self.ROC_curve_points (pd.DataFrame): a datafreame thae contains 
                cutoff, T_pos_rate, F_pos_rate and precition values per cutoff
        """
        # Spliting the cutoffs range into 3 parts: first and last parts are
        # exponentialy more dense (i.e. have more cutoffs), the middle part
        # is less dense and uses the defined step
        near_lower_bound_cutoffs =  [i for i in np.linspace(start=lower_lim, 
                    stop=step,num=int(1/(step*0.2)))]
        near_upper_bound_cutoffs = [i for i in np.linspace(start=upper_lim*0.9999, 
                            stop=upper_lim,num=int(1/(step*0.2)))]
        mid_cutoff_ls = [i for i in np.linspace(start=step+0.0001, 
                            stop=upper_lim*0.9998,num=int(1/step))]
        cutoff_ls = near_lower_bound_cutoffs + mid_cutoff_ls + near_upper_bound_cutoffs
        
        # creating an empty dataframe
        self.ROC_curve_points = pd.DataFrame(columns=['cutoff', 'T_pos_rate', 'F_pos_rate'], 
                    index=[i for i in range(len(cutoff_ls))]) 
        
        for inx,co in enumerate(cutoff_ls):
            self.confution_matrix_generation(cutoff=co)
            self.ROC_curve_points.loc[inx,'cutoff'] = co
            self.ROC_curve_points.loc[inx,'T_pos_rate'] = self.get_T_pos_rate()
            self.ROC_curve_points.loc[inx,'F_pos_rate'] = self.get_F_pos_rate()
        
        self.ROC_curve_points.loc[self.ROC_curve_points.index[0],'T_pos_rate'] = 1
        self.ROC_curve_points.loc[self.ROC_curve_points.index[0],'F_pos_rate'] = 1
        
        # pd.set_option('display.max_rows', None)
        # print(self.ROC_curve_points)

        
        
    def get_significant_cutoff(self, FPR: float) -> float:
        """
        gets a false positive rate (FPR) and reruens the lowest 
        cutoff that produce an FPR <= to the input FPR
        Args:
            FPR (float): maximum false positive rate desiered
        Returns:
            FPR_cutoff (float): the lowest cutoff that produce
                the desiered FPR
        """
        above_FPR_ls = self.ROC_curve_points.loc[
            ((self.ROC_curve_points['F_pos_rate'] <= FPR) & 
            (self.ROC_curve_points['T_pos_rate'])>0),'cutoff'].to_list()
        if len(above_FPR_ls)>0: return above_FPR_ls[0]
        else: return 'No satisfing cutoff for given false positive rate'

    def plot_ROC_curve(self, x_variable:str):
        """
        ploting the ROC curve itself based on the data from ROC_cuerve_data_generation
        Args:
            x_variable (str): a name of the parameter that'll be on the x-axis,
                    can be: False positive rate or precition.
        """
        print(self.ROC_curve_points[x_variable], 
                self.ROC_curve_points['T_pos_rate'])
        plt.plot(self.ROC_curve_points[x_variable], 
                self.ROC_curve_points['T_pos_rate'], '.r-')
        plt.plot([0, 1], [0, 1],ls="--", c=".3")
        plt.xlabel(x_variable)
        plt.ylabel("True Positive Rate")
        plt.title(f'{self.mechanism_name} ROC curve')

        plt.show()



