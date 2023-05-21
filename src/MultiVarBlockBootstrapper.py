"""
This module contains all the functions that toghether perform
a block bootstrap on a given dataframe for multiplr variables
"""

# import libreries

from statistics import mean
from typing import Optional, Union
from matplotlib.pyplot import get

import numpy as np
import pandas as pd
from pyparsing import col

from helper import type_checker

def blocks_generator(df: pd.DataFrame, block_size: int, 
            coordinate_col: str) -> pd.DataFrame:
        """
        This function taking the indel's coordinates and round it by the
        interval size in order to create binns along the chromosome that
        represents the position of groups of indels along the chromosome.
        Args:
            self.df, self.coordinate_col, self.block_size :
                docs are at __init__ docs.
        Returns:
            df (pd.DataFrame): a dataframe that contains the block column
        
        """
        df.loc[: ,'block'] = df.loc[: ,coordinate_col]/block_size
        df.loc[: ,'block'] = df.loc[: ,'block'].round(decimals=0).astype('int32')
        return df
        
        
def get_n_blocks(data: pd.DataFrame,
    force_n_block: Optional[int]=None) -> int:
    """
    Gets the maximum block number.
    Args:
        data (pd.DataFrame): a dataframe with the block column
        force_n_block (either None, or int): a flag to force the
            block numbers.
    Returns:
        n_block (int) The block number per data point.
    """
    if force_n_block != None:
            n_block = force_n_block +1
    else:
        n_block =  int(data.loc[:, 'block'].max())
    return n_block

def get_block_stat(data: pd.DataFrame, 
                    func: str, cond: str,
                    col_list: list,
                    n_blocks:int, agg_by:str='block') -> pd.DataFrame:
    """
    Calculate a statistic per block per mechanism.
    Args:
        data (pd.DataFrame): The data itself 
        func (str): a function to calculate the stat by 
        cond (str): a condition to devide the data by 
            (genic for example)
        col_list (list): a list with the mechanisms to bootstrapp
        agg_by (str): a column to aggregete by
        n_blocks (int): number of max blocks
    Returns: 
        agg_data (pd.DataFrame): a dataframe with the statistic
            per block per mechanism per cond.
    """
    # Calculating stats for cond == True
    pos_agg_data = data.loc[
            (data[cond] == True), :].groupby(agg_by).agg(func).copy()
    pos_agg_data[cond] = True
    pos_data = pd.DataFrame(index=range(0,n_blocks),
                            columns=pos_agg_data.columns)
    pos_data[cond] = True
    pos_data.loc[pos_agg_data.index, :] = pos_agg_data.loc[:,:]
    # Calculating stats for cond == False
    neg_agg_data = data.loc[
            (data[cond] == False), :].groupby(agg_by).agg(func).copy()
    neg_agg_data[cond] = False
    neg_data = pd.DataFrame(index=range(0,n_blocks),
                            columns=neg_agg_data.columns)
    neg_data[cond] = False    
    neg_data.loc[neg_agg_data.index, :] = neg_agg_data.loc[:,:]
    
    agg_data = pd.concat([pos_data, neg_data], axis=0)   
    return agg_data.loc[:,[cond] + col_list]
    
def block_stat_df_to_density(data: pd.DataFrame, 
                            norm_factors: list,
                            col_list: list,
                            cond: str) -> pd.DataFrame:
                            
    # normalizing by cond length in the genome
    data.loc[(data[cond] == True), col_list] = (data.loc[
        (data[cond] == True), col_list]/norm_factors[0])
    data.loc[(data[cond] == False), col_list] = (data.loc[
        (data[cond] == False), col_list]/norm_factors[1])
    return data
    
def block_stat_df_norm_to_partial_stat(
                        data: pd.DataFrame,
                        col_list: list,
                        cond: str) -> pd.DataFrame:
    # normalizing per by mechanism relative portion 
    # in the block
    row_sum = data.loc[:, col_list].sum(axis=1)
    row_sum.replace(to_replace=0,value=1, inplace=True)
    for col in col_list:
        data.loc[:, col] = (data.loc[:, col] / row_sum)

    return data
    
def sample_block_num_generator(n_block: int,
                                n_samples: int) -> pd.DataFrame:
        """
        Generates a dataframe with lampled block numbers as index
        """
        blocks = [i for i in range(0,n_block)]
        sample_block_num = np.random.choice(a = blocks , 
            size = (n_block* n_samples), replace = True)
        samp_block_num_df = pd.Series(index = sample_block_num, dtype='float64')
        # print(samp_block_num_df)
        return samp_block_num_df.index.tolist()

def add_samp_number(data: pd.DataFrame,n_block: int) -> pd.Series:
    """
    Adding a column with sample number to data.
    Args:
        data (pd.DataFrame): the bootstrapped data
        n_block (int): number of maximum blocks
    Return:
        sample_n (pd.Series): a series with the sample number.
    """
    s_i, e_i = 0, None
    sample_n = pd.Series([np.nan for i in range(0,int(round(len(data)/2, 0)))])
    for i in range(0,int(round(len(data)/n_block, 0)/2)):
        e_i = s_i + n_block 
        sample_n[s_i:e_i] = int(i)
        s_i = e_i
    sample_n = pd.concat([sample_n,sample_n], axis=0)
    sample_n.reset_index(drop=True, inplace=True)
    
    return sample_n

def get_sampled_data(data: pd.DataFrame,
                    sampled_blocks:list, cond:str,
                    n_block: int) -> pd.DataFrame:
    """
    Gets the sampled blocks according to the list of block numbers
    and associate them with the block's stats.
    Args:
        data (pd.DataFrame): the block statistic data
        sampled_blocks (list): a list with the sampled block numbers
        cond (str): a condition to split the data by
    Returns:
        samp_df (pd.DataFrame): a dataframe with the sampled blocks
            for cond == True and cond == False
    """
    pos_samp_df = data.loc[(data[cond] == True), :]
    neg_samp_df = data.loc[(data[cond] == False), :]
    
    samp_df = pd.concat([pos_samp_df.loc[sampled_blocks,:],
                        neg_samp_df.loc[sampled_blocks,:]], axis=0)
    samp_df.reset_index(inplace=True)
    samp_df.rename(columns={'index': 'block'}, inplace=True)
    samp_df['sample_n'] = add_samp_number(data=samp_df, n_block=n_block)
    return samp_df


def cond_sample_avg_block_stat(data: pd.DataFrame, cond: str, 
                            cols: list) -> pd.DataFrame:
    """
    Calculating a stat per re-sample per mechanism as:
    mean(blocks) per mechanism.
    Args:
        data (pd.DataFrame): data from get_sampled_data
        cond (str): a condition to split the data by
        cols (list): list of relevant columns
    Return:
        (pd.DataFrame)
    """
    samp_n = data.loc[:,'sample_n'].max()
    pos_comp_df = pd.DataFrame(columns=(cols + ['sample_n'] + [cond]))
    neg_comp_df = pd.DataFrame(columns=(cols + ['sample_n'] + [cond]))
    # iterating through the re-samples
    for samp in range(0, int(samp_n+1)):
        # calculating the stat: summing up blocks per mechanism
        # and recording it in pos/nef_comp_df
        pos_df = data.loc[((data[cond] == True) & 
                    (data['sample_n'] == samp)), cols].mean(axis=0)
        # print(pos_df/13)
        pos_comp_df.loc[samp, cols] = pos_df
        pos_comp_df.loc[samp, 'sample_n'] = samp
        pos_comp_df.loc[samp, cond] = True
    
        neg_df = data.loc[((data[cond] == False) & 
                    (data['sample_n'] == samp)), cols].mean(axis=0)
        neg_comp_df.loc[samp, cols] = neg_df
        neg_comp_df.loc[samp, 'sample_n'] = samp
        neg_comp_df.loc[samp, cond] = False
    
    return pd.concat([pos_comp_df,neg_comp_df])

def plot_data_formater(data: pd.DataFrame,
                       cols: list,
                       cond: str) -> pd.DataFrame :
    """
    Creating a dataframe with the right format for ploting
    with seaborn/matplotlib
    Args:
        data (pd.DataFrame): the data from cond_sample_avg_block_stat
        cols (list): list of relevant columns
        cond (str): a name of a column to split the data by
    Returns:
        boot_plt_df (pd.DataFrame): the reformated
            data.
    """
    boot_plt_df = pd.DataFrame()
    for col in cols:
        d = data.loc[:,[col,cond]]
        d['variable'] = col
        d.rename(columns={col: 'value'}, inplace=True)
        boot_plt_df = pd.concat([boot_plt_df,d],axis=0)
        boot_plt_df.reset_index(inplace=True, drop=True)
    return boot_plt_df

def cond_sample_sum_comparison(comp_stat_df: pd.DataFrame, 
                                cond: str, cols: list) -> pd.DataFrame:
    """
    Comparing if the sum of blocks were cond == True is greater
    to the onse were cond == False or not per mechanism (pairwise).
    Args:
        comp_stat_df (pd.DataFrame): the data from cond_sample_avg_block_stat 
        cond (str): a conditions to divide the data by
        cols (list): columns to compair pairwise
    Return:
        comp_df (pd.DataFrame): a list containing the comparision per sample
    """
    comp_df = (comp_stat_df.loc[comp_stat_df[cond] == True, cols] < 
                comp_stat_df.loc[comp_stat_df[cond] == False, cols])

    comp_df.replace(to_replace=[False, True], value=[0,1], inplace=True)
    return comp_df       


def get_p_val(comp_data: pd.DataFrame,
            BC: bool = False) -> pd.Series:
    """
    Calculating p-value of difference in proportion
    per mechanism based on the comparison data.
    Args:
        comp_data (pd.DataFrame): data from cond_sample_sum_comparison
        BC (bool): a flag to activate Bonferroni correction
    Returns:
        (dict): a dict with the p-vals per mechanism
    """      
    def p_val_comp(vec: pd.Series):
        if len(vec[vec == 1]) == 0: return 0
        else: return round(len(vec[vec == 1]) / len(vec), 4)
    
    if BC: 
        BC_p_val= dict(zip(comp_data.columns, 
            [(round(p_val_comp(comp_data[col])*len(comp_data.columns), 4)) 
                        for col in comp_data.columns]))
        
        # setting p values that are > 1 to be 'p-value > 1'
        for p in BC_p_val.keys():
            if BC_p_val[p] > 1: BC_p_val[p] = 'p-value > 1' 

        return BC_p_val
    else: return dict(zip(comp_data.columns, 
            [p_val_comp(comp_data[col]) for col in comp_data.columns]))
    


def get_CI(data: pd.DataFrame, CI: float,
            cols: list):
    """
    Calculates the confidence interals per mechanism
    Args:
        data (pd.DataFrame): data from get_sampled_data
        CI (float): confidence interval
        cols (list): list of the mechanisms
    Return:
        (dict): a dict containing upper and lower CI per
            mechanism
    """
    upper_ci = dict(zip([f'{col}_{CI}' for col in cols],
            [data.loc[:,col].quantile(q=CI) for col in cols]))
    lower_ci = dict(zip([f'{col}_{round((1-CI),3)}' for col in cols],
            [data.loc[:,col].quantile(q=round((1-CI),3)) for col in cols]))
    return {**upper_ci, **lower_ci}


def add_CI_to_avg_data(data: pd.DataFrame, data_to_append: pd.DataFrame, 
                        CI: float, cond: str) -> None:
    """
    Adding a confidence interval columns to a dataframe
    Args:
        data (pd.DataFrame): data to calc CI by
        data_to_append (pd.DataFrame): data to add the CI columns to
        CI (float): CI value
        cond (str): conditions to split the data by
    Returns:
        data_to_append (pd.DataFrame): a dataframe with the additional
            confidence interval
    """
    data_to_append[f'CI{CI}'] = np.nan
    data_to_append[f'CI{round((1-CI),3)}'] = np.nan
    
    for j in data_to_append['variable']:
        data_to_append.loc[((data_to_append[cond] == True) &
             (data_to_append['variable'] == j)),f'CI{CI}']= data.loc[
                 (data[cond] == True),j].quantile(q=CI)
        data_to_append.loc[((data_to_append[cond] == True) & 
            (data_to_append['variable'] == j)),f'CI{round((1-CI),3)}']= data.loc[
                (data[cond] == True),j].quantile(q=round((1-CI),3))
        
        data_to_append.loc[((data_to_append[cond] == False) &
             (data_to_append['variable'] == j)),f'CI{CI}']= data.loc[
                 (data[cond] == False),j].quantile(q=CI)
        data_to_append.loc[((data_to_append[cond] == False) & 
            (data_to_append['variable'] == j)),f'CI{round((1-CI),3)}']= data.loc[
                (data[cond] == False),j].quantile(q=round((1-CI),3))
    return data_to_append

def get_samples_avg(data: pd.DataFrame, col_list: list,
                     cond: str, comment: str = None) -> pd.DataFrame: 
    """
    Averaging the samples per mechanism per cond ==True/False
    Args:
        data (pd.DataFrame): data from cond_sample_avg_block_stat
        col_list (list): columns to average
        cond (str): condition to split the data by
        comment (str): an optional comment to add for convinient 
            (for example a comment on the input data such
             as indel length) 
    Returns:
        avg_df (pd.DataFrame): A dataframe that contains average values
            per mechanism per cond == True/False
    """
    pos_data = data.loc[(data[cond] == True), col_list].copy()
    neg_data = data.loc[(data[cond] == False), col_list].copy()
    pos_avg_df = pos_data.apply(mean, axis=0)
    neg_avg_df = neg_data.apply(mean, axis=0)
    
    avg_df = pd.DataFrame([pos_avg_df, neg_avg_df])
    avg_df[cond] = [True, False]
    if comment != None:
        avg_df['comment'] = comment
    return avg_df


class MultiVarBlockBootstrapper:
    def __init__(self, df: pd.DataFrame, block_size: int, 
            coordinate_col: str,
            n_sample: int, func: str, cond: str,
            col_list: list, 
            norm_factors: list,BC: bool = False,
            force_n_block: Union[int, None] = None,
            comment: str = None) -> None:
        
        self.df = df
        self.block_size = block_size
        self.coordinate_col = coordinate_col
        self.n_sample = n_sample
        self.func = func
        self.cond = cond
        self.col_list = col_list
        self.norm_factors = norm_factors
        self.BC = BC
        self.force_n_block = force_n_block
        self.comment = comment

        # type_checker(self.__dict__, self.__init__.__annotations__,mode='class')
        
        self.BootStrapper_multi_var()

    def BootStrapper_multi_var(self):
        """
        This is the main function that is basically the actuall pipline of the module.
        Args:
            df (pd.DataFrame): the original data 
            block_size (int): the block size in bp 
            coordinate_col (str): a columns name that contains 
                    the indel's coordinate
            n_sample (int): number of re-sampling from the data
            func (str): name of a function to calculate the stat by
            cond (str): a name of a column to split the data by
            col_list (list): a list with column names of the 
                variables to bootstraps 
            norm_factors (list): a list with length of 
                [cond==True, cond==False] that will be use to calc density
            BC (bool): a flag to activate Bonferroni correction
            force_n_block (Optional[int]): a flag to force the block number
        Returns:

        """
                
        #1. create blocks
        # genetating blocks
        data =  blocks_generator(self.df, 
                        self.block_size, self.coordinate_col)
        
        # get block number
        n_blocks = get_n_blocks(data=data, 
                        force_n_block=self.force_n_block)  
        
        #2. calc the block stat
        block_stat_df = get_block_stat(data, self.func, self.cond,
                                self.col_list, n_blocks=n_blocks)
        # print(f'block_stat_df (Stat = {self.func}):\n',block_stat_df)
        #3. Normalizing
        block_stat_density = block_stat_df_to_density(
                        data=block_stat_df, 
                        norm_factors=self.norm_factors,
                        col_list=self.col_list,
                        cond=self.cond)
        # print(f'Density ( that is, {self.func} / (for {self.cond} == True:{self.norm_factors[0]}, for {self.cond} == False:{self.norm_factors[1]}):\n',block_stat_norm_by_length)
        self.block_df = block_stat_df_norm_to_partial_stat(
                                    data=block_stat_density,
                                col_list=self.col_list, cond=self.cond)
        # print('Partial density (That is, mechaisms density / total mechanisms density):\n:' ,self.block_df)
        #4. generate sample block numbers
        self.sampled_blocks = sample_block_num_generator(n_blocks, self.n_sample)
        
        #5. take a sample set
        self.samp_df = get_sampled_data(self.block_df, self.sampled_blocks, 
                        cond=self.cond,n_block=n_blocks)
        add_samp_number(self.samp_df, n_block=n_blocks)
        # print('Sampled data:\n', self.samp_df)

        #6. compare the two groups (p-val,CI)
        self.comp_stat_df = cond_sample_avg_block_stat(data=self.samp_df, 
                                    cond=self.cond, cols=self.col_list)
        # print('Summing the partial density of all blocks per re-sample (to get density along the whole chromosome)\n',self.comp_stat_df)
        comp_df = cond_sample_sum_comparison(comp_stat_df=self.comp_stat_df,
                                    cond=self.cond, cols=self.col_list)
        
        self.boot_plt_df = plot_data_formater(data=self.comp_stat_df,
                        cols=self.col_list, cond=self.cond)
        # print(self.boot_plt_df)
        self.p_val_dict = get_p_val(comp_data=comp_df, BC=self.BC)
        
        self.ci_dict = get_CI(data=self.comp_stat_df, CI=0.95, 
                                            cols=self.col_list)
        
        # print(self.comp_stat_df)
        self.sample_avg_df = get_samples_avg(data=self.comp_stat_df, 
                        col_list=self.col_list, cond=self.cond,
                        comment=self.comment)
        # print(self.sample_avg_df)
        # print(self.ci_dict)

        self.boot_avg_plt_df = plot_data_formater(data=self.sample_avg_df,
                        cols=self.col_list, cond=self.cond)
        
        # Adding CI to boot_avg_plt_df
        self.boot_avg_plt_df = add_CI_to_avg_data(data=self.comp_stat_df, data_to_append= self.boot_avg_plt_df,
                                 CI=0.95, cond=self.cond)
        # print(self.boot_avg_plt_df)