# 03.03.2022
"""
This script contains all the functions that toghether perform
a block bootstrap on a given dataframe
"""

# import libreries
import re
from typing import Optional

import numpy as np
import pandas as pd

from helper import type_checker

class BlockBootstrap:
    def __init__(self, data: pd.DataFrame, _block_size: int, 
            _coordinate_col: str,
            _n_samples: int, force_n_block: Optional[int]=None):
        """
        Args:
            _df (pandas DataFrame): the data
            _coordinate_col (pandas Series): a colume that contains the indel's 
                coordinates (start/end)
            _block_size (int): size of the desiered blocks
            n_samples (int): Number of samples from the blocks
            _cond (str): the column name of the column that contains the
                category we want to check.
            _n_samples (int): number of resempeling.
        """
        
        ### Inputs Dtype chacking ###    
        # assert check_argument_types()
               
        # assining inputs as atributes
        self.df, self.block_size, self.coordinate_col = data, _block_size, _coordinate_col
        self.n_samples = _n_samples
        self.force_n_block = force_n_block
        self.coordinates_to_blocks()
      
    def coordinates_to_blocks(self):
        """
        This function taking the indel's coordinates and round it by the
        interval size in order to create binns along the chromosome that
        represents the position of groups of indels along the chromosome.
        Args:
            self.df, self.coordinate_col, self.block_size :
                docs are at __init__ docs.
        Returns:
            (int) The block number per data point.
        """
        self.df.loc[: ,'block'] = self.df.loc[: ,self.coordinate_col]/self.block_size
        self.df.loc[: ,'block'] = self.df.loc[: ,'block'].round(decimals=0).astype('int32')
        if self.force_n_block != None:
            self.n_block = self.force_n_block
        else:
            self.n_block =  int(self.df.loc[:, 'block'].max())

    def categorical_variable_block_stat(self, _func: str,
                                             _feature_col: str,
                                             _cond: str,
                                             _norm_factor: list):
        """
        Takes a dataframe , divide it to blocks with a given size and then calculate a
        statistic based on a categorical feature column.
        Args:
            _func (str): the function that is used to calculate the statistic
            _feature_col (pandas Series): a column that contains the feature
            _cond (str): the column name of the column that contains the
                category we want to check.

            self.df ,self.cond, self.feature_col
            docs are at __init__ docs.
        Returns:
            (pd.DataFrame) The statistic per block.
        """
        ### Inputs Dtype chacking ###    
        # assert check_argument_types()

        # Calculating the statistic per block
        # cond == True
        bootstrap_stat_cond = self.df.loc[
            (self.df[_cond] == True), ['block', _feature_col]].groupby(
            'block').agg(_func).copy()
        bootstrap_stat_cond.rename(
            columns={_feature_col:f'{_cond}_{_feature_col}'},
             inplace=True)
        # cond == false
        bootstrap_stat_non_cond = self.df.loc[
            (self.df[_cond] == False), ['block', _feature_col]].groupby(
            'block').agg(_func).copy()
        bootstrap_stat_non_cond.rename(
            columns={_feature_col:f'non_{_cond}_{_feature_col}'},
                                     inplace=True)
        
        block_stat_agg_df = pd.concat(
            [bootstrap_stat_cond ,bootstrap_stat_non_cond], axis=1)
        
        # correcting for the full block range
        block_stat_df = pd.DataFrame(columns=block_stat_agg_df.columns,
                index=[i for i in range(0,(self.n_block +1))])
        block_stat_df.rename_axis('block', inplace=True)
        
        block_stat_df.iloc[
            block_stat_agg_df.index, [0,1]] = block_stat_agg_df.iloc[:,[0,1]]       
        block_stat_df.fillna(value=0, inplace=True)
        
        return self.block_stat_normalizer(data=block_stat_df,
                                    _norm_factor=_norm_factor)


    def block_stat_normalizer(self, data: pd.DataFrame
                                  ,_norm_factor: list):
        """
        Taking a dataframe with statistic per block per pos/neg
        condition and normalize it by factors
        Args:
            data (pd.DataFrame): data to normalize with:
                cond == True as 1st column
                cond == False as 2nd column
            _norm_factor (list): 1st element is the normalization
                factor for the 1st column and 2nd element for the 
                2nd one.
        Returns:
            data (pd.DataFrame): Normalized data 
        """
        data.iloc[:,0] = data.iloc[:,0] / _norm_factor[0]
        data.iloc[:,1] = data.iloc[:,1] / _norm_factor[1]
        return data

    def block_bootstrap_stat_sampler(self, _func: str, _feature_col: str,
                                                _cond: str, _cond_type: str, 
                                                 _fake_block_sample: bool =False,
                                                 _fix_blocks_list: Optional[list]=None,
                                                 calc_ci: bool =False,
                                                 upper_ci: float = 0.95,
                                                 lower_ci: float = 0.05, 
                                                 norm_factor: Optional[list]=None):
        """
        Perform bootstrap from previosly calculated blocks data.

        Args:
            _func, _feature_col, _cond: docs are at categorical_variable_block_stat docs.
            _cond_type (str) : the variable type that we use as a 
                                condition (ether categorial or numerical)
            n (int): number of resamples.
            _fake_block_sample (bool): A flag to force a fake sample (block number 1 only)
                                This is to validate the method.
            _fix_blocks_list (bool): A flag to force a list of sampled blocks.
        Returns:
            self.sampled_stat_df (pd.DataFrame): A dataframe with all the sampled blocks
                    and thier values.
            self.sampled_stat_p_val (float): A p-value that is calculated as follow:
                1.for each sample we check and record weather:
                    statistic(_feature_col == True) > statistic(_feature_col == False).
                    True and False values are represented as 1 and 0 respectively.
                2. Calculate the ratio beween 1/0 counts, this is our p-value.
        """
        ### Inputs Dtype chacking ###    
        # assert check_argument_types()
        
        if _cond_type == 'categorial':
            # getting data frame with:
            # statistic(_feature_col == True) > statistic(_feature_col == False)
            #  per block. 
            block_stat_df = self.categorical_variable_block_stat(_func = _func,
                     _feature_col = _feature_col,
                    _cond = _cond, _norm_factor = norm_factor)
            # print(block_stat_df)
            # creating an empty dataframe to record the samples.
            sampled_df = pd.DataFrame(columns=[f'{_cond}_{_feature_col}',
                    f'non_{_cond}_{_feature_col}'])
            # creating lists to record the p-val and stat_ratio per
            # sampaling iteration
            p_val_record = []
            stat_ratio_ls = []
            self.sampled_blocks_list = []
            for i in range(0,self.n_samples):
                # sampeling the blocks numbers
                sampled_blocks = np.random.choice(
                    a = block_stat_df.index ,
                    size = self.n_block, replace = True)

                self.sampled_blocks_list.append(sampled_blocks)

                if _fix_blocks_list != None:
                    # print(_fix_blocks_list[i])
                    block_stat_sampled_df = block_stat_df.iloc[
                    _fix_blocks_list[i],:].copy()
                else: block_stat_sampled_df = block_stat_df.iloc[
                    sampled_blocks,:].copy()
                
                # calculating a statistic (for the p-val) to compare 
                # _cond == True vs _cond == False
                pos_cond_sum = (block_stat_sampled_df.loc[:,
                    f'{_cond}_{_feature_col}'].sum())

                neg_cond_sum = (block_stat_sampled_df.loc[:,
                    f'non_{_cond}_{_feature_col}'].sum())

                # comparing and recording _cond == True vs _cond == False
                if (pos_cond_sum > neg_cond_sum): n_stat = 1
                else: n_stat = 0
                p_val_record.append(n_stat)
                
                # Calculating the stat_ratio of 
                # sum(_cond == True)/(sum(_cond == False) + sum(_cond == True))
                # ## (Scalar) ##
                if (pos_cond_sum+neg_cond_sum) == 0: sumed_stat_ratio = 0
                else:
                    sumed_stat_ratio = (pos_cond_sum /
                                    (pos_cond_sum + neg_cond_sum))
                
                stat_ratio_ls.append(sumed_stat_ratio)
                sampled_df = pd.concat([sampled_df, block_stat_sampled_df])
            
            sampled_df.reset_index(inplace=True)
            sampled_df.rename(columns={'index': 'block'}, inplace=True)
            
            p_val_record = pd.Series(p_val_record)
            # calculating p-value as the ratio between 1 and 0 in p_val_record
            if len(p_val_record[p_val_record == 0]) == 0: p_val = 1
            else: p_val = len(p_val_record[p_val_record == 1]) / len(p_val_record)
            assert (len(p_val_record) == len(p_val_record[p_val_record == 1]) + len(p_val_record[p_val_record == 0]))
            
            self.sampled_stat_df = sampled_df
            self.sampled_stat_p_val = p_val
            self.stat_ratio = pd.Series(stat_ratio_ls)
            
            # Calculating ci if requiered
            if calc_ci != False:
                self.ci_u = self.stat_ratio.quantile(q=upper_ci)
                self.ci_d = self.stat_ratio.quantile(q=lower_ci)

    def sample_block_num_generator(self):
        """
        Generates a dataframe with lampled block numbers as index
        """
        blocks = [i for i in range(0,self.n_block+1)]
        sample_block_num = np.random.choice(a = blocks , 
            size = (self.n_block* self.n_samples), replace = True)
        samp_block_num_df = pd.DataFrame(index = sample_block_num)
        return samp_block_num_df

