# 03.03.2022
"""
This script contains all the functions that toghether perform
a block bootstrap on a given dataframe
"""

# import libreries
import sys
from unicodedata import name
import numpy as np
import pandas as pd

class BlockBootstrap:
    def __init__(self, data: pd.DataFrame, _block_size: int, 
            _coordinate_col: str,
            _n_samples: int,_cutoff = None):
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
            _cutoff (float): for numerical categories (_cond is numerical), 
                this is a cutoff to seperate the _cond column on.
        """
        
        # Checking that inputs are of the right type
        assert (type(data) == pd.DataFrame), f'Input of wrong type: _df = {data} (type: {type(data)})'
        assert (type(_block_size) == int), f'Input of wrong type: _block_size = {_block_size} (type: {type(_block_size)})'
        assert (type(_coordinate_col) == str), f'Input of wrong type: _coordinate_col = {_coordinate_col} (type: {type(_coordinate_col)})'
        assert (type(_n_samples) == int), f'Input of wrong type: _n_samples = {_n_samples} (type: {type(_n_samples)})'
        if _cutoff != None:
            assert (type(_cutoff) == float), f'Input of wrong type: _cutoff = {_cutoff} (type: {type(_cutoff)})'
        
        # assining inputs as atributes
        self.df, self.block_size, self.coordinate_col = data, _block_size, _coordinate_col
        self.n_samples = _n_samples
        if _cutoff != None:
            self.cutoff = _cutoff
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
        self.n_block =  self.df.loc[:, 'block'].nunique()


    def create_block_bootstrap_df(self):
        """
        An implementation of the block bootstrap technique: 
            Deviding the genome into block with fixed size -> calculate the the
            fiture that we want to sample -> sample the blocks with replacment
            N times.

        Args:
            self.df, self.n_block, self.n_samples
            docs are at __init__ docs.
        Returns:
            (pd.DataFrame) A dataframe (self.bootstrap_df) with all
            the bootstraped data.
        """

        sampled_blocks = np.random.choice(a = self.df.loc[:, 'block'] ,
                size = (self.n_block * self.n_samples), replace = True)
        # getting the index for all sampled_blocks enteries.
        indexes = []
        for b in sampled_blocks:
            indexes = indexes + self.df.loc[(self.df['block'] == b), :].index.tolist()
        
        self.bootstrap_df = self.df.loc[indexes,:]
        # print(self.bootstrap_df.info())
        return self.bootstrap_df


    def categorical_var_block_bootstrap_stat(self, _func: str, _feature_col: str,
                                                _cond: str):
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
        # validating input's types
        assert (type(_feature_col) == str), f'Input of wrong type: _feature_col = {_feature_col} (type: {type(_feature_col)})'
        assert (type(_func) == str), f'Input of wrong type: _func = {_func} (type: {type(_func)})'
        assert (type(_cond) == str), f'Input of wrong type: _cond = {_cond} (type: {type(_cond)})'

        # Calculating the statistic per block
        bootstrap_stat_pos = self.df.loc[self.df[_cond] == True, ['block', _feature_col]].groupby(
            'block').agg(_func).copy()
        bootstrap_stat_pos.rename(columns={_feature_col:f'{_cond}_{_feature_col}'}, inplace=True)
        bootstrap_stat_neg = self.df.loc[self.df[_cond] == False, ['block', _feature_col]].groupby(
            'block').agg(_func).copy()
        bootstrap_stat_neg.rename(columns={_feature_col:f'non_{_cond}_{_feature_col}'}, inplace=True)

        # print(bootstrap_stat_pos)
        # print(type(bootstrap_stat_neg))
        return pd.concat([bootstrap_stat_pos ,bootstrap_stat_neg], axis=1)
        

    def block_bootstrap_stat_sampler(self, _func: str, _feature_col: str,
                                                _cond: str, _cond_type: str, 
                                                 _fake_block_sample = False,
                                                 _fix_blocks_list = False):
        """
        Perform bootstrap from previosly calculated blocks data.

        Args:
            _func, _feature_col, _cond: docs are at categorical_var_block_bootstrap_stat docs.
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
        assert (type(_cond_type) == str), f'Input of wrong type: _cond_type = {_cond_type} (type: {type(_cond_type)})'
        
        if _cond_type == 'categorial':
            # getting data frame with:
            # statistic(_feature_col == True) > statistic(_feature_col == False)
            #  per block. 
            block_stat_df = self.categorical_var_block_bootstrap_stat(_func = _func,
                     _feature_col = _feature_col,
                    _cond = _cond)

            # creating an empty dataframe to record the samples.
            sampled_df = pd.DataFrame(columns=[f'{_cond}_{_feature_col}',
                    f'non_{_cond}_{_feature_col}'])
            p_val_record = []
            self.sampled_blocks_list = []
            for i in range(0,self.n_samples):
                # sampeling the blocks numbers
                sampled_blocks = np.random.choice(
                    a = self.df.loc[:, 'block'] ,
                    size = self.n_block, replace = True)
                self.sampled_blocks_list.append(sampled_blocks)
                # Checking if spliting data before sampaling is same as
                # splitinig after sampling
                # seting artifitial sample:
                if _fake_block_sample == True:
                    sampled_blocks = [1 for i in range(0,self.n_samples)]
                    # print(f'##{sampled_blocks}')
                if _fix_blocks_list != False:
                    # print(_fix_blocks_list[i])
                    block_stat_sampled_df = block_stat_df.iloc[
                    _fix_blocks_list[i],:].copy()

                else: block_stat_sampled_df = block_stat_df.iloc[
                    sampled_blocks,:].copy()
                
                # calculating a statistic to compare _cond == True vs _cond == False
                if _func == 'sum':
                    pos_cond_sum = block_stat_sampled_df.loc[:,
                    f'{_cond}_{_feature_col}'].sum()
                    neg_cond_sum = block_stat_sampled_df.loc[:,
                    f'non_{_cond}_{_feature_col}'].sum()
                if _func == 'mean':
                    pos_cond_sum = block_stat_sampled_df.loc[:,
                        f'{_cond}_{_feature_col}'].mean()
                    neg_cond_sum = block_stat_sampled_df.loc[:, 
                        f'non_{_cond}_{_feature_col}'].mean()
                
                # comparing and recording _cond == True vs _cond == False
                if (pos_cond_sum > neg_cond_sum): n_stat = 1
                else: n_stat = 0
                p_val_record.append(n_stat)
                sampled_df = pd.concat([sampled_df, block_stat_sampled_df])
            
            sampled_df.reset_index(inplace=True)
            sampled_df.rename(columns={'index': 'block'}, inplace=True)
            
            p_val_record = pd.Series(p_val_record)
            # calculating p-value as the ratio between 1 and 0 in p_val_record
            if len(p_val_record[p_val_record == 0]) == 0: p_val = 1
            else: p_val = (len(p_val_record[p_val_record == 1]) / 
                        (len(p_val_record[p_val_record == 0] + 
                        len(p_val_record[p_val_record == 1]))))
            self.sampled_stat_df = pd.DataFrame(sampled_df) 
            self.sampled_stat_p_val = p_val           









    # def numeric_var_block_bootstrap_stat(_df : pd.DataFrame, cond: str,
    #                 _columns : list, _feature_col: str ,block_size: int, coordinate_col : pd.Series,
    #                 _func:str, cutoff : float):
    #     """
    #     Takes a dataframe , divide it to blocks with a given size and then calculate a
    #     statistic based on a numerical feature column.
    #     Args:
    #         _df (pd.DataFrame): The data.
    #         cond (str): the column name of the column that contains the
    #             category.
    #         _columns (list): a list of relevant columns that one want to consider.
    #         _feature_col (str): the column name on which the statistic will be calculated.
    #         _block_size (int): size of the block
    #         coordinate_col (str): the column name of the coordinate column.
    #         _func (str): the function that is used to calculate the statistic
    #     Returns:
    #         (pd.DataFrame) The statistic per block.
    #     """
    #     # Assining blocks to data
    #     _df = _df.loc[:, _columns].copy()
    #     _df.loc[: ,'block'] = _df.apply(lambda x: 
    #         coordinates_to_blocks(x[coordinate_col], block_size), axis = 1)

    #     bootstrap_stat_pos = _df.loc[_df[cond] > cutoff, ['block', _feature_col]].groupby(
    #         'block').agg(_func).copy()
    #     bootstrap_stat_neg = _df.loc[_df[cond] < cutoff, ['block', _feature_col]].groupby(
    #         'block').agg(_func).copy()
    #     print(bootstrap_stat_pos)
    #     print(bootstrap_stat_neg)

    #     # n_blocks = _df.loc[:, 'block'].nunique()
    #     # sampled_blocks = np.random.choice(a = _df.loc[:, 'block'] ,
    #     #         size = (n_blocks * n_samples), replace = True)
