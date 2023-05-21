# 13.01.2022
"""
This file is a collection of functions that I often use
when analysing data.
"""
import gzip
import io
import os
import glob

import pandas as pd
import numpy as np
import scipy.stats as stats
# from seqfold import dg
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
# from typeguard import check_argument_types



def indel_bin(_df, _col, _rnages: list):
    """
    binning indel_length
    Args:
        _df (pandas DataFrame) : the data
        _col (pandas Series): coordinate column
        _rnages (list): a list of ranges to bin by

    """
    # set the ranges
    shorter_then_10 = range(_rnages[0],_rnages[1])
    _10to15 = range(_rnages[1],_rnages[2])
    _15to30 = range(_rnages[2],_rnages[3])
    longer_then_30 = range(_rnages[3],_rnages[4])
    # set conditions
    conditions =[
        _df[_col].isin(shorter_then_10),
        _df[_col].isin(_10to15),
        _df[_col].isin(_15to30),
        _df[_col].isin(longer_then_30)
    ] 
    # set choices
    choices = [
        _df[_col],
        '10<=Indel<15',
        '15<=Indel<30',
        '30<=Indel'
    ]
    # set new index
    new_index = ['2','3','4','5','6','7','8','9',
                 '10<=Indel<15',
        '15<=Indel<30',
        '30<=Indel']
    # assign the binned data to a new column
    binned_col = np.select(conditions, choices)
#     data['indel_len_bin'] = np.select(conditions, choices)
    return binned_col


def coordinates_converter(start: int, interval_size: int):
    """
    This function taking the indel's coordinates and round it by the
    interval size in order to create binns along the chromosome that
    represents the position of groups of indels along the chromosome.
    Args:
        start (int): the coordinate column.
        interval_size (int): the interval/block size.
    """
    return (int(round((start/interval_size), 0)))

def bin_avg_value(_df, coordinate_col, feature_col,
                func = 'mean',add_to_df = False, new_col_name = None,
                data_to_append = None, normalize = False, norm_by = False,
                keep_coordinate_ranges = False, manual_range = True):
    """
    This function takes a categorical column (coordinate_col) and a
    numerical column (feature_col) and generates a dataframe with the
    following structure: coordinate_col, avg_values_per_chategoty
    Then if 'add_to_df' == True, it adds the avg values to the df

    Args:
        _df (pandas DataFrame) : the data
        coordinate_col (pandas Series): coordinates of the data points
        feature_col (pandas Series): column with the feature that we want to bin
        func (str): name of the function to aggregate with
        add_to_df (bool): a flage that activate adding the new column to the _df
        new_col_name (str): if add_to_df == True -> this will be the column name
        data_to_append (pandas DataFrame): if add_to_df == True -> this will be the data
        normalize (bool): normalization flag
        norm_by (str) : if normalize == True -> this is the normalization method
        keep_coordinate_ranges (bool): flag to keep the original coordinate range
        manual_rnage (bool): inputing manual range

    """
    if keep_coordinate_ranges == True:
        r = range(_df.loc[:,coordinate_col].min(),(_df.loc[:,coordinate_col].max() + 1),1)
        d = pd.DataFrame(r, columns = [coordinate_col])

    if manual_range != True:
        if len (manual_range) == 2 : r = range(manual_range[0],(manual_range[1] + 1),1)
        if len (manual_range) == 4: r = ([i for i in range(manual_range[0],(manual_range[1] + 1),1)]
             + [i for i in range(manual_range[2],(manual_range[3] + 1),1)])
        d = pd.DataFrame(r, columns = [coordinate_col])
    
    avg_values_per_chategoty_df = _df.loc[:,[feature_col,coordinate_col]]
    coordinate_col_counts = (avg_values_per_chategoty_df[
        coordinate_col].value_counts())
    
    avg_values_per_chategoty_df = avg_values_per_chategoty_df.groupby(
        coordinate_col).agg(func)
    avg_values_per_chategoty_df['value_counts'] = coordinate_col_counts.sort_index()
    avg_values_per_chategoty_df.reset_index(inplace = True)
    if keep_coordinate_ranges == True:
        d = d.merge(avg_values_per_chategoty_df, how = 'left',left_on=coordinate_col, right_on=coordinate_col)
        
        d.fillna(value = 0, inplace = True)
        avg_values_per_chategoty_df = d

    if normalize == True:
        if norm_by == 'counts':
            ### TODO THIS SOULD NORM BY COUNTS PER COORDINATE (INDEL_LENGTH FOR EXAMPLE) ###
            avg_values_per_chategoty_df[feature_col] = (
                avg_values_per_chategoty_df[feature_col] / avg_values_per_chategoty_df[feature_col].sum()
            )    
        else: avg_values_per_chategoty_df[feature_col] = (
            avg_values_per_chategoty_df[feature_col] / avg_values_per_chategoty_df['value_counts'])
    
    if add_to_df == True:
        data_to_append[new_col_name] =  data_to_append.apply(lambda x: (
            avg_values_per_chategoty_df.loc[x[coordinate_col], feature_col]), axis = 1)
    
    
    return avg_values_per_chategoty_df

def GC_contant_calc(_seq, _seq_len):
    """
    This function calculates the GC content of a given sequence
    Args:
        _seq (str): the sequence
        _seq_len (int): the sequence length
    """
    return ((_seq.count('G') + _seq.count('C')) / _seq_len)



# def add_dg(ref_seq, alt_seq, window):
#     """
#     add_dg will calculate the dg (minimum free energy) using the dg function
#     from package seqfold. 
#     Args:
#         ref_seq -> the reference genome seq
#         alt_seq -> the seq with the indel
#         window -> am int that defin the window size on which the calculation
#             will take place (indel pos +/- window).
#             example: if window == 50 -> seq = -50->indel_pos -> +50
#                 len(seq) == 100
#     """
#     _ref_seq = ref_seq[int(round(len(ref_seq)/2)-window) : int(round(len(ref_seq)/2)+window)]
#     _alt_seq = alt_seq[int(round(len(alt_seq)/2)-window) : int(round(len(alt_seq)/2)+window)]
#     _ref_seq_dg = dg(_ref_seq, temp = 37.0)
#     _alt_seq_dg = dg(_alt_seq, temp = 37.0)
#     delta_dg = _ref_seq_dg - _alt_seq_dg
#     return [_ref_seq_dg, _alt_seq_dg, delta_dg]

def add_Tm_NN(_ref_seq, _window):
    """
    add_dg will calculate the melting temperature using the mt.Tm_NN
    from package MeltingTemp (from Bio). 
    Args:
        _ref_seq (str): the reference genome seq
        _window (int): defin the window size on which the calculation
            will take place (indel pos +/- window).
            example: if window == 50 -> seq = -50->indel_pos -> +50
                len(seq) == 100
    """
    seq = _ref_seq[int(round(len(_ref_seq)/2)-_window) : int(round(len(_ref_seq)/2)+_window)]
       
    return mt.Tm_NN(Seq(seq))
   
def annotate(data):
    """
    annotate will add N(sample size) to plots.
    Args:
        data: the data. 
    Returns:
        Add the sample size to the plot
    """
    n = len(data)
    ax = plt.gca()
    ax.text(0, 0.9, f"N = {n}", transform=ax.transAxes)

def add_certain_mechanism_col(data: pd.DataFrame, prob_col_list: list,
                              relative_prob_col_list: list, cols_list:list, prob_cutoff: float,
                              relative_prob_cutoff: float) -> None:
    for col, prob_col, rel_prob_col in zip(cols_list, prob_col_list, relative_prob_col_list):
        data.loc[:,f'{col}_certain'] = False
        
        data.loc[((data[prob_col] >= prob_cutoff) & 
                  (data[rel_prob_col] >= relative_prob_cutoff)),
                     f'{col}_certain'] = True


def mem_reporter():
    """
    Reports the RSS and VMS of a process
    No args.
    No return, just reporting VMS and RSS at a current time point
    """
    import os, psutil, datetime
    # from time import strftime, localtime
    # local_d = str(datetime.datetime.today()).split()[0]
    # local_t = strftime("%H:%M:%S", localtime())
    # print(f'# Time: {local_d}, {local_t}')
    print(f'# VMS(Mb):{psutil.Process(os.getpid()).memory_info().vms / 1024 ** 2}, # RSS(Mb):{psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2}')
    # print(f'# RSS(Mb):{psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2}')



def get_file_list(path_to_files: str, 
                file_format: str) -> list:
    """
    This function deemend importing os and glob.
    It changes the working directory and look for files
    with a given format
    Args:
        path_to_files (str): the folder that contains the files.
        file_format (str): the file format (csv for example)
    Returns:
        A list with all the files coresponding to the given format
        in the given location 
    """
    os.chdir(path_to_files)
    return glob.glob('*.{}'.format(file_format))

######## Statistic calculations ##########

def var_calc(vec: pd.Series) -> float:
    """
    This function calculates the variance of a vector
    using the following formula:
    sum((X - mean(vec)**2))/len(vec)
    Args:
        vec (pd.Series): vector of numeric values
    Returns:
        var (float): the variance of the vector
    """
    vec_mean = vec.mean()
    vec_as_df = pd.DataFrame([vec]).T
    vec_as_df['vec_minus_mean'] = (vec_as_df.iloc[:,0] - vec_mean)**2
    var = vec_as_df['vec_minus_mean'].sum() / len(vec)
    return var
    

def cov_calc(vec_a: pd.Series, vec_b: pd.Series) -> float:
    """
    This function calculate the covariance between two vectors
    using the following formula (Pearson):
    Corr = (E[((vec_a - mean(vec_a) * (vec_b - mean(vec_b))] / (var(vec_a) * var(vec_b))**0.5)**2
    Args:
        vec_a/b (pd.Series): vectors of numeric variables
    Returns:
        R_squere (float): the covariance between the two vectors 
    """
    cov_df = pd.DataFrame()
    cov_df['vec_a'], cov_df['vec_b'] = vec_a, vec_b
    cov_df['vec_a_minus_mean_a'] = cov_df.loc[:,'vec_a'] - cov_df.loc[:,'vec_a'].mean()
    cov_df['vec_a_minus_mean_b'] = cov_df.loc[:,'vec_b'] - cov_df.loc[:,'vec_b'].mean()
    # calculating the expectation
    E = sum(cov_df['vec_a_minus_mean_a'] * cov_df['vec_a_minus_mean_b']) / len(
        cov_df['vec_a_minus_mean_b'])

    corr = (E/((var_calc(vec_a) * var_calc(vec_b))**0.5))
    R_squere = corr**2
    return corr , R_squere

