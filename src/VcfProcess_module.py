# importing libraries

import io
import os
import pandas as pd
import re
import gzip
import numpy as np


class VcfProcess:
    def __init__(self, path: str, df_columns : list, df_dtypes: dict,
                 file_format : str, df_columns_drop = None,
                create_dataframe = None, data = None, context_window = None):
        if create_dataframe == False:
            self.data = pd.DataFrame(data)
            self.context_window = context_window
        else:
            # defining data path
            self.path = path    
            self.file_format = file_format
            # defining dataframe columns
            self.df_columns = df_columns
            self.df_dtypes = df_dtypes
            # defining dataframe columns to drop if needed
            self.df_columns_drop = df_columns_drop
            # calling create_dataframe module to create instance 
            self.data = self.create_dataframe()
        
        
        
    def create_dataframe(self) -> None:
        if self.file_format == 'vcf.gz':
            # opening the .bed.gz file as a pandas dataframe
            with gzip.open(self.path, 'rb') as file:
                bed_df = file.read()
            # transforming file to pandas datframe
            df = pd.read_csv(io.BytesIO(bed_df), sep=" ", header = None) 
            df.columns = self.df_columns # defining column names
             # droping empty columns (generated in the uinix script, 
             # during intersection)

            if self.df_columns_drop != None :
                df = df.drop(columns = self.df_columns_drop)    
        else:
            df = pd.read_csv(self.path, sep="\t", dtype = self.df_dtypes)
        
        return df


    """
    freq_extraction will take data and column as an argumants 
    and will return only the numerical value out of X:value format
    """    
    def calculate_frequancy_and_add_to_DataFrame(self, column: str):
        self.column = column
        freq_vector = self.data.loc[:, self.column].str.split(pat=':', 
                                                              expand = True)[1].astype('float32')
        self.data[f'{column}_value'] = freq_vector

    """
    indel_categorizer will take data, REF_freq_column, ALT_freq_column, 
    REF, ALT and will return an 'indel' column that identify
    the insertion and deletions
    """
    def add_indel_type(self, REF_freq_column: str, 
                       REF_freq_cutoff: float, 
                       ALT_freq_column: str, 
                       ALT_freq_cutoff: float, 
                       REF: str, ALT: str):
        REF_ancestral = self.data[REF_freq_column].astype('float') > REF_freq_cutoff
        ALT1_ancestral = self.data[ALT_freq_column].astype('float') > ALT_freq_cutoff
        REF_len = self.data[REF].str.len()
        ALT_len = self.data.loc[:, ALT].str.split(pat=',', expand = True)[0].str.len()
        indel = (ALT1_ancestral == True) & (ALT_len > 1)
        indel = indel.replace(to_replace = False, value = 'insertion')
        indel = indel.replace(to_replace = True, value = 'deletion')
        self.data['indel_type'] = indel
        
    """
    add_fasta_context_to_dataframe will get a location of a .fa file that
    has generated in a combination of bedtools-getfasta
    (see documentation on/home/labs/alevy/guyta/varitom_tabix/varitom_indels/
    analyzed/indel_data_to_python/chr06_analysis/chr06_getfasta50bp_bedtools/)
    and the output of transform_to_context_extraction_format, and will add 
    two columns to the dataframe: fasta_context_position, fasta_context_seq.
    """
    def add_fasta_context_to_dataframe(self, context_file: str, context_window: int):
        self.context_window = context_window
        with open(context_file, 'r') as fasta_context: # loading the file
            # create a dataframe out of the file
            context_df = pd.read_csv(fasta_context, header=None)
            # take only the rows that coresponse to 
            # context location in the fasta reference
            fasta_context_position = np.array(context_df.iloc[::2, 0])
            # take only the rows that coresponse to context sequance in the fasta reference
            fasta_context_seq = np.array(context_df.iloc[1::2, 0]) 
            # transforming into 1D vector
            fasta_context_position = pd.Series(fasta_context_position) 
            # transforming into 1D vector
            fasta_context_seq = pd.Series(fasta_context_seq) 
            # calculating REF position
            fasta_ref_position = int(len(fasta_context_seq[1])/2)
            # extracting the referance nucleotide from the fasta context
            fasta_context_checker = fasta_context_seq.str.get(fasta_ref_position)
            if self.file_format == 'vcf.gz':
                self.data['fasta_context_checker'] = fasta_context_checker
                self.data['fasta_context_position'] = fasta_context_position 
                self.data['fasta_context_seq'] = fasta_context_seq
            else:
                return (pd.DataFrame({'fasta_context_checker': fasta_context_checker,
                                       'fasta_context_position': fasta_context_position,
                                       'fasta_context_seq' : fasta_context_seq}))

        
    def extract_contaxt(self, sequance: str, startpos: int,
                        endpos: int, proximity_range: int): 
        contaxt = (sequance[(startpos - proximity_range) 
                            : (endpos + proximity_range)])
        assert 0 <= (startpos-proximity_range) , 'Out of lower range'
        assert len(seq) >= (endpos+proximity_range) , 'Out of upper range'
        return contaxt

    """
    transform_to_context_extraction_format will transform the 
    VcfProcess object and transform it into a format wich contain
    only the necessary columns for bedtools getfasta command to
    get the context of the indel based on its position
    """
    def transform_to_context_extraction_format(self, chr_column: str,
                                               start_pos: str, 
                                               end_pos: str, 
                                               proximity_range: int):   
        chr_num =self.data[chr_column].str.split(pat='ch',expand = True)[1].astype('int')
        self.data['chr_number'] = chr_num
        self.data['chr_ref_fasta_format'] = 'chr' + self.data['chr_number'].astype(str)
        self.data['context_start_pos'] = self.data['start_pos'] - proximity_range
        self.data['context_end_pos'] = self.data['end_pos'] + proximity_range
        
        return self.data.loc[:,['chr_ref_fasta_format', 
                                'context_start_pos', 
                                'context_end_pos']]

        
    """
    general_indels_stats will take only self and will calculate indels 
    counts and precentage, and will serve as a cheking point after 
    addition of the fasta context
    """
    def general_indels_stats(self, fasta_available = True):
        # calculating counts and presentage of deletions out of the data
        varitiom_vhr06_deletions = self.data['indel_type'] == 'deletion'
        varitiom_vhr06_deletions_counts = len(self.data[varitiom_vhr06_deletions].index)
        varitiom_vhr06_deletions_percentage = round(varitiom_vhr06_deletions_counts 
                                                    / len(self.data.index), 5)
        print('# number of deletions in the data: ', varitiom_vhr06_deletions_counts, 
              '(', varitiom_vhr06_deletions_percentage ,'% of the data)')
        # calculating counts and presentage of insertions out of the data
        varitiom_vhr06_insertion = self.data['indel_type'] == 'insertion'
        varitiom_vhr06_insertion_counts = len(self.data[varitiom_vhr06_insertion].index)
        varitiom_vhr06_insertion_percentage = round(varitiom_vhr06_insertion_counts 
                                                    / len(self.data.index), 5)
        print('# number of insertion in the data: ',varitiom_vhr06_insertion_counts, 
              '(', varitiom_vhr06_insertion_percentage,'% of the data)')
        if fasta_available == True:
        # comparing the original REF column to the one thet had calculated from 
        # the fasta file and state if they are the same or not
            print('# NO mismatches between the original REF column ',
                  'and the fasta reference (fasta_context_checker): ', 
                  any(self.data['REF'] == self.data['fasta_context_checker']))
        else: 
            print('Fasta referance not given')

                    

    """
    accession_context_generator function will generate the sequence of the accession
    as it should apear in the accession. for insertions: it basicaly take the context
    and insert the insertion in position self.context_window. for deletions: it finds the delerion 
    in the context sequence and omit it -> so only the context without the deletion left.
    """

    def accession_context_generator(self,indel_type : str, alt : str, ref : str):
#         self.data['indel_pos'] = 0
        # defining an inner function to parse deletions
        def deletion_context_generator(_context, _ref):
            full_inserion_context = (_context[ : (self.context_window - 1)] 
                                     + _context[(self.context_window -1 + len(_ref)) : ])
#             full_deletion_context = _context.replace('-', '')
            return full_inserion_context
        
        # defining an inner function to parse insertions
        def insertion_context_generator(_context, _alt):
#             marker = '*'
            full_inserion_context = (_context[ : (self.context_window + 1)] 
                                     + _alt 
                                     + _context[(self.context_window +1) : ]) 
            return full_inserion_context

        # defining two conditions in order to get pandas seriess with only
        # insertions and deletions
        _ins = self.data[indel_type] == 'INS'
        _del = self.data[indel_type] == 'DEL'
        
        # applying the pre-defined inner functions on the right rows
        accession_context_insetion = self.data[_ins].apply(lambda row:
                       insertion_context_generator(row['fasta_context_seq'],
                                                              row[alt]),
                                                           axis=1)
        accession_context_deletion = self.data[_del].apply(lambda row:
                        deletion_context_generator(row['fasta_context_seq'],
                                                             row[ref]),
                                                           axis=1)
        
        # droping NaN that are a result of the inner functions wen the conditions
        # are not met (i.e. when insertion_context_generator itereting over deletion
        # row and vice viersa)
        accession_context_deletion = accession_context_deletion.dropna()
        accession_context_insetion = accession_context_insetion.dropna()
        # creating one series from the deletion/ insertion
        accession_context = accession_context_deletion.append(accession_context_insetion)
        # sorting by index
        accession_context = accession_context.sort_index() 
        # concat the new series to the data
        self.data.loc[:, ['accession_context']] = accession_context
    
    
    
    """
    add_indel_pos_to_data will get the indel position inside the 
    reference context
    """
    def add_indel_pos_to_data(self, indel_type : str, ref : str):
        # memory allocation
        self.data['indel_pos'] = 0
        # insertions position == the window that have been used when
        # using - bedtools getfasta 
        self.data.loc[(self.data.indel_type == 'INS'), 'indel_pos'] = self.context_window + 1
        ### NEW 06.10.2021 ###
        self.data.loc[(self.data.indel_type == 'DEL'), 'indel_pos'] = (self.context_window - 1)
