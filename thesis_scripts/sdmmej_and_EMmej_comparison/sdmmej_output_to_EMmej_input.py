"""
This sctipt converts sdmmej output into a VCF
format that works with EMmel
"""
import os

# import regex as re
import pandas as pd
import numpy as np

pd.options.display.max_colwidth = 3500
# pd.set_option('display.max_rows', 1000)
pd.set_option("display.max_columns", None)
output_path = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input'
date = '20221008'
date = '20221030'


# ------------- R0_WT ------------- 
dir = os.path.join(f"{output_path}/R0_WT/{date}")
if not os.path.exists(dir):
    os.mkdir(dir)

ref = 'GATCCTAGGATCCTAGGATCCTAGGATCCTAGGATCCTAGGATCCTAGGATCCTAGGATCCTAGGATCCTAGGATCCTAGGAGGGAAAAAATTCGTACTTTGGAGTACGAAATGCGTCGTTTAGAGCAGCAGCCGAATTCGGTACATTACCCTGTTATCCCTAGCGGCCGCATAGGCCACTAGTGGATCTGGATCCTCTAGAGTCGACCTCGAACGTTAACGTTAACGTAACGTTAACTCGAGGCCCAGGACCCCAGGACCCCAGGACCCCAGGACCCCAGGACCCCAGGACCCCAGGACCCCAGGACCCCAGGACCCCAGGACC'
cutsite = 161
# setting the amplicon context in the plasmid
pre_insert_1000  = 'CTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAATGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGCCGGTTCAGACAGGATAAAGAGGAACGCAGAATGTTAGACAACACCCGCTTACGCATAGCTATTCAGAAATCAGGCCGTTTAAGCGATGATTCACGAGAATTGCTGGCCCGCTGCGGCATAAAAATTAATTTACACACTCAGCGCTGATGAATCCCCTAATGATTTTGGTAAAAATCATTAAGTTAAGGTGGACACACATCTTGTCATATGATTAAATGGTTTCGCGAAAAATCAATAATCAGACAACAAGATGTGCGAACTCGATATTTTACACGACTCTCTTTACCAATTCTGCCCCGAATTACACTTAAAACGACTCAACAGCTTAACGTTGGCTTGCCACGCATTACTTGACTGTAAAACTCTCACTCTTACCGAACTTGGCCGTAACCTGCCAACCAAAGCGAGAACAAAACATAACATCAAACGAATCGACCGATTGTTAGGTAATCGTCACCTGCAGGAAGGTTTAAACGCATTTAGGTGACACTATAGAAGTGTGTATCGCTCGAGGGATCCGAATTCAGGAGGTAAAAACCATGAT'
post_insert_1000 = 'CTGATAATAATTAATTAAGACGTCAGAATTCTCGAGGCGGCCGCATGTGCGTCTCCCTATAGTGAGTCGTATTAATTTCGCGGGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGATCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGA'
indel_pos_correction = 1000 - 79 # This number is a function of fow you generates the reference genome(+) and the padding(-)

# Deletions as appear in sdmmej output: Those are all deletion classified by sdmmej
R0_WT_outputs = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/R0_WT_data_from_Nick/raw_data_to_sdmmej_output/1023_20220903_Iw7/bwa_output/bwa_output_Final_with_headers_output'
path_to_sdmmej_del = f'{R0_WT_outputs}/table_outputs/bwa_output_Final_with_headers_deletions_table.csv'
sdmmej_del = pd.read_csv(path_to_sdmmej_del)

sdmmej_del['deltion_seq'] = sdmmej_del.apply(lambda row: ref[row['LEFT_DEL_INDEX']:row['RIGHT_DEL_INDEX']], axis=1)
sdmmej_del_vcf = sdmmej_del.copy()

sdmmej_del_vcf['REF'] = sdmmej_del.apply(lambda row: ref[(row['LEFT_DEL_INDEX']-1):row['RIGHT_DEL_INDEX']], axis=1)
sdmmej_del_vcf['ALT'] = sdmmej_del.apply(lambda row: ref[(row['LEFT_DEL_INDEX']-1):row['LEFT_DEL_INDEX']], axis=1)
sdmmej_del_vcf['POS'] = sdmmej_del_vcf['LEFT_DEL_INDEX']
sdmmej_del_vcf['CHR'] = [f'Iw7_D_{i}' for i in sdmmej_del_vcf.index]
sdmmej_del_vcf['del_seq'] = np.nan
sdmmej_del_vcf.rename(columns={'CONSISTENCY': 'consistency'}, inplace=True)

# Data for comparison
del_cols = ['CHR', 'POS', 'REF', 'ALT', 'RECONSTRUCTED_SEQ', 'REPAIR_TYPE', 'LEFT_DEL_INDEX',
                   'RIGHT_DEL_INDEX', 'consistency', 'CLASS', 'CLASS_final',
                   'READS', 'MICROHOMOLOGY', 'MH_Length']
sdmmej_del_vcf.loc[:, del_cols]


# insertions that was classified by HiFIBR as complex and then reclassified as ins by sdmmej
path_to_comp_ins = f'{R0_WT_outputs}/bwa_output_Final_with_headers_complex_insertion_consistency2.csv'
comp_ins_df = pd.read_csv(path_to_comp_ins)
comp_ins_df_vcf = comp_ins_df.copy()

comp_ins_df_vcf['ALT'] = comp_ins_df_vcf.apply(lambda row:row['RECONSTRUCTED_SEQ'][(row['left_del']-1): (row['right_del']-1)], axis=1)
comp_ins_df_vcf['REF'] = comp_ins_df_vcf.apply(lambda row:row['RECONSTRUCTED_SEQ'][(row['left_del']-1): row['left_del']], axis=1)
comp_ins_df_vcf['POS'] = comp_ins_df_vcf['left_del'] # 0 based
comp_ins_df_vcf = comp_ins_df_vcf.loc[comp_ins_df_vcf['POS'] != 0, :]
comp_ins_df_vcf['CHR'] = [f'Iw7_CI_{i}' for i in comp_ins_df_vcf.index]

# insertions that was classified by HiFIBR as insertion and then reclassified as ins by sdmmej
path_to_ins_ins = f'{R0_WT_outputs}/bwa_output_Final_with_headers_insertion_insertion_consistency2.csv'
ins_ins = pd.read_csv(path_to_ins_ins)

ins_ins_vcf = ins_ins.copy()
ins_ins_vcf['ALT'] = ins_ins.apply(lambda row:row['RECONSTRUCTED_SEQ'][(row['left_del']-1): (row['right_del']-1)], axis=1)
ins_ins_vcf['REF'] = ins_ins.apply(lambda row:row['RECONSTRUCTED_SEQ'][(row['left_del']-1): row['left_del']], axis=1)
ins_ins_vcf['POS'] = ins_ins['left_del'] # 0 based
ins_ins_vcf['CHR'] = [f'Iw7_II_{i}' for i in ins_ins_vcf.index]

del_cols = ['CHR', 'POS', 'REF', 'ALT','ID', 'RECONSTRUCTED_SEQ','del_seq', 'REPAIR_TYPE', 'consistency', 'CLASS','CLASS_final','READS', 'MICROHOMOLOGY', 'MH_Length', 'LEFT_DEL_INDEX', 'RIGHT_DEL_INDEX']
ins_cols = ['CHR', 'POS', 'REF', 'ALT','ID', 'RECONSTRUCTED_SEQ', 'del_seq','consistency','left_del','right_del', 'DRmotif_length', 'RCmotif_length', 'Loop-out', 'Snap-back']
comb_df = pd.concat([sdmmej_del_vcf.loc[:, del_cols], 
          ins_ins_vcf.loc[:, ins_cols],
         comp_ins_df_vcf.loc[:,ins_cols]]).reset_index(drop=True)
# Exporing data
comb_df['POS'] = comb_df['POS'] + indel_pos_correction
comb_df.rename(columns={'CHR':'#CHR'}, inplace=True)

comb_df.to_csv(f'{output_path}/R0_WT/{date}/{date}_R0_WT_sdmmej_output_EMmej_input.tsv',
                index=False, sep='\t')

# Exporing the VCF format for EMmej
comb_df.loc[:, ['#CHR', 'POS', 'REF', 'ALT']].to_csv(f'{output_path}/R0_WT/{date}/{date}_R0_WT_sdmmej_output_EMmej_input.vcf',
                index=False, sep='\t')

# # ----------- Creating a 'reference genome' ----------- 
with open(f'{output_path}/R0_WT/{date}/{date}_R0_WT_sdmmej_output_EMmej_input_ref.fa', 'w') as f:
    for i in comb_df.loc[:198, :].index:
        f.write(f">\t{comb_df.loc[i, '#CHR']}\n")
        f.write(f'{pre_insert_1000}{ref[79:-79]}{post_insert_1000}\n')
    for i in comb_df.loc[199:, :].index:
        f.write(f">\t{comb_df.loc[i, '#CHR']}\n")
        f.write(f"{pre_insert_1000}{comb_df.loc[i, 'del_seq'][79:-79]}{post_insert_1000}\n")

# Writing a REDME file
README = f"""This folder contains the output of: 
/home/labs/alevy/guyta/guy_master_project/scripts/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input.py
The inputs are:
    sdmmej deletions : {path_to_sdmmej_del}
    sdmmej complex insertions: {path_to_comp_ins}
    sdmmej insertions: {path_to_ins_ins}

The outputs are:
    {date}_R0_WT_sdmmej_output_EMmej_input.vcf : The data in a VCF format for EMmej
    {date}_R0_WT_sdmmej_output_EMmej_input.tsv : The data with informarive columns for later use in the comparison
    {date}_R0_WT_sdmmej_output_EMmej_input_ref.fa : A reference genome for the data

"""
with open(f'{output_path}/R0_WT/{date}/README.md', 'w') as f:
    f.write(README)




# ------------- R0_POLQ ------------- 
dir = os.path.join(f"{output_path}/R0_POLQ/{date}")
if not os.path.exists(dir):
    os.mkdir(dir)

# Deletions as appear in sdmmej output: Those are all deletion classified by sdmmej
R0_POLQ_outputs = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/R0_POLQ/raw_reads_to_sdmmej_output/1600_20220929_Iw7/bwa_output/bwa_output_Final_with_headers_output'
path_to_sdmmej_del = f'{R0_POLQ_outputs}/table_outputs/bwa_output_Final_with_headers_deletions_table.csv'
sdmmej_del = pd.read_csv(path_to_sdmmej_del)

sdmmej_del['deltion_seq'] = sdmmej_del.apply(lambda row: ref[row['LEFT_DEL_INDEX']:row['RIGHT_DEL_INDEX']], axis=1)
sdmmej_del_vcf = sdmmej_del.copy()

sdmmej_del_vcf['REF'] = sdmmej_del.apply(lambda row: ref[(row['LEFT_DEL_INDEX']-1):row['RIGHT_DEL_INDEX']], axis=1)
sdmmej_del_vcf['ALT'] = sdmmej_del.apply(lambda row: ref[(row['LEFT_DEL_INDEX']-1):row['LEFT_DEL_INDEX']], axis=1)
sdmmej_del_vcf['POS'] = sdmmej_del_vcf['LEFT_DEL_INDEX']
sdmmej_del_vcf['CHR'] = [f'Iw7_D_{i}' for i in sdmmej_del_vcf.index]
sdmmej_del_vcf['del_seq'] = np.nan
sdmmej_del_vcf.rename(columns={'CONSISTENCY': 'consistency'}, inplace=True)

# Data for comparison
del_cols = ['CHR', 'POS', 'REF', 'ALT', 'RECONSTRUCTED_SEQ', 'REPAIR_TYPE', 'LEFT_DEL_INDEX',
                   'RIGHT_DEL_INDEX', 'consistency', 'CLASS', 'CLASS_final',
                   'READS', 'MICROHOMOLOGY', 'MH_Length']
sdmmej_del_vcf.loc[:, del_cols]


# insertions that was classified by HiFIBR as complex and then reclassified as ins by sdmmej
path_to_comp_ins = f'{R0_POLQ_outputs}/bwa_output_Final_with_headers_complex_insertion_consistency2.csv'
comp_ins_df = pd.read_csv(path_to_comp_ins)
comp_ins_df_vcf = comp_ins_df.copy()

comp_ins_df_vcf['ALT'] = comp_ins_df_vcf.apply(lambda row:row['RECONSTRUCTED_SEQ'][(row['left_del']-1): (row['right_del']-1)], axis=1)
comp_ins_df_vcf['REF'] = comp_ins_df_vcf.apply(lambda row:row['RECONSTRUCTED_SEQ'][(row['left_del']-1): row['left_del']], axis=1)
comp_ins_df_vcf['POS'] = comp_ins_df_vcf['left_del'] # 0 based
comp_ins_df_vcf = comp_ins_df_vcf.loc[comp_ins_df_vcf['POS'] != 0, :]
comp_ins_df_vcf['CHR'] = [f'Iw7_CI_{i}' for i in comp_ins_df_vcf.index]

# insertions that was classified by HiFIBR as insertion and then reclassified as ins by sdmmej
path_to_ins_ins = f'{R0_POLQ_outputs}/bwa_output_Final_with_headers_insertion_insertion_consistency2.csv'
ins_ins = pd.read_csv(path_to_ins_ins)

ins_ins_vcf = ins_ins.copy()
ins_ins_vcf['ALT'] = ins_ins.apply(lambda row:row['RECONSTRUCTED_SEQ'][(row['left_del']-1): (row['right_del']-1)], axis=1)
ins_ins_vcf['REF'] = ins_ins.apply(lambda row:row['RECONSTRUCTED_SEQ'][(row['left_del']-1): row['left_del']], axis=1)
ins_ins_vcf['POS'] = ins_ins['left_del'] # 0 based
ins_ins_vcf['CHR'] = [f'Iw7_II_{i}' for i in ins_ins_vcf.index]

del_cols = ['CHR', 'POS', 'REF', 'ALT','ID', 'RECONSTRUCTED_SEQ','del_seq', 'REPAIR_TYPE', 'consistency', 'CLASS','CLASS_final','READS', 'MICROHOMOLOGY', 'MH_Length']
ins_cols = ['CHR', 'POS', 'REF', 'ALT','ID', 'RECONSTRUCTED_SEQ', 'del_seq','consistency','left_del','right_del', 'DRmotif_length', 'RCmotif_length', 'Loop-out', 'Snap-back']
comb_df = pd.concat([sdmmej_del_vcf.loc[:, del_cols], 
          ins_ins_vcf.loc[:, ins_cols],
         comp_ins_df_vcf.loc[:,ins_cols]]).reset_index(drop=True)

# Exporing data
comb_df['POS'] = comb_df['POS'] + indel_pos_correction
comb_df.rename(columns={'CHR':'#CHR'}, inplace=True)

comb_df.to_csv(f'{output_path}/R0_POLQ/{date}/{date}_R0_POLQ_sdmmej_output_EMmej_input.tsv',
                index=False, sep='\t')

# Exporing the VCF format for EMmej
comb_df.loc[:, ['#CHR', 'POS', 'REF', 'ALT']].to_csv(f'{output_path}/R0_POLQ/{date}/{date}_R0_POLQ_sdmmej_output_EMmej_input.vcf',
                index=False, sep='\t')

# # ----------- Creating a 'reference genome' ----------- 
with open(f'{output_path}/R0_POLQ/{date}/{date}_R0_POLQ_sdmmej_output_EMmej_input_ref.fa', 'w') as f:
    for i in comb_df.loc[:0, :].index:
        f.write(f">\t{comb_df.loc[i, '#CHR']}\n")
        f.write(f'{pre_insert_1000}{ref[79:-79]}{post_insert_1000}\n')
    for i in comb_df.loc[1:, :].index:
        f.write(f">\t{comb_df.loc[i, '#CHR']}\n")
        f.write(f"{pre_insert_1000}{comb_df.loc[i, 'del_seq'][79:-79]}{post_insert_1000}\n")

# Writing a REDME file
README = f"""This folder contains the output of: 
/home/labs/alevy/guyta/guy_master_project/scripts/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input.py
The inputs are:
    sdmmej deletions : {path_to_sdmmej_del}
    sdmmej complex insertions: {path_to_comp_ins}
    sdmmej insertions: {path_to_ins_ins}

The outputs are:
    {date}_R0_POLQ_sdmmej_output_EMmej_input.vcf : The data in a VCF format for EMmej
    {date}_R0_POLQ_sdmmej_output_EMmej_input.tsv : The data with informarive columns for later use in the comparison
    {date}_R0_POLQ_sdmmej_output_EMmej_input_ref.fa : A reference genome for the data

"""
with open(f'{output_path}/R0_POLQ/{date}/README.md', 'w') as f:
    f.write(README)

# ------------- R0_Lig4 ------------- 
dir = os.path.join(f"{output_path}/R0_Lig4/{date}")
if not os.path.exists(dir):
    os.mkdir(dir)

# Deletions as appear in sdmmej output: Those are all deletion classified by sdmmej
R0_Lig4_outputs = '/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/R0_Lig4/raw_data_to_sdmmej_output/1626_20220929_Iw7/bwa_output/bwa_output_Final_with_headers_output'
path_to_sdmmej_del = f'{R0_Lig4_outputs}/table_outputs/bwa_output_Final_with_headers_deletions_table.csv'
sdmmej_del = pd.read_csv(path_to_sdmmej_del)

sdmmej_del['deltion_seq'] = sdmmej_del.apply(lambda row: ref[row['LEFT_DEL_INDEX']:row['RIGHT_DEL_INDEX']], axis=1)
sdmmej_del_vcf = sdmmej_del.copy()

sdmmej_del_vcf['REF'] = sdmmej_del.apply(lambda row: ref[(row['LEFT_DEL_INDEX']-1):row['RIGHT_DEL_INDEX']], axis=1)
sdmmej_del_vcf['ALT'] = sdmmej_del.apply(lambda row: ref[(row['LEFT_DEL_INDEX']-1):row['LEFT_DEL_INDEX']], axis=1)
sdmmej_del_vcf['POS'] = sdmmej_del_vcf['LEFT_DEL_INDEX']
sdmmej_del_vcf['CHR'] = [f'Iw7_D_{i}' for i in sdmmej_del_vcf.index]
sdmmej_del_vcf['del_seq'] = np.nan
sdmmej_del_vcf.rename(columns={'CONSISTENCY': 'consistency'}, inplace=True)

# Data for comparison
del_cols = ['CHR', 'POS', 'REF', 'ALT', 'RECONSTRUCTED_SEQ', 'REPAIR_TYPE', 'LEFT_DEL_INDEX',
                   'RIGHT_DEL_INDEX', 'consistency', 'CLASS', 'CLASS_final',
                   'READS', 'MICROHOMOLOGY', 'MH_Length']
sdmmej_del_vcf.loc[:, del_cols]


# insertions that was classified by HiFIBR as complex and then reclassified as ins by sdmmej
path_to_comp_ins = f'{R0_Lig4_outputs}/bwa_output_Final_with_headers_complex_insertion_consistency2.csv'
comp_ins_df = pd.read_csv(path_to_comp_ins)
comp_ins_df_vcf = comp_ins_df.copy()

comp_ins_df_vcf['ALT'] = comp_ins_df_vcf.apply(lambda row:row['RECONSTRUCTED_SEQ'][(row['left_del']-1): (row['right_del']-1)], axis=1)
comp_ins_df_vcf['REF'] = comp_ins_df_vcf.apply(lambda row:row['RECONSTRUCTED_SEQ'][(row['left_del']-1): row['left_del']], axis=1)
comp_ins_df_vcf['POS'] = comp_ins_df_vcf['left_del'] # 0 based
comp_ins_df_vcf = comp_ins_df_vcf.loc[comp_ins_df_vcf['POS'] != 0, :]
comp_ins_df_vcf['CHR'] = [f'Iw7_CI_{i}' for i in comp_ins_df_vcf.index]


# insertions that was classified by HiFIBR as insertion and then reclassified as ins by sdmmej
path_to_ins_ins = f'{R0_Lig4_outputs}/bwa_output_Final_with_headers_insertion_insertion_consistency2.csv'
ins_ins = pd.read_csv(path_to_ins_ins)

ins_ins_vcf = ins_ins.copy()
ins_ins_vcf['ALT'] = ins_ins.apply(lambda row:row['RECONSTRUCTED_SEQ'][(row['left_del']-1): (row['right_del']-1)], axis=1)
ins_ins_vcf['REF'] = ins_ins.apply(lambda row:row['RECONSTRUCTED_SEQ'][(row['left_del']-1): row['left_del']], axis=1)
ins_ins_vcf['POS'] = ins_ins['left_del'] # 0 based
ins_ins_vcf['CHR'] = [f'Iw7_II_{i}' for i in ins_ins_vcf.index]

del_cols = ['CHR', 'POS', 'REF', 'ALT','ID', 'RECONSTRUCTED_SEQ','del_seq', 'REPAIR_TYPE', 'consistency', 'CLASS','CLASS_final','READS', 'MICROHOMOLOGY', 'MH_Length']
ins_cols = ['CHR', 'POS', 'REF', 'ALT','ID', 'RECONSTRUCTED_SEQ', 'del_seq','consistency','left_del','right_del', 'DRmotif_length', 'RCmotif_length', 'Loop-out', 'Snap-back']

comb_df = pd.concat([sdmmej_del_vcf.loc[:, del_cols], 
          ins_ins_vcf.loc[:, ins_cols],
         comp_ins_df_vcf.loc[:,ins_cols]]).reset_index(drop=True)
# Exporing data
comb_df['POS'] = comb_df['POS'] + indel_pos_correction
comb_df.rename(columns={'CHR':'#CHR'}, inplace=True)

comb_df.to_csv(f'{output_path}/R0_Lig4/{date}/{date}_R0_Lig4_sdmmej_output_EMmej_input.tsv',
                index=False, sep='\t')

# Exporing the VCF format for EMmej
comb_df.loc[:, ['#CHR', 'POS', 'REF', 'ALT']].to_csv(f'{output_path}/R0_Lig4/{date}/{date}_R0_Lig4_sdmmej_output_EMmej_input.vcf',
                index=False, sep='\t')

# # ----------- Creating a 'reference genome' ----------- 
with open(f'{output_path}/R0_Lig4/{date}/{date}_R0_Lig4_sdmmej_output_EMmej_input_ref.fa', 'w') as f:
    for i in comb_df.loc[:50, :].index:
        f.write(f">\t{comb_df.loc[i, '#CHR']}\n")
        f.write(f'{pre_insert_1000}{ref[79:-79]}{post_insert_1000}\n')
    for i in comb_df.loc[51:, :].index:
        f.write(f">\t{comb_df.loc[i, '#CHR']}\n")
        f.write(f"{pre_insert_1000}{comb_df.loc[i, 'del_seq'][79:-79]}{post_insert_1000}\n")


# Writing a REDME file
README = f"""This folder contains the output of: 
/home/labs/alevy/guyta/guy_master_project/scripts/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input.py
The inputs are:
    sdmmej deletions : {path_to_sdmmej_del}
    sdmmej complex insertions: {path_to_comp_ins}
    sdmmej insertions: {path_to_ins_ins}

The outputs are:
    {date}_R0_Lig4_sdmmej_output_EMmej_input.vcf : The data in a VCF format for EMmej
    {date}_R0_Lig4_sdmmej_output_EMmej_input.tsv : The data with informarive columns for later use in the comparison
    {date}_R0_Lig4_sdmmej_output_EMmej_input_ref.fa : A reference genome for the data

"""
with open(f'{output_path}/R0_Lig4/{date}/README.md', 'w') as f:
    f.write(README)

    