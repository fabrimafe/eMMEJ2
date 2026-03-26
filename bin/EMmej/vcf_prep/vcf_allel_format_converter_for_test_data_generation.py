"""
This script converts the following vcf format:
REF     ALT
AC      A

into:
REF     ALT
A       .

This is only to create a fake dataset so I can write
a script that does the exact opposite thing
"""

import pandas as pd

path_to_df = '/home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter/GDL_Indels_regular_monoalelic_format_no_headers.vcf'
columns = ['CHR', 'POS', '.', 'REF', 'ALT'] + [i for i in range(90)]
df = pd.read_csv(path_to_df, sep='\t', names=columns)

# Transformation part
def allele_representation_converter(ref:str, alt:str):
    if len(ref) < len(alt):
        ref = '-'
        alt = alt[1:]
    if len(ref) > len(alt):
        alt = '-'
        ref = ref[1:]
    return (ref,alt)

transformed_allels = df.apply(
    lambda row: allele_representation_converter(ref=row['REF'], alt=row['ALT']), 
    axis=1, result_type='expand')
transformed_allels.rename(columns={0: 'REF_transformed', 1: 'ALT_transformed'}, inplace=True)

df.loc[:, ['REF_transformed', 'ALT_transformed']] = transformed_allels.loc[:, ['REF_transformed', 'ALT_transformed']]
df['CHR'] = 'NT_033779.5'
path_to_output = '/home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter/GDL_Indels_reformated_no_headers.vcf'
df.loc[:,['CHR', 'POS', '.', 'REF_transformed', 'ALT_transformed']+[i for i in range(90)]].to_csv(path_to_output, 
        sep='\t', index=False, header=False)