# 24.03.2022
"""
This is a test file to make sure that the BlockBootstraper
work as expected
"""

# import libreries
import sys
sys.path.append('/home/labs/alevy/guyta/guy_master_project/src')

import pandas as pd
import numpy as np

from data_exploration_util import *
from BlockBootstraper import *
mem_reporter()
# setting pandas display options
pd.options.display.max_colwidth = 2200
pd.set_option("display.max_columns", None)

# creating a fake dataframe
df_size = 10
fake_df = pd.DataFrame(columns=['coordinate', 'feature','value'],
                     index=[i for i in range(0,df_size)])
fake_df.loc[:,'coordinate'] = [i for i in range(0,(df_size*2),2)]

print(fake_df)
