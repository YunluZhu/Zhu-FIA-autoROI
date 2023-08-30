
#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.doQC_getSlope import doQC_getSlope_wTimedAvgAmp
from scipy.signal import savgol_filter


#%%

def get_peakTimingCat(df):
    '''
    get peak timing category, assign a new column
    '''
    df_oneTiming = df.groupby(['ROI_id']).apply(
        lambda x: x.loc[(x['cond_num'].isin([1])) & (x['nsti'].isin([2,3])), 'peak_time_smval'].mean()
    )
    df_oneTiming_cat = pd.cut(df_oneTiming, bins=[-1, 1, 100], labels = ['1fast', '2slow']).to_dict()
    df = df.assign(
        peakTiming_cat = df['ROI_id'].map(df_oneTiming_cat)
    )
    return df, 'peakTiming_cat'


#%%


