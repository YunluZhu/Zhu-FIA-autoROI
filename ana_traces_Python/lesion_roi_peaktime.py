#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.plt_functions import plt_categorical_grid2
from functions.doQC_getSlope import doQC_getSlope_4lesion
import matplotlib.pyplot as plt
from sklearn import preprocessing as pre

#%%

root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_lesion"

#%%
amp_QC, amp_goodFit, amp_goodTuning = doQC_getSlope_4lesion(root)

# %%
STIMULUS_EXT = [0,5,10,20,30]
fig_dir = f"{root}/figures_timing"
try:
    os.makedirs(fig_dir)
except:
    pass

amp_QC = amp_QC.assign(
    exp_cond_ordered = amp_QC['cond_num'].astype(str) + amp_QC['exp_cond']
)
df = amp_goodTuning

df = df.sort_values(by=['fish_id', 'ROI','exp_cond_ordered']).reset_index(drop=True)
# df = df.assign(
#     peak_cat_cond1 = df.groupby(['ROI_id'])['peak_cat'].transform(lambda g: g.iloc[0])
# )
df.rename(columns={'last3sec':'amp_0'}, inplace=True)
df_toplt = pd.wide_to_long(df, stubnames='amp', i=['ROI_id', 'cond_num'], j='nsti', sep='_').reset_index()
sti_map = dict([(ii, sti) for ii, sti  in enumerate(STIMULUS_EXT)])
df_toplt = df_toplt.assign(
    stimulus = df_toplt['nsti'].map(sti_map),
)

df_toplt.loc[df_toplt['cond_num']==2, 'exp_cond_ordered'] = '2lesion'

#%%
df_toplt = df_toplt.assign(
    amp_normed = df_toplt.groupby(['which_exp', 'ROI_id','cond_num'])['amp'].transform(
        lambda x: (pre.MinMaxScaler().fit_transform(x.values.reshape(-1,1))).ravel().tolist()
        ),
    peak_cat_sti_peakTime = df_toplt.groupby(['ROI_id','cond_num'])['peak_time'].transform(
        'mean'
    )
)

df_toplt = df_toplt.assign(
    peak_cat_avg = pd.cut(df_toplt['peak_cat_sti_peakTime'], bins=[-1,0.4,3])
)
# %%

sns.relplot(
    kind='line',
    data=df_toplt.loc[df_toplt['cond_num'].isin([1,2])],
    x='stimulus',
    y='amp',
    row='which_exp',
    col='peak_cat_avg',
    hue='exp_cond_ordered',
    height=3
)
plt.savefig(f"{fig_dir}/goodTuning avgLine ampXstimulus.pdf", format='PDF')

# %%
