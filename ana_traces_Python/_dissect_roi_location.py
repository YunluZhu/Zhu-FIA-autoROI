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

amp_QC, amp_goodFit, amp_goodTuning = doQC_getSlope_4lesion(root)
ROI_metadata = pd.read_hdf(f'{root}/res_concatenated.h5', key='roi_metadata')

STIMULUS_EXT = [0,5,10,20,30]
fig_dir = f"{root}/figures"
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

df_toplt = df_toplt.sort_values(by=['ROI_id','cond_num','nsti']).reset_index(drop=True)

cond1_df = df_toplt.query("cond_num == 1")
cond2_df = df_toplt.query("cond_num == 2")
df_change = cond1_df.copy()
df_change = df_change.assign(
    amp_chg = cond2_df['amp'].values - cond1_df['amp'].values,
    amp_chg_ratio = (cond2_df['amp'].values - cond1_df['amp'].values)/ cond1_df['amp'].values,
    amp_chg_norm = (cond2_df['amp'].values - cond1_df['amp'].values)/ (cond2_df['amp'].values + cond1_df['amp'].values),
)
df_change = df_change.query('nsti != 0')
exclude_for_plotting1 = df_change.loc[df_change['amp_chg_norm'].abs() > 1].ROI_id.unique()
exclude_for_plotting2 = cond1_df.loc[(cond1_df['nsti'] > 0) & (cond1_df['amp'] < 0)].ROI_id.unique()
exclude_for_plotting = np.union1d(exclude_for_plotting1, exclude_for_plotting2)

df_change = df_change.loc[~df_change['ROI_id'].isin(exclude_for_plotting)]
# %%
ROI_metadata = ROI_metadata.assign(
    ROI_id = ROI_metadata['fish_id'] + '_' + ROI_metadata['id'].astype(str)
)
# %%



# anatomical?
df_change = df_change.assign(
    amp_chg_cat = pd.cut(df_change['amp_chg_norm'], bins=[-1, -0.1, 0.1, 1], labels = ['reduced', 'noChg', 'increased'])
)
amp_chg_mata = df_change.merge(ROI_metadata, on='ROI_id')
# %
sns.scatterplot(
    data=amp_chg_mata.query("nsti==4"),
    hue='amp_chg_cat',
    x='xCenter',
    y='yCenter'
)
plt.axis('equal')

# %%

# %%

# %%
