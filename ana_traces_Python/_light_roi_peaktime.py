#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.plt_functions import plt_categorical_grid2
from functions.doQC_getSlope import doQC_getSlope_4LD, doQC_getSlope_4lesion
import matplotlib.pyplot as plt

#%%

root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_light"

#%%
all_amp_long, amp_QC, amp_goodFit, amp_goodTuning = doQC_getSlope_4LD(root)

# %%
STIMULUS_EXT = [0,5,10,20,30]
fig_root = f"/Users/yunluzhu/Documents/Lab2/caiman/Volumetric_code/YZ_nMLF_speed/figures"
fig_folder_name = "lightPeakTime"
fig_dir = os.path.join(fig_root, fig_folder_name)

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

# %%
# timing

df = all_amp_long#.loc[all_amp_long['ROI_id'].isin(amp_goodTuning.ROI_id.unique())]
df = df.assign(
    exp_cond_ordered = df['cond_num'].astype(str) + df['exp_cond']
)
df = df.sort_values(by=['ROI_id','exp_cond_ordered','nsti','repeat']).reset_index(drop=True)

df_oneTimePersti = df.groupby(['which_exp','ROI_id','cond_num','exp_cond_ordered','nsti'])[['peak_time_onTrialAvg','peak_time_onAvgTrial']].mean().reset_index()

df_oneTimePersti = df_oneTimePersti.assign(
    peak_time_adj = (df_oneTimePersti['peak_time_onTrialAvg']+df_oneTimePersti['peak_time_onAvgTrial'])/2
)
#%% plot timing per sti
sns.catplot(
    kind='point',
    data=df_oneTimePersti,
    x='nsti',
    y='peak_time_adj',
    row='which_exp',
    col='exp_cond_ordered',
    hue='exp_cond_ordered',
    height=3,
)
plt.savefig(f"{fig_dir}/peak timing.pdf", format='PDF')

#%% use sti 10,20 for timing calculation

df_oneTimePersti_sel = df_oneTimePersti.loc[df_oneTimePersti['nsti'].isin([2,3])]

# determine peak timing category
df_oneTimePersti_sel_oneTimePerROI = df_oneTimePersti_sel.groupby(['ROI_id','cond_num', 'nsti'])['peak_time_adj'].mean().reset_index()
df_oneTimePersti_sel_oneTimePerROI = df_oneTimePersti_sel_oneTimePerROI.loc[df_oneTimePersti_sel_oneTimePerROI['cond_num'].isin([1,2])]
df_oneTimePersti_sel_oneTimePerROI = df_oneTimePersti_sel_oneTimePerROI.groupby(['ROI_id', 'nsti'])['peak_time_adj'].max().reset_index()
ROI_id = df_oneTimePersti_sel_oneTimePerROI['ROI_id'].values
ROI_cat = pd.cut(df_oneTimePersti_sel_oneTimePerROI['peak_time_adj'], bins=[-1,0.5,5], labels=['1slow','2fast']).values

ROI_cat_mapper = dict(zip(ROI_id, ROI_cat))
#%%
df_toplt = df_toplt.loc[df_toplt['cond_num'].isin([1,2])]

df_toplt = df_toplt.assign(
    peak_cat = df_toplt['ROI_id'].map(ROI_cat_mapper)
)

df_toplt = df_toplt.sort_values(by=['ROI_id','exp_cond_ordered','stimulus','peak_cat']).reset_index(drop=True)
# %%

sns.relplot(
    kind='line',
    data=df_toplt,
    x='stimulus',
    y='amp',
    row='which_exp',
    col='peak_cat',
    hue='exp_cond_ordered',
    height=3,
    # units='ROI_id',
    # estimator=None,
    # alpha=0.2
)
plt.savefig(f"{fig_dir}/amp by peak time.pdf", format='PDF')

# %%
