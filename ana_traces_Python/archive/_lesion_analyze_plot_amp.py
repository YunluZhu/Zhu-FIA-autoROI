#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.plt_functions import plt_categorical_grid2
from functions.doQC_getSlope import doQC_getSlope_4lesion
import matplotlib.pyplot as plt

#%%

root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_lesion"

#%%
amp_QC, amp_goodFit, amp_goodTuning = doQC_getSlope_4lesion(root)

# %%
STIMULUS_EXT = [0,5,10,20,30]
fig_root = f"/Users/yunluzhu/Documents/Lab2/caiman/Volumetric_code/YZ_nMLF_speed/figures"
fig_folder_name = "lesionAmp"
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

df_toplt.loc[df_toplt['cond_num']==2, 'exp_cond_ordered'] = '2lesion'
# %% plot
sns.relplot(
    kind='line',
    data=df_toplt.loc[df_toplt['cond_num'].isin([1,2])],
    x='stimulus',
    y='amp',
    units='ROI_id',
    estimator=None,
    row='which_exp',
    col='exp_cond_ordered',
    alpha=0.1,
    height=3
)
plt.savefig(f"{fig_dir}/goodTuning rawLine ampXstimulus.pdf", format='PDF')


sns.relplot(
    kind='line',
    data=df_toplt.loc[df_toplt['cond_num'].isin([1,2])],
    x='stimulus',
    y='amp',
    row='which_exp',
    hue='exp_cond_ordered',
    height=3
)
plt.savefig(f"{fig_dir}/goodTuning avgLine ampXstimulus.pdf", format='PDF')


# %%

# %%
g = plt_categorical_grid2(
    data=df_toplt.loc[df_toplt['cond_num'].isin([1,2])],
    gridcol='stimulus',
    y_name='amp',
    gridrow='which_exp',
    x_name='exp_cond_ordered',
    units='ROI_id',
    alpha=0.1,
    aspect=0.7
)
g.set(ylim=[-0.4, np.percentile(df_toplt.amp, 99.9)])
plt.savefig(f"{fig_dir}/goodTuning ampXcond.pdf", format='PDF')

# %%
g = plt_categorical_grid2(
    data=df_toplt.loc[(df_toplt['cond_num'].isin([1,2])) & (df_toplt['nsti'].isin([1,2,3])),:],
    gridcol='stimulus',
    y_name='amp',
    gridrow='which_exp',
    x_name='exp_cond_ordered',
    units='ROI_id',
    alpha=0.1,
    aspect=0.7
)
plt.savefig(f"{fig_dir}/goodTuning lowAngSel ampXcond.pdf", format='PDF')

# %%
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

for y_name in ['amp_chg', 'amp_chg_ratio', 'amp_chg_norm']:
    x_name='which_exp'

    g = plt_categorical_grid2(
        data=df_change,
        y_name=y_name,
        x_name=x_name,
        gridcol='stimulus',
        # gridrow='stimulus',
        units='ROI_id',
        aspect=0.7
    )
    if y_name != 'amp_chg_norm':
        g.set(ylim=[np.percentile(df_change[y_name], 0.1),np.percentile(df_change[y_name], 99)])
    plt.savefig(f"{fig_dir}/AmpChg {y_name}_{x_name}.pdf", format='PDF')

# %%
