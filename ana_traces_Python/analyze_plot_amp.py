#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from plot_functions.plt_functions import plt_categorical_grid2
from doQC_getSlope import doQC_getSlope
import matplotlib.pyplot as plt

#%%

root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT"

#%%
amp_QC, amp_goodFit, amp_goodTuning = doQC_getSlope(root)

# %%
STIMULUS = [0,5,10,20,30]
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
sti_map = dict([(ii, sti) for ii, sti  in enumerate(STIMULUS)])
df_toplt = df_toplt.assign(
    stimulus = df_toplt['nsti'].map(sti_map),
)
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

dark_df = df_toplt.query("exp_cond_ordered == '1dark'")
light_df = df_toplt.query("exp_cond_ordered == '2light'")
df_change = dark_df.copy()
df_change = df_change.assign(
    amp_chg = light_df['amp'].values - dark_df['amp'].values,
    amp_chg_ratio = (light_df['amp'].values - dark_df['amp'].values)/ dark_df['amp'].values,
    amp_chg_norm = (light_df['amp'].values - dark_df['amp'].values)/ (light_df['amp'].values + dark_df['amp'].values),
)
df_change = df_change.query('nsti != 0')

y_name='amp_chg_norm'
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
# g.set(ylim=[-1,5])
plt.savefig(f"{fig_dir}/AmpChg {y_name}_{x_name}.pdf", format='PDF')

# %%
