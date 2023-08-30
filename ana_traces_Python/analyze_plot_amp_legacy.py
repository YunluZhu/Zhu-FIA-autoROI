#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.plt_functions import plt_categorical_grid2
from functions.doQC_getSlope import doQC_getSlope, sel_exp
import matplotlib.pyplot as plt

#%%
sel_dir = 'lesion' # lesion or light
sel_qc = 'good_tuning'

#%%
fig_root = f"/Users/yunluzhu/Documents/Lab2/caiman/Volumetric_code/YZ_nMLF_speed/figures"
fig_folder_name = sel_dir
STIMULUS_EXT = [0,5,10,20,30]
fig_dir = os.path.join(fig_root, fig_folder_name)
try:
    os.makedirs(fig_dir)
except:
    pass

root, cond4qc = sel_exp(sel_dir)

# %%
_, slope = doQC_getSlope(root, cond4qc)
slope = slope.reset_index(drop=True)
# %%

if 'lesion' in sel_dir:
    slope.loc[slope['cond_num']==2, 'exp_cond_ordered'] = '2lesion'
    lesion_df = slope.loc[slope['cond_num']==2].copy()
    values_to_norm = slope.query("cond_num == 2")['last3sec'].values - slope.query("cond_num == 1")['last3sec'].values
    
    for amp in ['amp_1', 'amp_2', 'amp_3', 'amp_4','last3sec']:
        normed_values = (lesion_df[amp].values - values_to_norm)/slope.query("cond_num == 2")['last3sec'].values
        lesion_df.loc[:,amp] = normed_values
    lesion_df = lesion_df.assign(
        cond_num = 2.5,
        exp_cond = 'postProxBi_adj',
        area = 2.5,
        exp_cond_ordered = '3lesion_adj'
    )
    slope = pd.concat([slope.query("exp_cond != 'postProxBi_adj'"), lesion_df], ignore_index=True)

# %%
df = slope.loc[slope[sel_qc]]

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
df_toplt = df_toplt.loc[df_toplt['amp'].abs() < 20]
# %% plot
sns.relplot(
    kind='line',
    data=df_toplt.loc[df_toplt['cond_num']<3],
    x='stimulus',
    y='amp',
    units='ROI_id',
    estimator=None,
    row='which_exp',
    col='exp_cond_ordered',
    alpha=0.1,
    height=3
)
plt.savefig(f"{fig_dir}/rawLine ampXstimulus.pdf", format='PDF')


sns.relplot(
    kind='line',
    data=df_toplt.loc[df_toplt['cond_num']<3],
    x='stimulus',
    y='amp',
    row='which_exp',
    hue='exp_cond_ordered',
    height=3
)
plt.savefig(f"{fig_dir}/avgLine ampXstimulus.pdf", format='PDF')

# %%
g = plt_categorical_grid2(
    data=df_toplt.loc[df_toplt['cond_num']<3],
    gridcol='stimulus',
    y_name='amp',
    gridrow='which_exp',
    x_name='exp_cond_ordered',
    units='ROI_id',
    alpha=0.1,
    aspect=0.7
)
g.set(ylim=[-0.4, np.percentile(df_toplt.amp, 99.9)])
plt.savefig(f"{fig_dir}/ampXcond.pdf", format='PDF')

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
plt.savefig(f"{fig_dir}/lowAngSel ampXcond.pdf", format='PDF')

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
        g.set(ylim=[np.percentile(df_change[y_name], 0.02),np.percentile(df_change[y_name], 99)])
    plt.savefig(f"{fig_dir}/AmpChg {y_name}_{x_name}.pdf", format='PDF')

# %%
###### separate by timing ########
sns.relplot(
    kind='line',
    data=df_toplt.loc[(df_toplt['cond_num'].isin([1,2]))],
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

# %%
