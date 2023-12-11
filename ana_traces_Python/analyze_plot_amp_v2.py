#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.plt_functions import plt_categorical_grid2
from functions.doQC_getSlope import doQC_getSlope_wTimedAvgAmp, sel_exp
from functions.doCategorization import get_peakTimingCat
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
_, amp_smval, slope = doQC_getSlope_wTimedAvgAmp(root, cond4qc)
amp_smval = amp_smval.reset_index(drop=True)
which_amp = 'amp_smval'
# slope = slope.reset_index(drop=True)
# %%

if 'lesion' in sel_dir:
    amp_smval.loc[amp_smval['cond_num']==2, 'exp_cond_ordered'] = '2lesion'
    lesion_df = amp_smval.loc[amp_smval['cond_num']==2].copy()
    values_to_norm = amp_smval.query("cond_num == 2 & stimulus == 0")['amp_smval'].values - amp_smval.query("cond_num == 1 & stimulus == 0")['amp_smval'].values
    values_to_norm = np.repeat(values_to_norm, len(lesion_df.nsti.unique()))
    cond2_base = np.repeat(lesion_df.query("stimulus == 0")['amp_smval'].values, len(lesion_df.nsti.unique()))
    amp_normed = (lesion_df['amp_smval'] - values_to_norm)/cond2_base 
    amp_smval = amp_smval.assign(
        amp_normed = amp_smval['amp_smval']
    )
    amp_smval.loc[amp_smval['cond_num']==2,'amp_normed'] = amp_normed.values
    which_amp = 'amp_normed'

# %%
df = amp_smval.loc[amp_smval[sel_qc]]

df = df.sort_values(by=['ROI_id','exp_cond_ordered']).reset_index(drop=True)

# %%
# categorize by smval peak timing
df, _ = get_peakTimingCat(df)
# %% plot
df_toplt = df

sns.relplot(
    kind='line',
    data=df_toplt.loc[df_toplt['cond_num']<3],
    x='stimulus',
    y=which_amp,
    units='ROI_id',
    estimator=None,
    row='which_exp',
    col='exp_cond_ordered',
    alpha=0.1,
    height=3
)
plt.savefig(f"{fig_dir}/rawLine {which_amp}Xstimulus.pdf", format='PDF')


sns.relplot(
    kind='line',
    data=df_toplt.loc[df_toplt['cond_num']<3],
    x='stimulus',
    y=which_amp,
    row='which_exp',
    hue='exp_cond_ordered',
    height=3
)
plt.savefig(f"{fig_dir}/avgLine {which_amp}Xstimulus.pdf", format='PDF')

g = plt_categorical_grid2(
data=df_toplt.loc[df_toplt['cond_num']<3],
gridcol='stimulus',
y_name=which_amp,
gridrow='which_exp',
x_name='exp_cond_ordered',
units='ROI_id',
alpha=0.1,
aspect=0.7
)
g.set(ylim=[-0.4, np.percentile(df_toplt[which_amp], 99.9)])
plt.savefig(f"{fig_dir}/{which_amp}Xcond.pdf", format='PDF')


#%% check timing
for which_time in ['peak_time_smval', 'half_decay_time']:
    sns.relplot(
        kind='line',
        data=df_toplt.loc[df_toplt['cond_num']<3],
        x='stimulus',
        y=which_time,
        units='ROI_id',
        estimator=None,
        row='which_exp',
        col='exp_cond_ordered',
        alpha=0.1,
        height=3
    )
    plt.savefig(f"{fig_dir}/rawLine {which_amp}Xstimulus.pdf", format='PDF')


    sns.relplot(
        kind='line',
        data=df_toplt.loc[df_toplt['cond_num']<3],
        x='stimulus',
        y=which_time,
        row='which_exp',
        hue='exp_cond_ordered',
        height=3
    )
    plt.savefig(f"{fig_dir}/avgLine {which_amp}Xstimulus.pdf", format='PDF')


# %%


# %%
df_toplt = df
for which_amp in ['amp_raw', 'amp_smval']:
    df_toplt = df_toplt.sort_values(by=['ROI_id','cond_num','nsti']).reset_index(drop=True)

    cond1_df = df_toplt.query("cond_num == 1")
    cond2_df = df_toplt.query("cond_num == 2")
    df_change = cond1_df.copy()
    df_change = df_change.assign(
        amp_chg = cond2_df[which_amp].values - cond1_df[which_amp].values,
        amp_chg_ratio = (cond2_df[which_amp].values - cond1_df[which_amp].values)/ cond1_df[which_amp].values,
        amp_chg_norm = (cond2_df[which_amp].values - cond1_df[which_amp].values)/ (cond2_df[which_amp].values + cond1_df[which_amp].values),
    )
    df_change = df_change.query('nsti != 0')
    exclude_for_plotting1 = df_change.loc[df_change['amp_chg_norm'].abs() > 1].ROI_id.unique()
    exclude_for_plotting2 = cond1_df.loc[(cond1_df['nsti'] > 0) & (cond1_df[which_amp] < 0)].ROI_id.unique()
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
        plt.savefig(f"{fig_dir}/{which_amp}Chg {y_name}_{x_name}.pdf", format='PDF')

#%%
df_toplt = df
###### separate by timing ########
for which_amp in ['amp_raw', 'amp_smval']:
    sns.relplot(
        kind='line',
        data=df_toplt.loc[(df_toplt['cond_num'].isin([1,2]))],
        x='stimulus',
        y=which_amp,
        row='which_exp',
        col='peakTiming_cat',
        hue='exp_cond_ordered',
        height=3,
        # units='ROI_id',
        # estimator=None,
        # alpha=0.2
    )
    plt.savefig(f"{fig_dir}/amp by peak time.pdf", format='PDF')


# %%

# %%
