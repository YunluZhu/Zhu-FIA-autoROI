#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.plt_functions import plt_categorical_grid2
import matplotlib.pyplot as plt

#%%

root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_light"


#%%
amp_long = pd.read_hdf(f'{root}/res_shuffled.h5', key='amp')
traces_avg = pd.read_hdf(f'{root}/res_shuffled.h5', key='long_data')
ROI_metadata = pd.read_hdf(f'{root}/res_shuffled.h5', key='roi_metadata')
slope = pd.read_hdf(f'{root}/res_shuffled.h5', key='slope')

# %%
all_roi = slope.ROI_id.unique()
# QC
max_serial_diff_amp_threshold = 0.3
# max_baseline_amp = 5
peakAmp_threshold = 0.5
endRes_threshold = 0.7
r_threshold = 0.3
# sq_residual_avg_threshold = 100
slope_threshold = 0.08

# 1. how does peak amp change across repeats? -------------------
# exclude ROIs that increase responses as the trial goes on
qc_repeatDiff = amp_long.groupby(['cond_num','ROI_id','nsti']).apply(
    lambda g: g['dFF'].diff().mean()
).reset_index()
qc_repeatDiff.columns = ['cond_num', 'ROI_id', 'nsti', 'amp_trial_diff']
qc_repeatDiff = qc_repeatDiff.loc[qc_repeatDiff['cond_num'].isin([1,2])]
qc_repeatIncrease = qc_repeatDiff.groupby(['ROI_id','cond_num']).apply(
    lambda x: (x['amp_trial_diff'] > 0).sum() == len(x)
).reset_index()
qc_repeatIncrease.columns = ['ROI_id', 'cond_num', 'if_increase']
ROI_increaseAmp = qc_repeatIncrease.loc[qc_repeatIncrease['if_increase']].ROI_id.unique()
qc_repeatDiffToFilter = qc_repeatDiff.loc[qc_repeatDiff.ROI_id.isin(ROI_increaseAmp)].groupby(['cond_num','ROI_id'])['amp_trial_diff'].mean().reset_index()
ROI_2exclude = qc_repeatDiffToFilter.query("amp_trial_diff > @max_serial_diff_amp_threshold").ROI_id.unique()
ROI_passAwaken = list(all_roi)
for ele in ROI_2exclude:
    ROI_passAwaken.remove(ele)

print(f'Ratio of ROIs pass awakening QC: {len(ROI_passAwaken)}/{len(all_roi)}')


# 2. baseline amp median -------------------
# qc_amp_baseline = amp_long.groupby(['cond_num','ROI_id'])['dFF'].median().reset_index()
# qc_amp_baseline = qc_amp_baseline.loc[qc_amp_baseline['cond_num'].isin([3,4])]
# qc_amp_baseline = qc_amp_baseline.groupby(['ROI_id'])['dFF'].mean().reset_index()
# ROI_base_baseTrial_2exclude = qc_amp_baseline.query('dFF > @max_baseline_amp').ROI_id.unique()

# ROI_passBaseTrial = ROI_passAwaken.copy()
# for ele in np.intersect1d(ROI_base_baseTrial_2exclude, ROI_passBaseTrial):
#     ROI_passBaseTrial.remove(ele)

# and last3sec of each response
slope_sel_for_endRes = slope.loc[slope['cond_num'].isin([1,2])].groupby(['ROI_id'])['last3sec'].median().reset_index()
ROI_goodEnd = slope_sel_for_endRes.loc[slope_sel_for_endRes['last3sec'] < endRes_threshold].ROI_id.unique()

ROI_afterBase = np.intersect1d(ROI_passAwaken, ROI_goodEnd)
print(f'Ratio of ROIs pass baseline QC: {len(ROI_afterBase)}/{len(ROI_passAwaken)}')

# 3. peak amplitude mean across sti and trials under Dark and Light across ALL sti and repeat -----------
amp1mean = amp_long.loc[amp_long["cond_num"].isin([1,2])].groupby(['ROI_id','nsti','cond_num'])['dFF'].median().reset_index()
amp1meanmax = amp1mean.groupby(['ROI_id'])['dFF'].max().reset_index()
ROI_goodRes = amp1meanmax.query("dFF > @peakAmp_threshold").ROI_id.values

ROI_regQC_pass = np.intersect1d(ROI_afterBase, ROI_goodRes)
print(f'Ratio of ROIs pass peak Amplitude QC: {len(ROI_regQC_pass)}/{len(ROI_afterBase)}')

amp_QC = slope.loc[slope['ROI_id'].isin(ROI_regQC_pass)]

# % filter by regression -------
# 1. regression reliability

sel_fitted_forQC1 = amp_QC.loc[amp_QC['cond_num'].isin([1,2]),:].groupby('ROI_id')['r'].min().reset_index()
ROI_goodFit1 = sel_fitted_forQC1.query("r > @r_threshold").ROI_id.unique()

ROI_goodFit = ROI_goodFit1

amp_goodFit = amp_QC.loc[amp_QC['ROI_id'].isin(ROI_goodFit)]

amp_goodFit = amp_goodFit.assign(
    exp_cond_ordered = amp_goodFit['cond_num'].astype(str) + amp_goodFit['exp_cond']
)
print(f'Ratio of ROIs pass regression QC: {len(ROI_goodFit)}/{len(amp_QC.ROI_id.unique())}')

# 2. tuned or not
sel_fitted_forQC = amp_goodFit.loc[amp_goodFit['cond_num'].isin([1, 2]),:]
ROI_passTuning = sel_fitted_forQC.groupby(['ROI_id'])['slopeAll_rawAmp'].max() > slope_threshold
ROI_goodTuning = ROI_passTuning.loc[ROI_passTuning].index

amp_goodTuning = amp_goodFit.loc[amp_goodFit['ROI_id'].isin(ROI_goodTuning)]
print(f'Ratio of ROIs are tuned: {len(ROI_goodTuning)}/{len(amp_goodFit.ROI_id.unique())}')

#%% ---------------------PLOTTING STARTS HERE------------------
STIMULUS_EXT = [0,5,10,20,30]
fig_root = f"/Users/yunluzhu/Documents/Lab2/caiman/Volumetric_code/YZ_nMLF_speed/figures"
fig_folder_name = "lightShuffle"
fig_dir = os.path.join(fig_root, fig_folder_name)

try:
    os.makedirs(fig_dir)
except:
    pass

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


# %%  ------ plot ---------

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
exclude_for_plotting1 = df_change.loc[df_change['amp_chg_norm'].abs() > 1].ROI_id.unique()
exclude_for_plotting2 = dark_df.loc[(dark_df['nsti'] > 0) & (dark_df['amp'] < 0)].ROI_id.unique()
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


#%% ------ peak time ---------


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

df = amp_long.loc[amp_long['ROI_id'].isin(amp_goodTuning.ROI_id.unique())]
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
df_oneTimePersti_sel_oneTimePerROI = df_oneTimePersti_sel_oneTimePerROI.groupby(['ROI_id', 'nsti'])['peak_time_adj'].min().reset_index()
ROI_id = df_oneTimePersti_sel_oneTimePerROI['ROI_id'].values
ROI_cat = pd.cut(df_oneTimePersti_sel_oneTimePerROI['peak_time_adj'], bins=[-1,.5,5], labels=['1fast','2slow']).values

ROI_cat_mapper = dict(zip(ROI_id, ROI_cat))
#%%
df_toplt = df_toplt.loc[df_toplt['cond_num'].isin([1,2])]

df_toplt = df_toplt.assign(
    peak_cat = df_toplt['ROI_id'].map(ROI_cat_mapper)
)

df_toplt = df_toplt.sort_values(by=['peak_cat','ROI_id','exp_cond_ordered','stimulus']).reset_index(drop=True)
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
