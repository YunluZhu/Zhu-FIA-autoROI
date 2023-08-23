#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.getAmp_fitDualSlope import getAmp_fitDualSlope_wBaseTrials
from functions.plt_functions import plt_categorical_grid2
# %%
root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT"
traces_avg = pd.read_hdf(f"{root}/res_concatenated.h5", key='long_data')
slope = pd.read_hdf(f"{root}/res_concatenated.h5", key='slope')
amp_long = pd.read_hdf(f"{root}/res_concatenated.h5", key='amp')


# %%
total_roi_num = len(slope.ROI_id.unique())
# QC
max_serial_diff_amp_threshold = 0.5
max_baseline_amp = 2
endRes_threshold = 0.2

peakAmp_threshold = 0.5

# 1. how does peak amp change across repeats? -------------------
# exclude ROIs that increase responses as the trial goes on
qc_repeatDiff = amp_long.groupby(['cond_num','ROI_id','nsti']).apply(
    lambda g: g['dF_ksDensity'].diff().mean()
).reset_index()
qc_repeatDiff.columns = ['cond_num', 'ROI_id', 'nsti', 'amp_trial_diff']
qc_repeatDiff = qc_repeatDiff.loc[qc_repeatDiff['cond_num'].isin([1,2])]
qc_repeatDiff_avg = qc_repeatDiff.groupby(['cond_num', 'ROI_id'])['amp_trial_diff'].mean().reset_index()
ROI_2exclude1 = qc_repeatDiff_avg.query('amp_trial_diff > @max_serial_diff_amp_threshold').ROI_id.values
ROI_passAwaken = slope.loc[~slope.ROI_id.isin(ROI_2exclude1)].ROI_id.unique()

print(f'Ratio of ROIs pass awakening QC: {len(ROI_passAwaken)/total_roi_num}')


# 2. baseline amp median -------------------
qc_amp_baseline = amp_long.groupby(['cond_num','ROI_id'])['dF_ksDensity'].median().reset_index()
qc_amp_baseline = qc_amp_baseline.loc[qc_amp_baseline['cond_num'].isin([3,4])]
qc_amp_baseline_avg = qc_amp_baseline.groupby("ROI_id")['dF_ksDensity'].mean().reset_index()
ROI_base_baseTrial = qc_amp_baseline_avg.query('dF_ksDensity < @max_baseline_amp').ROI_id.unique()

# and last3sec of each response
slope_sel_for_endRes = slope.query("cond_num == 1")
ROI_goodEnd = slope_sel_for_endRes.loc[slope_sel_for_endRes['last3sec'].abs() < endRes_threshold].ROI_id.unique()

ROI_passBase = np.intersect1d(ROI_base_baseTrial, ROI_goodEnd)

ROI_baseQC_pass = np.intersect1d(ROI_passBase, ROI_passAwaken)
print(f'Ratio of ROIs pass baseline QC: {len(ROI_baseQC_pass)/len(ROI_passAwaken)}')

# 3. peak amplitude mean across sti and trials under Dark 
amp1mean = amp_long.query("cond_num == 1").groupby(['ROI_id'])['dF_ksDensity'].mean().reset_index()
ROI_goodRes = amp1mean.query("dF_ksDensity > @peakAmp_threshold").ROI_id.values

ROI_regQC_pass = np.intersect1d(ROI_baseQC_pass, ROI_goodRes)
print(f'Ratio of ROIs pass Dark peak Amplitude QC: {len(ROI_regQC_pass)/len(ROI_baseQC_pass)}')

slope_QC = slope.loc[slope['ROI_id'].isin(ROI_regQC_pass)]


# %% filter by regression
# 1. regression reliability
r_threshold = 0.7 # 0.6-0.8
sq_residual_avg_threshold = 5
slope_threshold = 0.08

sel_fitted_forQC = slope_QC.loc[slope_QC['cond_num']==1,:]
ROI_goodFit1 = sel_fitted_forQC.query("rLowAng > @r_threshold").ROI_id.unique()

ROI_goodFit2 = sel_fitted_forQC.query("sq_residualLowAng_avg < @sq_residual_avg_threshold").ROI_id.values

ROI_goodFit = np.intersect1d(ROI_goodFit1, ROI_goodFit2)

amp_goodFit = slope_QC.loc[slope_QC['ROI_id'].isin(ROI_goodFit)]

amp_goodFit = amp_goodFit.assign(
    exp_cond_ordered = amp_goodFit['cond_num'].astype(str) + amp_goodFit['exp_cond']
)
print(f'Ratio of ROIs pass fitting QC: {len(ROI_goodFit)/len(sel_fitted_forQC.ROI_id.unique())}')

# 2. tuned or not
sel_fitted_forQC = amp_goodFit.loc[amp_goodFit['cond_num'].isin([1,2])]
sel_mean_slope = sel_fitted_forQC.groupby(['ROI_id'])['slopeLowAng'].max().reset_index()
ROI_goodTuning = sel_mean_slope.query("slopeLowAng > @slope_threshold").ROI_id.unique()
amp_goodTuning = amp_goodFit.loc[amp_goodFit['ROI_id'].isin(ROI_goodTuning)]
print(f'Ratio of ROIs with good tuning: {len(ROI_goodTuning)/len(amp_goodFit.ROI_id.unique())}')


# %% plot

toplt = amp_goodFit.loc[amp_goodFit['cond_num'].isin([1,2])]

x_name='exp_cond_ordered'
y_name='slopeLowAng'
units = 'ROI_id'

p = plt_categorical_grid2(
    gridrow='NS',
    data=toplt,
    x_name=x_name,
    y_name=y_name,
    units=units,
    height=3,
    aspect=1,
)


# %%
toplt = amp_goodTuning.loc[amp_goodTuning['cond_num'].isin([1,2])]

x_name='exp_cond_ordered'
y_name='slopeLowAng'
units = 'ROI_id'

p = plt_categorical_grid2(
    gridrow='NS',
    data=toplt,
    x_name=x_name,
    y_name=y_name,
    units=units,
    height=3,
    aspect=1,
)
# %%
