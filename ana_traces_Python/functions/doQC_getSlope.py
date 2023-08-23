
#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.plt_functions import plt_categorical_grid2

#%%

# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT"

def doQC_getSlope_4LD(root):
    # only use condition 1 & 2 for quality control bc cond 3 & 4 have different baseline activities

    amp_long = pd.read_hdf(f'{root}/res_concatenated.h5', key='amp')
    traces_avg = pd.read_hdf(f'{root}/res_concatenated.h5', key='long_data')
    ROI_metadata = pd.read_hdf(f'{root}/res_concatenated.h5', key='roi_metadata')
    slope = pd.read_hdf(f'{root}/res_concatenated.h5', key='slope')

    # %%
    all_roi = slope.ROI_id.unique()
    # QC
    max_serial_diff_amp_threshold = 0.3
    # max_baseline_amp = 5
    peakAmp_threshold = 0.5
    endRes_threshold = 0.7
    r_threshold = 0.5
    # sq_residual_avg_threshold = 100
    slope_threshold = 0.1

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
    
    return amp_QC, amp_goodFit, amp_goodTuning


def doQC_getSlope_4lesion(root):
    # only use condition 1 for quality control bc cond 2 is lesion
    
    amp_long = pd.read_hdf(f'{root}/res_concatenated.h5', key='amp')
    traces_avg = pd.read_hdf(f'{root}/res_concatenated.h5', key='long_data')
    ROI_metadata = pd.read_hdf(f'{root}/res_concatenated.h5', key='roi_metadata')
    slope = pd.read_hdf(f'{root}/res_concatenated.h5', key='slope')

    all_roi = slope.ROI_id.unique()
    # QC

    max_serial_diff_amp_threshold = 0.3
    # max_baseline_amp = 5
    peakAmp_threshold = 0.5
    endRes_threshold = 0.7
    r_threshold = 0.5
    # sq_residual_avg_threshold = 100
    slope_threshold = 0.1


    # 1. how does peak amp change across repeats? -------------------
    # exclude ROIs that increase responses as the trial goes on
    qc_repeatDiff = amp_long.groupby(['cond_num','ROI_id','nsti']).apply(
        lambda g: g['dFF'].diff().mean()
    ).reset_index()
    qc_repeatDiff.columns = ['cond_num', 'ROI_id', 'nsti', 'amp_trial_diff']
    qc_repeatDiff = qc_repeatDiff.loc[qc_repeatDiff['cond_num'].isin([1])]
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
    slope_sel_for_endRes = slope.loc[slope['cond_num'].isin([1])].groupby(['ROI_id'])['last3sec'].median().reset_index()
    ROI_goodEnd = slope_sel_for_endRes.loc[slope_sel_for_endRes['last3sec'] < endRes_threshold].ROI_id.unique()

    ROI_afterBase = np.intersect1d(ROI_passAwaken, ROI_goodEnd)
    print(f'Ratio of ROIs pass baseline QC: {len(ROI_afterBase)}/{len(ROI_passAwaken)}')

    # 3. peak amplitude mean across sti and trials under Dark and Light across ALL sti and repeat -----------
    amp1mean = amp_long.loc[amp_long["cond_num"].isin([1])].groupby(['ROI_id','nsti','cond_num'])['dFF'].median().reset_index()
    amp1meanmax = amp1mean.groupby(['ROI_id'])['dFF'].max().reset_index()
    ROI_goodRes = amp1meanmax.query("dFF > @peakAmp_threshold").ROI_id.values

    ROI_regQC_pass = np.intersect1d(ROI_afterBase, ROI_goodRes)
    print(f'Ratio of ROIs pass peak Amplitude QC: {len(ROI_regQC_pass)}/{len(ROI_afterBase)}')

    amp_QC = slope.loc[slope['ROI_id'].isin(ROI_regQC_pass)]

    # % filter by regression -------
    # 1. regression reliability

    sel_fitted_forQC1 = amp_QC.loc[amp_QC['cond_num'].isin([1]),:].groupby('ROI_id')['r'].min().reset_index()
    ROI_goodFit1 = sel_fitted_forQC1.query("r > @r_threshold").ROI_id.unique()

    # sel_fitted_forQC2 = amp_QC.loc[amp_QC['cond_num'].isin([1,2]),:].groupby('ROI_id')['sq_residualLowAng_avg'].max().reset_index()
    # ROI_goodFit2 = sel_fitted_forQC2.query("sq_residualLowAng_avg < @sq_residual_avg_threshold").ROI_id.values

    ROI_goodFit = ROI_goodFit1

    amp_goodFit = amp_QC.loc[amp_QC['ROI_id'].isin(ROI_goodFit)]

    amp_goodFit = amp_goodFit.assign(
        exp_cond_ordered = amp_goodFit['cond_num'].astype(str) + amp_goodFit['exp_cond']
    )
    print(f'Ratio of ROIs pass regression QC: {len(ROI_goodFit)}/{len(amp_QC.ROI_id.unique())}')

    # 2. tuned or not
    sel_fitted_forQC = amp_goodFit.loc[amp_goodFit['cond_num'].isin([1]),:]
    ROI_passTuning = sel_fitted_forQC.groupby(['ROI_id'])['slopeAll_rawAmp'].max() > slope_threshold
    ROI_goodTuning = ROI_passTuning.loc[ROI_passTuning].index

    amp_goodTuning = amp_goodFit.loc[amp_goodFit['ROI_id'].isin(ROI_goodTuning)]
    print(f'Ratio of ROIs are tuned: {len(ROI_goodTuning)}/{len(amp_goodFit.ROI_id.unique())}')

    return amp_QC, amp_goodFit, amp_goodTuning
