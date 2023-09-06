
#%%
import os,glob
import pandas as pd
import numpy as np 


#%%

def sel_exp(sel_exp:str):
    exp_dict = {
        #'name of exp': ("dir", [condition number for quality control])
        'light': ("/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_light", [1,2]),
        'lesion': ("/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_lesion", [1]),
        'local': ("/Users/yunluzhu/Documents/Lab2/Manuscripts/2023 nMLF Speed Navigation/data/data_copy_localUse/analyzed_light", [1,2])
    }
    path, cond_4qc = exp_dict[sel_exp]
    return path, cond_4qc

#%%
def doQC_getSlope(root:str, sel_cond_4qc:list):
    #%%
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
    grp = amp_long.groupby(['cond_num','fish_id','ROI_id','nsti'])
    dFF_diff = grp['dFF'].transform(
        lambda g: g.diff()
    )
    qc_repeatDiff_all = amp_long.assign(
        amp_trial_diff = dFF_diff
    )
    qc_repeatDiff = qc_repeatDiff_all.groupby(['cond_num','fish_id','ROI_id','repeat'])['amp_trial_diff'].mean().reset_index()
    qc_repeatDiff.dropna(inplace=True)
    qc_repeatDiff = qc_repeatDiff.loc[qc_repeatDiff['cond_num'].isin(sel_cond_4qc)]
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
    slope_sel_for_endRes = slope.loc[slope['cond_num'].isin(sel_cond_4qc)].groupby(['ROI_id'])['last3sec'].median().reset_index()
    ROI_goodEnd = slope_sel_for_endRes.loc[slope_sel_for_endRes['last3sec'] < endRes_threshold].ROI_id.unique()

    ROI_afterBase = np.intersect1d(ROI_passAwaken, ROI_goodEnd)
    print(f'Ratio of ROIs pass baseline QC: {len(ROI_afterBase)}/{len(ROI_passAwaken)}')

    # 3. peak amplitude mean across sti and trials 
    amp1mean = amp_long.loc[amp_long["cond_num"].isin(sel_cond_4qc)].groupby(['ROI_id','nsti','cond_num'])['dFF'].median().reset_index()
    amp1meanmax = amp1mean.groupby(['ROI_id'])['dFF'].max().reset_index()
    ROI_goodRes = amp1meanmax.query("dFF > @peakAmp_threshold").ROI_id.values

    ROI_regQC_pass = np.intersect1d(ROI_afterBase, ROI_goodRes)
    print(f'Ratio of ROIs pass peak Amplitude QC: {len(ROI_regQC_pass)}/{len(ROI_afterBase)}')

    amp_QC = slope.loc[slope['ROI_id'].isin(ROI_regQC_pass)]

    # % filter by regression -------
    # 1. regression reliability

    sel_fitted_forQC1 = amp_QC.loc[amp_QC['cond_num'].isin(sel_cond_4qc),:].groupby('ROI_id')['r'].min().reset_index()
    ROI_goodFit1 = sel_fitted_forQC1.query("r > @r_threshold").ROI_id.unique()

    # sel_fitted_forQC2 = amp_QC.loc[amp_QC['cond_num'].isin(sel_cond_4qc),:].groupby('ROI_id')['sq_residualLowAng_avg'].max().reset_index()
    # ROI_goodFit2 = sel_fitted_forQC2.query("sq_residualLowAng_avg < @sq_residual_avg_threshold").ROI_id.values

    ROI_goodFit = ROI_goodFit1
    amp_goodFit = amp_QC.loc[amp_QC['ROI_id'].isin(ROI_goodFit)]

    print(f'Ratio of ROIs pass regression QC: {len(ROI_goodFit)}/{len(amp_QC.ROI_id.unique())}')

    # 2. tuned or not
    sel_fitted_forQC = amp_goodFit.loc[amp_goodFit['cond_num'].isin(sel_cond_4qc),:]
    ROI_passTuning = sel_fitted_forQC.groupby(['ROI_id'])['slopeAll_rawAmp'].max() > slope_threshold
    ROI_goodTuning = ROI_passTuning.loc[ROI_passTuning].index

    print(f'Ratio of ROIs are tuned: {len(ROI_goodTuning)}/{len(ROI_goodFit)}')
    
    amp_long = amp_long.assign(
        exp_cond_ordered = amp_long['cond_num'].astype(str) + amp_long['exp_cond'],
        good_fit = False,
        good_tuning = False,
    )
    amp_long.loc[amp_long['ROI_id'].isin(ROI_goodFit),'good_fit'] = True
    amp_long.loc[amp_long['ROI_id'].isin(ROI_goodTuning),'good_tuning'] = True
    
    amp_QC = amp_QC.assign(
        exp_cond_ordered = amp_QC['cond_num'].astype(str) + amp_QC['exp_cond'],
        good_fit = False,
        good_tuning = False,
    )
    amp_QC.loc[amp_QC['ROI_id'].isin(ROI_goodFit),'good_fit'] = True
    amp_QC.loc[amp_QC['ROI_id'].isin(ROI_goodTuning),'good_tuning'] = True
    
    ########## PEAK TIMING ##############
    df = amp_long#.loc[all_amp_long['ROI_id'].isin(amp_goodTuning.ROI_id.unique())]

    df = amp_long.sort_values(by=['ROI_id','exp_cond_ordered','nsti','repeat']).reset_index(drop=True)

    df_oneTimePersti = df.groupby(['which_exp','ROI_id','cond_num','exp_cond_ordered','nsti'])[['peak_time_onTrialAvg','peak_time_onAvgTrial']].mean().reset_index()

    df_oneTimePersti = df_oneTimePersti.assign(
        peak_time_adj = (df_oneTimePersti['peak_time_onTrialAvg']+df_oneTimePersti['peak_time_onAvgTrial'])/2
    )

    #% use sti 10,20 for timing calculation ---- 

    df_oneTimePersti_sel = df_oneTimePersti.loc[df_oneTimePersti['nsti'].isin([2,3])]

    # determine peak timing category
    df_oneTimePersti_sel_oneTimePerROI = df_oneTimePersti_sel.groupby(['ROI_id','cond_num', 'nsti'])['peak_time_adj'].mean().reset_index()
    df_oneTimePersti_sel_oneTimePerROI = df_oneTimePersti_sel_oneTimePerROI.loc[df_oneTimePersti_sel_oneTimePerROI['cond_num'].isin(sel_cond_4qc)]
    df_oneTimePersti_sel_oneTimePerROI = df_oneTimePersti_sel_oneTimePerROI.groupby(['ROI_id', 'nsti'])['peak_time_adj'].max().reset_index()
    ROI_id = df_oneTimePersti_sel_oneTimePerROI['ROI_id'].values
    ROI_cat = pd.cut(df_oneTimePersti_sel_oneTimePerROI['peak_time_adj'], bins=[-1,0.5,5], labels=['1slow','2fast']).values

    ROI_cat_mapper = dict(zip(ROI_id, ROI_cat))

    amp_long = amp_long.assign(
        peak_cat = amp_long['ROI_id'].map(ROI_cat_mapper)
    )
    
    amp_QC = amp_QC.assign(
        peak_cat = amp_QC['ROI_id'].map(ROI_cat_mapper)
    )
    slope_QC = amp_QC
    
    return amp_long, slope_QC

#%%
def doQC_getSlope_wTimedAvgAmp(root:str, sel_cond_4qc:list):
    #%%
    amp_long = pd.read_hdf(f'{root}/res_concatenated.h5', key='amp')
    amp_avg = pd.read_hdf(f'{root}/res_concatenated.h5', key='amp_avg')
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
    grp = amp_long.groupby(['cond_num','fish_id','ROI_id','nsti'])
    dFF_diff = grp['dFF'].transform(
        lambda g: g.diff()
    )
    qc_repeatDiff_all = amp_long.assign(
        amp_trial_diff = dFF_diff
    )
    qc_repeatDiff = qc_repeatDiff_all.groupby(['cond_num','fish_id','ROI_id','repeat'])['amp_trial_diff'].mean().reset_index()
    qc_repeatDiff.dropna(inplace=True)
    qc_repeatDiff = qc_repeatDiff.loc[qc_repeatDiff['cond_num'].isin(sel_cond_4qc)]
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
    slope_sel_for_endRes = slope.loc[slope['cond_num'].isin(sel_cond_4qc)].groupby(['ROI_id'])['last3sec'].median().reset_index()
    ROI_goodEnd = slope_sel_for_endRes.loc[slope_sel_for_endRes['last3sec'] < endRes_threshold].ROI_id.unique()

    ROI_afterBase = np.intersect1d(ROI_passAwaken, ROI_goodEnd)
    print(f'Ratio of ROIs pass baseline QC: {len(ROI_afterBase)}/{len(ROI_passAwaken)}')

    # 3. peak amplitude mean across sti and trials 
    amp1mean = amp_long.loc[amp_long["cond_num"].isin(sel_cond_4qc)].groupby(['ROI_id','nsti','cond_num'])['dFF'].median().reset_index()
    amp1meanmax = amp1mean.groupby(['ROI_id'])['dFF'].max().reset_index()
    ROI_goodRes = amp1meanmax.query("dFF > @peakAmp_threshold").ROI_id.values

    ROI_regQC_pass = np.intersect1d(ROI_afterBase, ROI_goodRes)
    print(f'Ratio of ROIs pass peak Amplitude QC: {len(ROI_regQC_pass)}/{len(ROI_afterBase)}')

    amp_QC = slope.loc[slope['ROI_id'].isin(ROI_regQC_pass)]

    # % filter by regression -------
    # 1. regression reliability

    sel_fitted_forQC1 = amp_QC.loc[amp_QC['cond_num'].isin(sel_cond_4qc),:].groupby('ROI_id')['r'].min().reset_index()
    ROI_goodFit1 = sel_fitted_forQC1.query("r > @r_threshold").ROI_id.unique()

    # sel_fitted_forQC2 = amp_QC.loc[amp_QC['cond_num'].isin(sel_cond_4qc),:].groupby('ROI_id')['sq_residualLowAng_avg'].max().reset_index()
    # ROI_goodFit2 = sel_fitted_forQC2.query("sq_residualLowAng_avg < @sq_residual_avg_threshold").ROI_id.values

    ROI_goodFit = ROI_goodFit1
    amp_goodFit = amp_QC.loc[amp_QC['ROI_id'].isin(ROI_goodFit)]

    print(f'Ratio of ROIs pass regression QC: {len(ROI_goodFit)}/{len(amp_QC.ROI_id.unique())}')

    # 2. tuned or not
    sel_fitted_forQC = amp_goodFit.loc[amp_goodFit['cond_num'].isin(sel_cond_4qc),:]
    ROI_passTuning = sel_fitted_forQC.groupby(['ROI_id'])['slopeAll_rawAmp'].max() > slope_threshold
    ROI_goodTuning = ROI_passTuning.loc[ROI_passTuning].index

    print(f'Ratio of ROIs are tuned: {len(ROI_goodTuning)}/{len(ROI_goodFit)}')
    
    amp_long = amp_long.assign(
        exp_cond_ordered = amp_long['cond_num'].astype(str) + amp_long['exp_cond'],
        good_fit = False,
        good_tuning = False,
    )
    amp_long.loc[amp_long['ROI_id'].isin(ROI_goodFit),'good_fit'] = True
    amp_long.loc[amp_long['ROI_id'].isin(ROI_goodTuning),'good_tuning'] = True
    
    amp_QC = amp_QC.assign(
        exp_cond_ordered = amp_QC['cond_num'].astype(str) + amp_QC['exp_cond'],
        good_fit = False,
        good_tuning = False,
    )
    amp_QC.loc[amp_QC['ROI_id'].isin(ROI_goodFit),'good_fit'] = True
    amp_QC.loc[amp_QC['ROI_id'].isin(ROI_goodTuning),'good_tuning'] = True
    
    ########## PEAK TIMING ##############
    df = amp_long#.loc[all_amp_long['ROI_id'].isin(amp_goodTuning.ROI_id.unique())]

    df = amp_long.sort_values(by=['ROI_id','exp_cond_ordered','nsti','repeat']).reset_index(drop=True)

    df_oneTimePersti = df.groupby(['which_exp','ROI_id','cond_num','exp_cond_ordered','nsti'])[['peak_time_onTrialAvg','peak_time_onAvgTrial']].mean().reset_index()

    df_oneTimePersti = df_oneTimePersti.assign(
        peak_time_adj = (df_oneTimePersti['peak_time_onTrialAvg']+df_oneTimePersti['peak_time_onAvgTrial'])/2
    )

    #%% use sti 10,20 for timing calculation ---- 

    df_oneTimePersti_sel = df_oneTimePersti.loc[df_oneTimePersti['nsti'].isin([2,3])]

    # determine peak timing category
    df_oneTimePersti_sel_oneTimePerROI = df_oneTimePersti_sel.groupby(['ROI_id','cond_num', 'nsti'])['peak_time_adj'].mean().reset_index()
    df_oneTimePersti_sel_oneTimePerROI = df_oneTimePersti_sel_oneTimePerROI.loc[df_oneTimePersti_sel_oneTimePerROI['cond_num'].isin(sel_cond_4qc)]
    df_oneTimePersti_sel_oneTimePerROI = df_oneTimePersti_sel_oneTimePerROI.groupby(['ROI_id', 'nsti'])['peak_time_adj'].max().reset_index()
    ROI_id = df_oneTimePersti_sel_oneTimePerROI['ROI_id'].values
    ROI_cat = pd.cut(df_oneTimePersti_sel_oneTimePerROI['peak_time_adj'], bins=[-1,0.5,5], labels=['1slow','2fast']).values

    ROI_cat_mapper = dict(zip(ROI_id, ROI_cat))

    amp_long = amp_long.assign(
        peak_cat = amp_long['ROI_id'].map(ROI_cat_mapper)
    )
    
    amp_QC = amp_QC.assign(
        peak_cat = amp_QC['ROI_id'].map(ROI_cat_mapper)
    )
    slope_QC = amp_QC
    
    #%% use amp_avg dataframe which contains peak calculated on smoothed median traces and timing for peak and half decay
    amp_smval = amp_avg.assign(
        exp_cond_ordered = amp_avg['cond_num'].astype(str) + amp_avg['exp_cond'],
        good_awakened = False,
        good_endRes = False,
        good_fit = False,
        good_tuning = False,
        
    )
    amp_smval.loc[amp_smval['ROI_id'].isin(ROI_passAwaken),'good_awakened'] = True
    amp_smval.loc[amp_smval['ROI_id'].isin(ROI_afterBase),'good_endRes'] = True
    amp_smval.loc[amp_smval['ROI_id'].isin(ROI_goodFit),'good_fit'] = True
    amp_smval.loc[amp_smval['ROI_id'].isin(ROI_goodTuning),'good_tuning'] = True
    
    
    return amp_long, amp_smval, slope_QC