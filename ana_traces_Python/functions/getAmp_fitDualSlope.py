'''
FINISHED
'''

    # %%
import pandas as pd
import numpy as np
import seaborn as sns
import os
import scipy.io
import matplotlib.pyplot as plt
import itertools
import plotly.express as px
import plotly.graph_objs as go
import time
import random

def pairwise_correlation(A, B):
    am = A - np.mean(A, axis=0, keepdims=True)
    bm = B - np.mean(B, axis=0, keepdims=True)
    return am.T @ bm /  (np.sqrt(
        np.sum(am**2, axis=0,
               keepdims=True)).T * np.sqrt(
        np.sum(bm**2, axis=0, keepdims=True)))
    
def np_pearson_cor(x, y):
    xv = x - x.mean(axis=0)
    yv = y - y.mean(axis=0)
    xvss = (xv * xv).sum(axis=0)
    yvss = (yv * yv).sum(axis=0)
    result = np.matmul(xv.transpose(), yv) / np.sqrt(np.outer(xvss, yvss))
    # bound the values to -1 to 1 in the event of precision issues
    return np.maximum(np.minimum(result, 1.0), -1.0)

    # %%

def getAmp_fitDualSlope_kdeBaseCond1base(root, STIMULUS=[5, 10, 20, 30], if_shuffle=False, shuffle_cond=[1,2]):
    # %%
    DATA_FILE='dF_ksDensity'
    fig_dir = f"{root}/figures"
    
    try:
        os.makedirs(fig_dir)
    except:
        pass
        
    nsti = len(STIMULUS)

    raw_df = pd.read_csv(os.path.join(root, "rawF_df.csv"))
    mat = scipy.io.loadmat(os.path.join(root, DATA_FILE+'.mat'))[DATA_FILE]
    traces = pd.DataFrame(data=mat)

    # find baseline file
    baseline_file = 'baseline' + '_' + DATA_FILE.split('_')[1]
    if DATA_FILE.endswith('adj'):
        baseline_file = baseline_file + '_adj'
    mat = scipy.io.loadmat(os.path.join(root, baseline_file+'.mat'))

    baseline_file_name = list(mat.keys())[-1]
    baseline_wide = pd.DataFrame(data=mat[baseline_file_name])
    baseline_cond1 = baseline_wide.iloc[0,:].values
    dFF = pd.DataFrame(data=np.divide(traces.values, baseline_cond1))
    
    # baseline_diff = baseline_wide.diff().iloc[1,:].T
    # ROI_resonable_baseline = list(baseline_diff.loc[baseline_diff>0].index + 1)

    if len(raw_df.columns) > len(dFF.columns):
        col_diff = len(raw_df.columns) - len(dFF.columns) 
        dFF = pd.concat([dFF, raw_df.iloc[:,-col_diff:]], axis=1)

    dFF.columns = [sub.replace('rawF', 'roi_') for sub in raw_df.columns]

    nreps = dFF.groupby(['area','repeat']).ngroups
    trial_frames = int(len(raw_df)/nreps)

    dFF = dFF.assign(
        frames = list(np.arange(trial_frames)) * nreps
    )
    
    if if_shuffle:
        # SHUFFLE area repeats together
        dFF_1 = dFF.loc[dFF['area'].isin(shuffle_cond)]
        dFF_2 = dFF.loc[~dFF['area'].isin(shuffle_cond)]
        groups = [dFF_1 for _, dFF_1 in dFF_1.groupby(['area','repeat','frames'])]
        random.shuffle(groups)
        dFF1_sh = pd.concat(groups).reset_index(drop=True)
        dFF_1 = dFF_1.assign(
            area = dFF1_sh.area.values,
            repeat = dFF1_sh.repeat.values,
            frames = dFF1_sh.frames.values,
        )
        dFF = pd.concat([dFF_1, dFF_2], ignore_index=True)
    
    # find fish info
    fish_info = root.split("/")[-1].split(" ")
    exp_date, fishNum = fish_info[0].split("_")
    exp_date = int(exp_date)
    fishNum = fishNum.split("fish")[1]
    fishNum = int(fishNum)

    print(fish_info)

    # find area names (which is also saved in metadata files)
    area_cond_dict = {}

    for file in os.listdir(root):
        d = os.path.join(root, file)
        if os.path.isdir(d):
            if str.startswith(file, "area"):
                area, cond = file.split("_")
                area_number = area.split("area",1)[1]
                area_cond_dict[area_number] = cond

    # this is cumbersome, but read fish metadata to get frame rate

    ROI_metadata = pd.read_csv(os.path.join(root, "ROI_metadata.csv"))

    fish_metaDATA_FILE = []
    parent_dir = os.path.abspath(os.path.join(root, os.pardir))
    for file in os.listdir(parent_dir):
        if str.endswith(file, " metadata.csv"):
            fish_metaDATA_FILE.append(os.path.join(parent_dir, file))

    frame_rate = pd.Series(dtype='float64')
    fish_metaDATA_FILE.sort()
    for fish_metaDATA_FILE_sel in fish_metaDATA_FILE:
        if frame_rate.empty:
            fish_metadata = pd.read_csv(fish_metaDATA_FILE_sel)
            try:
                frame_rate = fish_metadata.loc[
                    (fish_metadata['exp_date']==exp_date) & (fish_metadata['fish_num']==fishNum), 'frame_rate'
                    ]
            except:
                pass
        else:
            break

    frame_rate = frame_rate.iloc[0]
    vol_rate = frame_rate / (ROI_metadata['zPos'].max() + 1)
    time_stamp = np.arange(0,20,1/vol_rate)[1:]

    sti_frames = trial_frames / nsti

    #%%
    dFF_long = pd.wide_to_long(dFF.reset_index(), stubnames='roi', i='index', j='ROI', sep='_').reset_index()
    dFF_long = dFF_long.rename(columns={'roi':DATA_FILE})
    dFF_long = dFF_long.assign(
        exp_cond = dFF_long['area'].astype(str).map(area_cond_dict),
        fish_id = fish_info[0],
        fish_info = fish_info[1],
        time = list(time_stamp) * nsti * nreps * len(ROI_metadata),
        nsti = list(np.repeat(np.arange(nsti)+1, sti_frames)) * nreps * len(ROI_metadata)
    )

    ROI_metadata = ROI_metadata.assign(
        fish_id = fish_info[0],
        fish_info = fish_info[1],
    )

    # %%
    # clean up ROI
    # Get paires of ROIs that are on adjacent slices and are close with eachother in x y

    # to do: save a new column into ROI metadata: roi_cond containing standarized description of the condition of each ROI.
    # read that via metadata file
    # map to the response result dataframe
    ROI_RESPONSE_CORR_THRESHOLD = 0.8
    
    # print(f"- {len(ROI_resonable_baseline)}/{len(ROI_metadata)} ROIs with reasonable baseline")

    dFF_long_roi_corrected = dFF_long#.loc[dFF_long['ROI'].isin(ROI_resonable_baseline)]

    n_combined = 1000
    all_combined = 0
    
    while n_combined > 0:
        unique_roi = dFF_long_roi_corrected['ROI'].unique()
        similar_roi = []
        
        if n_combined == 1000: # if first run. go through all ROIs
            for id1, id2 in itertools.combinations(unique_roi, 2):
                roi1 = ROI_metadata.query("id == @id1").iloc[0]
                roi2 = ROI_metadata.query("id == @id2").iloc[0]
                if_adjacent_z = abs(roi1.zPos - roi2.zPos) <= 1
                if_close_x = abs(roi1.xCenter - roi2.xCenter) <= 2
                if_close_y = abs(roi1.yCenter - roi2.yCenter) <= 2
                if_spacially_same = if_adjacent_z & if_close_x & if_close_y
                if if_spacially_same:
                    similar_roi.append((id1, id2))
        else: # otherwise, recalculate ROIs that have been changed
            for id1, id2 in itertools.combinations(unique_roi, 2):
                # only do the rois that have been modified last run
                if (id1 in roi_modified_lastRun) | (id2 in roi_modified_lastRun):
                    roi1 = ROI_metadata.query("id == @id1").iloc[0]
                    roi2 = ROI_metadata.query("id == @id2").iloc[0]
                    if_adjacent_z = abs(roi1.zPos - roi2.zPos) <= 1
                    if_close_x = abs(roi1.xCenter - roi2.xCenter) <= 2
                    if_close_y = abs(roi1.yCenter - roi2.yCenter) <= 2
                    if_spacially_same = if_adjacent_z & if_close_x & if_close_y
                    if if_spacially_same:
                        similar_roi.append((id1, id2))
                
        n_combined = 0
                
        r_out = []   
        roi_modified_lastRun = []
        # response correlation
        for id1, id2 in similar_roi:
            res1 = dFF_long_roi_corrected.query("ROI == @id1")
            res2 = dFF_long_roi_corrected.query("ROI == @id2")
            # this_corr = np.corrcoef(res1[DATA_FILE].values, res2[DATA_FILE].values)
            # this_r = this_corr[0,1]
            this_corr = pairwise_correlation(res1[DATA_FILE].values, res2[DATA_FILE].values)
            this_r = this_corr[0]
            # this_corr = np_pearson_cor(res1[DATA_FILE].values, res2[DATA_FILE].values)
            # this_r = this_corr[0,0]
            
            r_out.append(this_r)

        for i, this_r in enumerate(r_out):
            id1, id2 = similar_roi[i]
            if this_r > ROI_RESPONSE_CORR_THRESHOLD:  # if two nearby ROIs have similar response, combine them
                unique_roi = dFF_long_roi_corrected['ROI'].unique()
                if (id1 in unique_roi) & (id2 in unique_roi):
                    n_combined += 1
                    res1 = dFF_long_roi_corrected.query("ROI == @id1")
                    res2 = dFF_long_roi_corrected.query("ROI == @id2")
                    res_combined = res1.copy()
                    res_combined[DATA_FILE] = (res1[DATA_FILE].values + res2[DATA_FILE].values)/2
                    dFF_long_roi_corrected.loc[dFF_long_roi_corrected['ROI']==id1, DATA_FILE] = res_combined[DATA_FILE]
                    dFF_long_roi_corrected = dFF_long_roi_corrected.drop(index=res2.index)
                    roi_modified_lastRun.append(id1)
        all_combined += n_combined
        
    # end_time = time.time()
    # print(end_time-start_time)
    print(f'- {all_combined} ROIs combined')




    # %%
    # Extract peaks and make a new dataframe of peaks
    # the dFF_long_roi_corrected should contain time column
   # find the idx of time of response during a certain winder after each sti
    amp_cols = ['amp_' + str(sti) for sti in np.arange(nsti)+1]
    time_for_amp = 4 #sec
    frame_for_amp = int(np.floor(time_for_amp * vol_rate))

    sel_unique_roi_level = ['area','repeat','ROI']
    idx_for_amp = dFF_long_roi_corrected.groupby(sel_unique_roi_level + ['nsti']).head(frame_for_amp).index
    rows_for_amp = dFF_long_roi_corrected.loc[idx_for_amp]
    amp_long = rows_for_amp.groupby(sel_unique_roi_level + ['nsti'])[DATA_FILE].max().reset_index()

    # %%
    # QC

    dFF_THRESHOLD = 0.3  # mean peak amp of dFF
    # dFF_stdMean_THRESHOLD = 10
    # BASE_THRESHOLD = 0.8  # std/mean
    LAST3s_THRESHOLD = 1 # dFF

    # 1. find ROIs with mean peak amp across all stimulus greater than threshold 
    amp_long_ampfilter = amp_long.loc[amp_long['area'].isin([1,2])].groupby(['area','ROI', 'nsti'])[DATA_FILE].median().reset_index() # mdian across repeats for each sti of each ROI
    amp_long_ampfilter_max = amp_long_ampfilter.groupby(['area','ROI'])[DATA_FILE].max().reset_index() # max of the amp to all sti should pass threshold
    ROI_pass1 = amp_long_ampfilter_max.loc[amp_long_ampfilter_max[DATA_FILE] > dFF_THRESHOLD].ROI.unique()
    print(f"- {len(ROI_pass1)}/{len(dFF_long_roi_corrected.ROI.unique())} ROIs passed dFF threshold")

    # 2. find ROIs with relatively stable responses across repeats
    # amp_stability = amp_long.query("area == 1").groupby(['ROI','nsti'])[DATA_FILE].apply(
    #     lambda x: x.std()/x.mean()
    # ).reset_index()
    # amp_stability = amp_stability.groupby('ROI').mean().reset_index()
    # amp_stability_th = np.percentile(amp_stability[DATA_FILE],95)
    # ROI_pass2 = amp_stability.loc[amp_stability[DATA_FILE]<amp_stability_th].ROI.unique()
    # print(f"- {len(ROI_pass2)}/{len(ROI_metadata)} ROIs passed dFF amp stability threshold")


    # # 3. find ROIs with response down to baseline at the end of each sti
    # idx_for_last3sec = dFF_long_roi_corrected.groupby(sel_unique_roi_level + ['nsti']).tail(int(np.floor(3 * vol_rate))).index
    # rows_for_last3sec = dFF_long_roi_corrected.loc[idx_for_last3sec]
    # last3sec_long = rows_for_last3sec.groupby(sel_unique_roi_level + ['nsti'])[DATA_FILE].mean(numeric_only=True).reset_index() # mean of the values of last 3 sec
    # last3sec_long = last3sec_long.query("area == 1")
    # last3sec_long_groupMean = last3sec_long.groupby(['ROI','area'])[DATA_FILE].median() # median across repeats
    # ROI_pass3 = last3sec_long_groupMean.loc[last3sec_long_groupMean < LAST3s_THRESHOLD].reset_index().ROI.unique()
    # print(f"- {len(ROI_pass3)}/{len(ROI_metadata)} ROIs passed baseline threshold")


    # ROI_passQC = np.intersect1d(ROI_pass1, ROI_pass3)
    # ROI_passQC = np.intersect1d(ROI_pass12, ROI_pass3)
    ROI_passQC = ROI_pass1

    # print(f"> {len(ROI_passQC)} passed QC")
    amp_selected = amp_long.loc[amp_long['ROI'].isin(ROI_passQC)]
    dFF_long_roi_selected = dFF_long_roi_corrected.loc[dFF_long_roi_corrected['ROI'].isin(ROI_passQC)]

    ############# look at different raw traces #############
    # peak time
    dFF_long_roi_averaged = dFF_long_roi_selected.groupby(['frames','exp_cond','fish_id','nsti','area','ROI',])[DATA_FILE].median().reset_index()
    df_peak_time_perRep = dFF_long_roi_selected.groupby(['fish_id','exp_cond','area','ROI','nsti','repeat']).head(int(5*vol_rate)).groupby(['fish_id','exp_cond','area','ROI','nsti','repeat'])[DATA_FILE].apply(np.argmax).reset_index()
    df_peak_time = df_peak_time_perRep.groupby(['fish_id','exp_cond','area','ROI','nsti'])[DATA_FILE].mean().reset_index()
    df_peak_time = df_peak_time.assign(
        peak_time = df_peak_time[DATA_FILE]/vol_rate
    )

    # categorize responses by their peak time in selected stimulus
    df_peak_time_cat = df_peak_time.loc[df_peak_time.nsti.isin([2,3,4])].groupby(['fish_id','exp_cond','area','ROI'])['peak_time'].median().reset_index()
    df_peak_time_cat = df_peak_time_cat.assign(
        peak_cat = pd.cut(df_peak_time_cat.peak_time, bins=[-np.inf, 1.5, 5.5, np.inf], labels=['fast', 'slow', 'late'])
    )

    df_peak_time_cat_area1 = df_peak_time_cat.query("area==1")

    res_traces_avg = dFF_long_roi_averaged.merge(df_peak_time_cat[['ROI','fish_id','peak_cat','peak_time','area']], on=['fish_id', 'ROI','area'])
    res_traces_avg = res_traces_avg.assign(
        fish_id = fish_info[0],
        fish_info = fish_info[1],
    )
    # dFF_long_roi_averaged = dFF_long_roi_averaged.sort_values(by=['fish_id','area','ROI','nsti','frames'])
    # %%
    # change of response to individual stimuli

    ############### % check raw traces ##############
    df_toplt = dFF_long_roi_averaged.merge(df_peak_time_cat_area1[['ROI','fish_id','peak_cat','peak_time']], on=['fish_id', 'ROI'])
    df_toplt.peak_cat = df_toplt.peak_cat.astype(str)
    g = sns.relplot(row='exp_cond', col='peak_cat', kind='line', data=df_toplt, x='frames',y=DATA_FILE,units='ROI', estimator=None, alpha=0.1, aspect=3, height=1.8,
                    col_order=['fast','slow']
                    )
    g.set(ylim=[-1,12])
    plt.savefig(f"{fig_dir}/traces.pdf", format='PDF')
    plt.close()


    # %% 
    amp_wide = pd.pivot(
        amp_selected,
        index=sel_unique_roi_level,
        columns='nsti',
        values=DATA_FILE,
    ).reset_index()

    amp_wide.columns.name = None
    amp_wide.columns = sel_unique_roi_level + amp_cols
    amp_wide = amp_wide.rename(columns={'area':'cond_num'})
    # amp_wide = amp_wide.merge(ROI_metadata, how='left', left_on='ROI', right_on='id')

    amp_wide = amp_wide.assign(
        exp_cond = amp_wide['cond_num'].astype(str).map(area_cond_dict),
    )

    idx_for_last3sec = dFF_long_roi_corrected.groupby(sel_unique_roi_level + ['nsti']).tail(int(np.floor(3 * vol_rate))).index
    rows_for_last3sec = dFF_long_roi_corrected.loc[idx_for_last3sec]
    rows_for_last3sec = rows_for_last3sec.loc[rows_for_last3sec['nsti'].isin([1,2])]
    last3sec_avg = rows_for_last3sec.loc[rows_for_last3sec['ROI'].isin(ROI_passQC)].groupby(['ROI','repeat','area'])[DATA_FILE].median()
    last3sec_avg.name = 'last3sec' 
    last3sec_avg = last3sec_avg.reset_index()
    
    amp_wide = amp_wide.merge(last3sec_avg, left_on=['cond_num', 'repeat', 'ROI'], right_on=['area', 'repeat', 'ROI'])


# %%
    ################## SLOPE on raw amplitudes #################################

    slopeAll = []
    sq_residual = []
    corr = []
    sq_residual_avg = []

    slopeLowAng = []
    sq_residualLowAng = []
    sq_residualLowAng_avg = []
    corrLowAng = []

    X = np.array(STIMULUS)
    X = np.append(X, [0])
    XlowAng = np.append(X[0:2], [0])
    amp_colslowAng = amp_cols[0:2]+ ['last3sec']

    for (cond_num, roi), group in amp_wide.groupby(['cond_num', 'ROI']):
        g_num_of_repeat = len(group.repeat)
        g_X = np.vstack((np.ones(len(X)), X)).T
        g_X = np.tile(g_X, (g_num_of_repeat, 1))
        g_Y = np.float64(group[amp_cols + ['last3sec']].values).flatten()
        b = (np.linalg.inv(g_X.T @ g_X) @ (g_X.T) @ g_Y.T)
        slopeAll.append(b[1])
        e = np.sum(np.square(g_Y - g_X.dot(b)))
        sq_residual.append(e)
        sq_residual_avg.append(e/g_num_of_repeat)
        r = pairwise_correlation(g_X[:,1], g_Y)
        corr.append(r[0])
        
        g_XlowAng = np.vstack((np.ones(len(XlowAng)), XlowAng)).T
        g_XlowAng = np.tile(g_XlowAng, (g_num_of_repeat, 1))
        g_YlowAng = np.float64(group[amp_colslowAng]).flatten()
        b = (np.linalg.inv(g_XlowAng.T @ g_XlowAng) @ (g_XlowAng.T) @ g_YlowAng.T)
        slopeLowAng.append(b[1])
        e = np.sum(np.square(g_YlowAng - g_XlowAng.dot(b)))
        sq_residualLowAng.append(e)
        sq_residualLowAng_avg.append(e/g_num_of_repeat)
        r = pairwise_correlation(g_XlowAng[:,1], g_YlowAng)
        corrLowAng.append(r[0])

    amp_avg = amp_wide.groupby(['cond_num', 'ROI','exp_cond'])[amp_cols + ['last3sec']].median().reset_index()
    amp_avg = amp_avg.assign(
        slopeAll_rawAmp = slopeAll,
        sq_residual = sq_residual,
        sq_residual_avg = sq_residual_avg,
        r = corr,
        slopeLowAng = slopeLowAng,
        sq_residualLowAng = sq_residualLowAng,
        sq_residualLowAng_avg = sq_residualLowAng_avg,
        rLowAng = corrLowAng,
        fish_id = fish_info[0],
        fish_info = fish_info[1],
    )

    res_slope = amp_avg.merge(df_peak_time_cat[['ROI', 'peak_time', 'area', 'peak_cat']], left_on=['ROI','cond_num'], right_on=['ROI','area'])


    # %% look at raw data
    sti_map = dict([(ii+1, sti) for ii, sti  in enumerate(STIMULUS)])

    res_amp_long = amp_selected.assign(
        stimulus = amp_selected['nsti'].map(sti_map),
        exp_cond = amp_selected['area'].astype(str).map(area_cond_dict),
        fish_id = fish_info[0],
        fish_info = fish_info[1],

    ).rename(columns={'area':'cond_num'})

    # %%
    filename = 'dFF_analyzed.h5'
    
    if if_shuffle:
        filename = 'dFF_shuffled.h5'
        
    res_traces_avg.rename(columns={DATA_FILE:'dFF'}, inplace=True)
    res_amp_long.rename(columns={DATA_FILE:'dFF'}, inplace=True)
    res_slope.rename(columns={DATA_FILE:'dFF'}, inplace=True)

    res_amp_long.to_hdf(f'{root}/{filename}', key='amp', mode='w', format='table')
    res_traces_avg.to_hdf(f'{root}/{filename}', key='traces', format='table')
    ROI_metadata.to_hdf(f'{root}/{filename}', key='roi_metadata', format='table')
    res_slope.to_hdf(f'{root}/{filename}', key='slope', format='table')

    return res_traces_avg, res_amp_long, res_slope, ROI_metadata, STIMULUS
