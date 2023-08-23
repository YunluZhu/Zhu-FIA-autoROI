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

# %%
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
root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_lesion/230622_fish4 prox_bi"
STIMULUS=[5, 10, 20, 30]

# %%
# read dF data 
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

# use cond 1 baseline for division
baseline_cond1 = baseline_wide.iloc[0,:].values
dFF = pd.DataFrame(data=np.divide(traces.values, baseline_cond1))

if len(raw_df.columns) > len(dFF.columns):
    col_diff = len(raw_df.columns) - len(dFF.columns) 
    dFF = pd.concat([dFF, raw_df.iloc[:,-col_diff:]], axis=1)

dFF.columns = [sub.replace('rawF', 'roi_') for sub in raw_df.columns]

nreps = dFF.groupby(['area','repeat']).ngroups
trial_frames = int(len(raw_df)/nreps)

dFF = dFF.assign(
    frames = list(np.arange(trial_frames)) * nreps
)

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

dFF_long_roi_corrected = dFF_long

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
time_for_amp = 3 #sec
frame_for_amp = int(np.floor(time_for_amp * vol_rate))

sel_unique_roi_level = ['area','repeat','ROI']
idx_for_amp = dFF_long_roi_corrected.groupby(sel_unique_roi_level + ['nsti']).head(frame_for_amp).index
rows_for_amp = dFF_long_roi_corrected.loc[idx_for_amp]
amp_long = rows_for_amp.groupby(sel_unique_roi_level + ['nsti'])[DATA_FILE].max().reset_index()

# %%
# QC

dFF_THRESHOLD = 0.5  # mean peak amp of dFF
# dFF_stdMean_THRESHOLD = 10
# BASE_THRESHOLD = 0.8  # std/mean
LAST3s_THRESHOLD = 1 # dFF

# 1. find ROIs with mean peak amp across all stimulus greater than threshold 
amp_long_cond1_group_median = amp_long.query("area == 1").groupby(sel_unique_roi_level + ['nsti'])[DATA_FILE].median().reset_index() # mdian across repeats for each sti of each ROI
amp_long_cond1_group_max = amp_long_cond1_group_median.groupby(sel_unique_roi_level)[DATA_FILE].max().reset_index() # max of the amp to all sti should pass threshold
ROI_pass1 = amp_long_cond1_group_max.loc[amp_long_cond1_group_max[DATA_FILE] > dFF_THRESHOLD].ROI.unique()
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
df_peak_time = dFF_long_roi_averaged.groupby(['fish_id','exp_cond','area','ROI','nsti']).head(int(5*vol_rate)).groupby(['fish_id','exp_cond','area','ROI','nsti'])[DATA_FILE].apply(np.argmax).reset_index()
df_peak_time = df_peak_time.assign(
    peak_time = df_peak_time[DATA_FILE]/vol_rate
)

# categorize responses by their peak time in selected two stimulus
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
# plt.close()


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

res_traces_avg.rename(columns={DATA_FILE:'dFF'}, inplace=True)
res_amp_long.rename(columns={DATA_FILE:'dFF'}, inplace=True)
res_slope.rename(columns={DATA_FILE:'dFF'}, inplace=True)

res_amp_long.to_hdf(f'{root}/dFF_analyzed.h5', key='amp', mode='w', format='table')
res_traces_avg.to_hdf(f'{root}/dFF_analyzed.h5', key='traces', format='table')
ROI_metadata.to_hdf(f'{root}/dFF_analyzed.h5', key='roi_metadata', format='table')
res_slope.to_hdf(f'{root}/dFF_analyzed.h5', key='slope', format='table')


# %% ---------------------------------------------------
traces_avg = res_traces_avg.assign(
    ROI_id = res_traces_avg['fish_id'] + '_' + res_traces_avg['ROI'].astype(str),
    which_exp = 'nMLF',
    stimulus = res_traces_avg['nsti'].map(sti_map),
)
traces_avg.loc[traces_avg['fish_info'].str.contains('S'), 'which_exp'] = 'TAN'

slope = res_slope.assign(
    ROI_id = res_slope['fish_id'] + '_' + res_slope['ROI'].astype(str),
    which_exp = 'nMLF',
)
slope.loc[slope['fish_info'].str.contains('S'), 'which_exp'] = 'TAN'

amp_long = res_amp_long.assign(
    ROI_id = res_amp_long['fish_id'] + '_' + res_amp_long['ROI'].astype(str),
    NS = 'nMLF',
    stimulus = res_amp_long['nsti'].map(sti_map),
)
amp_long.loc[amp_long['fish_info'].str.contains('S'), 'which_exp'] = 'TAN'

#%% ---------------- FUNCTION STOPS HERE -----------------
#%% get slope QC

all_roi = slope.ROI_id.unique()
# QC

max_serial_diff_amp_threshold = 0.3
# max_baseline_amp = 5
peakAmp_threshold = 0.5
endRes_threshold = 0.7
r_threshold = 0.5
# sq_residual_avg_threshold = 100
slope_threshold = 0.08


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


# %% plot
STIMULUS_EXT = [0,5,10,20,30]
amp_QC = amp_QC.assign(
    exp_cond_ordered = amp_QC['cond_num'].astype(str) + amp_QC['exp_cond']
)
df = amp_goodFit

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


# %%
