'''
nMLF pipeline visualization
prerequisite: fluorescent traces extracted by MATLAB script. 
HARD CODED for DURATION of STIMULUS 20s
for details of the whole pipeline, see: https://benchling.com/s/etr-GCcRV7cmhRg7u4JQuDfr?m=slm-n9tOEQ0aekst3JyASOkw

for dF/dFF calculated using KSDensity method or raw F
reads .mat files

'''

# %%
import pandas as pd
import numpy as np
import seaborn as sns
import os
import scipy.io
import matplotlib.pyplot as plt
import itertools
from functions.plt_functions import plt_categorical_grid2
import plotly.express as px
import plotly.graph_objs as go


DATA_FILE = 'dF_ksDensity' # dFF_adj dFF_ksDensity rawF
if_plot = False
# sel_area = [1,2]
STIMULUS = [5, 10, 20, 30]
nsti = len(STIMULUS)

# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/lightanalyzed_UP/230630_fish1 ST"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/lightanalyzed_UP/230714_fish2 NT"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze/230618_fish1 prox_R"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/analyzed_UP/230620_fish1 prox_bi"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/analyzed_UP/230622_fish2 distal_bi"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/analyzed_UP/230624_fish3 distal_bi"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/analyzed_UP/230626_fish1 distal_bi"

# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/analyzed_UP/230714_fish1 prox_bi"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/analyzed_DN/230725_fish1 NOSE DOWN"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/analyzed_UP/230727_fish2 distal_bi"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT/230727_fish1 NT"

root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT/ana/230810_fish4 NTY"
# %%
raw_df = pd.read_csv(os.path.join(root, "rawF_df.csv"))
mat = scipy.io.loadmat(os.path.join(root, DATA_FILE+'.mat'))[DATA_FILE]
dFF = pd.DataFrame(data=mat)

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

dFF_long_roi_corrected = dFF_long.copy()

n_combined = 1000
all_combined = 0

while n_combined > 0:
    unique_roi = dFF_long_roi_corrected['ROI'].unique()
    similar_roi = []
    n_combined = 0
    for id1, id2 in itertools.combinations(unique_roi, 2):
        roi1 = ROI_metadata.query("id == @id1").iloc[0]
        roi2 = ROI_metadata.query("id == @id2").iloc[0]
        if_adjacent_z = abs(roi1.zPos - roi2.zPos) == 1
        if_close_x = abs(roi1.xCenter - roi2.xCenter) <= 2
        if_close_y = abs(roi1.yCenter - roi2.yCenter) <= 2
        if_spacially_same = if_adjacent_z & if_close_x & if_close_y
        if if_spacially_same:
            similar_roi.append((id1, id2))
            
    r_out = []   
    # response correlation
    for id1, id2 in similar_roi:
        res1 = dFF_long_roi_corrected.query("ROI == @id1")
        res2 = dFF_long_roi_corrected.query("ROI == @id2")
        this_corr = np.corrcoef(res1[DATA_FILE].values, res2[DATA_FILE].values)
        this_r = this_corr[0,1]
        r_out.append(this_r)

    for i, this_r in enumerate(r_out):
        id1, id2 = similar_roi[i]
        if this_r > 0.85:  # if two nearby ROIs have similar response, combine them
            unique_roi = dFF_long_roi_corrected['ROI'].unique()
            if (id1 in unique_roi) & (id2 in unique_roi):
                n_combined += 1
                res1 = dFF_long_roi_corrected.query("ROI == @id1")
                res2 = dFF_long_roi_corrected.query("ROI == @id2")
                res_combined = res1.copy()
                res_combined[DATA_FILE] = (res1[DATA_FILE].values + res2[DATA_FILE].values)/2
                dFF_long_roi_corrected.loc[dFF_long_roi_corrected['ROI']==id1, DATA_FILE] = res_combined[DATA_FILE]
                dFF_long_roi_corrected = dFF_long_roi_corrected.drop(index=res2.index)
    
    all_combined += n_combined
    
print(f'{all_combined} ROIs combined')

# %%
if if_plot:
    df = df = dFF_long_roi_corrected.groupby(['frames','exp_cond','fish_id','area','ROI'])[DATA_FILE].mean().reset_index()
    g = sns.relplot(row='exp_cond', kind='line', data=df,x='frames',y=DATA_FILE,units='ROI', estimator=None, alpha=0.1, aspect=3)
    g.set(ylim=[0,20])



# %%
# Extract peaks and make a new dataframe of peaks
# the dFF_long_roi_corrected should contain time column

# find the idx of time of response during a certain winder after each sti
amp_cols = ['amp_' + str(sti) for sti in np.arange(nsti)+1]
time_for_amp = 2 #sec
frame_for_amp = int(np.floor(time_for_amp * vol_rate))

sel_unique_roi_level = ['area','repeat','ROI']
idx_for_amp = dFF_long_roi_corrected.groupby(sel_unique_roi_level + ['nsti']).head(frame_for_amp).index
rows_for_amp = dFF_long_roi_corrected.loc[idx_for_amp]
amp_long = rows_for_amp.groupby(sel_unique_roi_level + ['nsti']).mean(numeric_only=True).reset_index()

# %%
# QC

dFF_THRESHOLD = 0.5  # mean peak amp of dFF
dFF_stdMean_THRESHOLD = 1
BASE_THRESHOLD = 0.8  # std/mean
LAST3s_THRESHOLD = 0.75 # dFF
# 1. find ROIs with mean peak amp across all stimulus greater than threshold 
amp_long_group = amp_long.groupby(sel_unique_roi_level)[DATA_FILE]
ROI_pass1 = amp_long[(amp_long_group.transform('mean') > dFF_THRESHOLD) &  (amp_long_group.transform('std')/amp_long_group.transform('mean') < dFF_stdMean_THRESHOLD)].query("area == 1").ROI.unique()
print(f"{len(ROI_pass1)}/{len(ROI_metadata)} ROIs passed dFF threshold")

# # 2. find ROIs with relatively stable baseline after each sti
# baseline_file = 'baseline' + '_' + DATA_FILE.split('_')[1]
# if DATA_FILE.endswith('adj'):
#     baseline_file = baseline_file + '_adj'
# mat = scipy.io.loadmat(os.path.join(root, baseline_file+'.mat'))
# baseline_file_name = list(mat.keys())[-1]

# baseline_wide = pd.DataFrame(data=mat[baseline_file_name])
# bb = baseline_wide.std()/baseline_wide.mean()
# ROI_pass2 = bb.loc[bb < BASE_THRESHOLD].index + 1
# print(f"{len(ROI_pass2)}/{len(ROI_metadata)} ROIs passed stable baseline threshold")

# 3. find ROIs with response down to baseline at the end of each sti
idx_for_last3sec = dFF_long_roi_corrected.groupby(sel_unique_roi_level + ['nsti']).tail(int(np.floor(3 * vol_rate))).index
rows_for_last3sec = dFF_long_roi_corrected.loc[idx_for_last3sec]
last3sec_long = rows_for_last3sec.groupby(sel_unique_roi_level + ['nsti'])[DATA_FILE].median(numeric_only=True).reset_index()
last3sec_long_groupMean = last3sec_long.groupby(['ROI','area'])[DATA_FILE].median()
ROI_pass3 = last3sec_long_groupMean.loc[last3sec_long_groupMean < LAST3s_THRESHOLD].reset_index().ROI.unique()
print(f"{len(ROI_pass3)}/{len(ROI_metadata)} ROIs passed baseline threshold")


ROI_passQC = np.intersect1d(ROI_pass1, ROI_pass3)
# ROI_passQC = np.intersect1d(ROI_passQC, ROI_pass3)

print(f"> {len(ROI_passQC)} passed QC")
amp_selected = amp_long.loc[amp_long['ROI'].isin(ROI_passQC)]

# %% 
################## SLOPE on average amplitudes #################################
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

slope = []
sq_residual = []
corr = []

amp_avg = amp_wide.groupby(['cond_num', 'ROI','exp_cond'])[amp_cols].median().reset_index()

X = np.array(STIMULUS)
X = np.vstack((np.ones(len(X)), X)).T

for index, row in amp_avg.iterrows():
    Y = np.float64(row[amp_cols].values)
    b = (np.linalg.inv(X.T @ X) @ (X.T) @ Y.T)
    slope.append(b[1])
    e = np.sum(np.square(Y - X.dot(b)))
    sq_residual.append(e)
    r = np.corrcoef(X[:,1], Y)
    corr.append(r[0,1])

amp_avg = amp_avg.assign(
    slope_avgAmp = slope,
    sq_residual_avgAmp = sq_residual,
    r = corr
)

# # use squared error to find ROIs with tuning 

r_threshold = 0.8
ratio_pass = sum(amp_avg.loc[amp_avg['cond_num']==1,'r'] > r_threshold)/len(amp_avg.ROI.unique())
print(f'Ratio of ROIs with > {r_threshold} correlation of amplitudes X stimulus: {ratio_pass}')
ROI_goodFit = amp_avg.loc[amp_avg['cond_num']==1,:].query("r > @r_threshold").ROI.values
amp_goodFit_SlopeOnAveragedAmp = amp_avg.loc[amp_avg['ROI'].isin(ROI_goodFit)]

amp_goodFit_SlopeOnAveragedAmp = amp_goodFit_SlopeOnAveragedAmp.assign(
    exp_cond_ordered = amp_goodFit_SlopeOnAveragedAmp['cond_num'].astype(str) + amp_goodFit_SlopeOnAveragedAmp['exp_cond']
)
x_name='exp_cond_ordered'
y_name='slope_avgAmp'
units = 'ROI'

p = plt_categorical_grid2(
    gridcol=None,
    gridrow=None,
    data=amp_goodFit_SlopeOnAveragedAmp,
    x_name=x_name,
    y_name=y_name,
    units=units,
    height=3,
    aspect=1.3,
)

# %% plot good tuning ROI
slope_threshold = np.percentile(amp_goodFit_SlopeOnAveragedAmp.slope_avgAmp, 75)
print(f"Slope thresholded {slope_threshold}")
ROI_goodTuning = amp_goodFit_SlopeOnAveragedAmp.query("cond_num == 1 & slope_avgAmp > @slope_threshold").ROI.values
amp_goodTuning = amp_goodFit_SlopeOnAveragedAmp.loc[amp_goodFit_SlopeOnAveragedAmp.ROI.isin(ROI_goodTuning)]

sns.lineplot(data=amp_goodTuning, 
             x='cond_num', 
             y='slope_avgAmp')

x_name='exp_cond_ordered'
y_name='slope_avgAmp'
units = 'ROI'

p = plt_categorical_grid2(
    gridcol=None,
    gridrow=None,
    data=amp_goodTuning,
    x_name=x_name,
    y_name=y_name,
    units=units,
    height=3,
    aspect=1.3
)


# # %% Check one ROI tuning

# selROI = amp_goodFit_SlopeOnAveragedAmp.sample(1).ROI.values[0]
# ppdf = amp_goodFit_SlopeOnAveragedAmp.query("ROI == @selROI")
# sns.lineplot(data=ppdf, x='cond_num', y='slope_avgAmp')

# %% look at raw data
sti_map = dict([(ii+1, sti) for ii, sti  in enumerate(STIMULUS)])

tuning_goodFit_longdf = amp_selected.loc[amp_selected.ROI.isin(amp_goodFit_SlopeOnAveragedAmp.ROI.unique())]
tuning_goodFit_longdf = tuning_goodFit_longdf.assign(
    stimulus = tuning_goodFit_longdf['nsti'].map(sti_map),
    exp_cond = tuning_goodFit_longdf['area'].astype(str).map(area_cond_dict),

)

average_tuning_goodFit = tuning_goodFit_longdf.groupby(['exp_cond','ROI','stimulus']).median().reset_index()
sns.lineplot(
    data=average_tuning_goodFit,
    x='stimulus',
    y=DATA_FILE,
    hue='exp_cond',
    palette=sns.color_palette("tab10"),
)


# %%

ROI_metadata = ROI_metadata.assign(
    if_goodFit=0,
    if_goodTuning=0, 
    slope_area1=0,
)
ROI_metadata.loc[ROI_metadata['id'].isin(ROI_goodFit), 'if_goodFit'] = 1
ROI_metadata.loc[ROI_metadata['id'].isin(ROI_goodTuning), 'if_goodTuning'] = 1
ROI_metadata = ROI_metadata.merge(amp_avg.query("cond_num == 1")[['ROI', 'amp_4', 'slope_avgAmp']], how='left',left_on='id', right_on='ROI')

dFF_long_roi_corrected.to_hdf(f'{root}/dFF_analyzed.h5', key='long_data', mode='w', format='table')
ROI_metadata.to_hdf(f'{root}/dFF_analyzed.h5', key='roi_metadata', format='table')

# %%
# look at diff of slopes between conditions

df_sel = amp_goodTuning
ctrl_df = df_sel.query("cond_num == 1").reset_index(drop=True)
cond_df = df_sel.query("cond_num == 2").reset_index(drop=True)
res = ctrl_df.copy()
res = res.assign(
    diff = cond_df['slope_avgAmp'] - ctrl_df['slope_avgAmp'], 
)
res = res.assign(
    ratio = res['diff'] / res['slope_avgAmp']
)
plt.figure()
sns.scatterplot(data=res,
                x='slope_avgAmp',
                y='diff')
plt.figure()
sns.histplot(data=res,
                y='diff',
                bins=4)
# %%
# look at location of good tuning neurons


df = ROI_metadata.query("if_goodFit == 1")
fig = px.scatter_3d(df, x='xCenter', y='yCenter', z='zPos',
              size='slope_avgAmp', color='amp_4')
fig.update_layout(
    scene=dict(
           aspectmode='data', #this string can be 'data', 'cube', 'auto', 'manual'
           )
)
fig.show()

# %%
