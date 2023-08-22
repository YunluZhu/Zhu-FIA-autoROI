'''
nMLF pipeline visualization
prerequisite: fluorescent traces extracted by MATLAB script. 
HARD CODED for DURATION of STIMULUS 20s
for details of the whole pipeline, see: https://benchling.com/s/etr-GCcRV7cmhRg7u4JQuDfr?m=slm-n9tOEQ0aekst3JyASOkw

NOTE reads dF_ksDensity and uses Dark KS baseline (which includes DarkBase) to calculate dFF for all conditions
NOTE slopeAll on raw data with repeats
-----todo-----
NOTE under light, response to lower angles are reduced. Compare amplituds
'''

# %%
import pandas as pd
import numpy as np
import seaborn as sns
import os
import scipy.io
import matplotlib.pyplot as plt
import itertools
from plot_functions.plt_functions import plt_categorical_grid2
import plotly.express as px
# import plotly.graph_objs as go


# if_plot = False

STIMULUS = [5, 10, 20, 30]
nsti = len(STIMULUS)

# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT/230727_fish1 NT"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT/230808_fish1 TNY"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT/230808_fish2 TNY"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT/230809_fish1 TNY"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT/230810_fish3 NTY"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT/230810_fish4 NTY"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT/230814_fish2 SY"
# root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT/new/230814_fish4 SY"

root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT/230808_fish1 TNY"

# %%
fig_dir = f"{root}/figures"
try:
    os.makedirs(fig_dir)
except:
    print('Notes: re-writing old figures')
    
DATA_FILE = 'dF_ksDensity' # dFF_adj dFF_ksDensity rawF
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
# initial filtering
# do responses increase as trials ++? if so, fish has not been properly waked up before imaging started


dFF_THRESHOLD = 0.5  # mean peak amp of dFF
# dFF_stdMean_THRESHOLD = 10
# BASE_THRESHOLD = 0.8  # std/mean
LAST3s_THRESHOLD = 0.5 # dFF



# 1. find ROIs with mean peak amp across all stimulus greater than threshold 
amp_long_cond1_group = amp_long.query("area == 1").groupby(sel_unique_roi_level)[DATA_FILE]
ROI_pass1 = amp_long.query("area == 1")[(amp_long_cond1_group.transform('median') > dFF_THRESHOLD)].ROI.unique()
print(f"{len(ROI_pass1)}/{len(ROI_metadata)} ROIs passed dFF threshold")

# 2. find ROIs with relatively stable responses across repeats
amp_stability = amp_long.query("area == 1").groupby(['ROI','nsti'])[DATA_FILE].apply(
    lambda x: x.std()/x.mean()
).reset_index()
amp_stability = amp_stability.groupby('ROI').mean().reset_index()
amp_stability_th = np.percentile(amp_stability[DATA_FILE],95)
ROI_pass2 = amp_stability.loc[amp_stability[DATA_FILE]<amp_stability_th].ROI.unique()
print(f"{len(ROI_pass2)}/{len(ROI_metadata)} ROIs passed dFF amp stability threshold")


# 3. find ROIs with response down to baseline at the end of each sti
idx_for_last3sec = dFF_long_roi_corrected.groupby(sel_unique_roi_level + ['nsti']).tail(int(np.floor(3 * vol_rate))).index
rows_for_last3sec = dFF_long_roi_corrected.loc[idx_for_last3sec]
last3sec_long = rows_for_last3sec.groupby(sel_unique_roi_level + ['nsti'])[DATA_FILE].mean(numeric_only=True).reset_index()
last3sec_long = last3sec_long.query("area == 1")
last3sec_long_groupMean = last3sec_long.groupby(['ROI','area'])[DATA_FILE].mean()
ROI_pass3 = last3sec_long_groupMean.loc[last3sec_long_groupMean < LAST3s_THRESHOLD].reset_index().ROI.unique()
print(f"{len(ROI_pass3)}/{len(ROI_metadata)} ROIs passed baseline threshold")


ROI_pass12 = np.intersect1d(ROI_pass1, ROI_pass2)
ROI_passQC = np.intersect1d(ROI_pass12, ROI_pass3)

print(f"> {len(ROI_passQC)} passed QC")
amp_selected = amp_long.loc[amp_long['ROI'].isin(ROI_passQC)]
dFF_long_roi_selected = dFF_long_roi_corrected.loc[dFF_long_roi_corrected['ROI'].isin(ROI_passQC)]

#%%
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
    peak_cat = pd.cut(df_peak_time_cat.peak_time, bins=[-np.inf, 1.5, 3, 5.5, np.inf], labels=['fast', 'mid', 'slow', 'late'])
)

df_peak_time_cat_area1 = df_peak_time_cat.query("area==1")

dFF_long_roi_averaged = dFF_long_roi_averaged.merge(df_peak_time_cat_area1[['ROI','fish_id','peak_cat','peak_time']], on=['fish_id', 'ROI'])

# %%
# change of response to individual stimuli

############### % check raw traces ##############
df_toplt = dFF_long_roi_averaged.loc[dFF_long_roi_averaged['peak_cat'].isin(['fast','mid'])]
df_toplt.peak_cat = df_toplt.peak_cat.astype(str)
g = sns.relplot(row='exp_cond', col='peak_cat', kind='line', data=df_toplt, x='frames',y=DATA_FILE,units='ROI', estimator=None, alpha=0.1, aspect=3, height=1.8,
                col_order=['fast','mid'])
g.set(ylim=[-1,12])
plt.savefig(f"{fig_dir}/traces.pdf", format='PDF')

# df_toplt = dFF_long_roi_selected.loc[dFF_long_roi_selected['ROI'] == ROI_passQC[k]]
# df_toplt = df_toplt.groupby(['frames','exp_cond','fish_id','nsti','area','ROI',])[DATA_FILE].median().reset_index()
# g = sns.relplot(row='exp_cond', kind='line', data=df_toplt, x='frames',y=DATA_FILE, alpha=.5, aspect=3, height=1.5,
#                 # row_order=['dark','light'],
#                 )
# g.set(ylim=[0,10])

# k += 1


############## RAW trace analysis done ##################

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

last3sec_avg = rows_for_last3sec.loc[rows_for_last3sec['ROI'].isin(ROI_passQC)].groupby(['ROI','repeat','area'])[DATA_FILE].median()
last3sec_avg.name = 'last3sec' 
last3sec_avg = last3sec_avg.reset_index()
 
amp_wide = amp_wide.merge(last3sec_avg, left_on=['cond_num', 'repeat', 'ROI'], right_on=['area', 'repeat', 'ROI'])


# %%
################## slopeAll on raw amplitudes #################################

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
    r = np.corrcoef(g_X[:,1], g_Y)
    corr.append(r[0,1])
    
    g_XlowAng = np.vstack((np.ones(len(XlowAng)), XlowAng)).T
    g_XlowAng = np.tile(g_XlowAng, (g_num_of_repeat, 1))
    g_YlowAng = np.float64(group[amp_colslowAng]).flatten()
    b = (np.linalg.inv(g_XlowAng.T @ g_XlowAng) @ (g_XlowAng.T) @ g_YlowAng.T)
    slopeLowAng.append(b[1])
    e = np.sum(np.square(g_YlowAng - g_XlowAng.dot(b)))
    sq_residualLowAng.append(e)
    sq_residualLowAng_avg.append(e/g_num_of_repeat)
    r = np.corrcoef(g_XlowAng[:,1], g_YlowAng)
    corrLowAng.append(r[0,1])

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
)

tocalc_slopeAllDiff = amp_avg.loc[amp_avg['cond_num'].isin([1, 2])]
tocalc_slopeAllDiff = tocalc_slopeAllDiff.sort_values(by=['ROI'])
slopeAll_cond2 = tocalc_slopeAllDiff.query("cond_num == 2")
slopeAll_cond1 = tocalc_slopeAllDiff.query("cond_num == 1")
slopeAll_chg = slopeAll_cond1[['cond_num', 'ROI', 'exp_cond']].copy()
slopeAll_chg = slopeAll_chg.assign(
    slopeAll_chg = slopeAll_cond2['slopeAll_rawAmp'].values - slopeAll_cond1['slopeAll_rawAmp'].values,
    slopeAll_chg_ratio = (slopeAll_cond2['slopeAll_rawAmp'].values - slopeAll_cond1['slopeAll_rawAmp'].values)/slopeAll_cond1['slopeAll_rawAmp'].values
)
amp_avg = amp_avg.merge(df_peak_time_cat[['fish_id', 'ROI', 'peak_time', 'area', 'peak_cat']], left_on=['ROI','cond_num'], right_on=['ROI','area'])

# %%
# # use squared error to find ROIs with tuning 

r_threshold = 0.6
# sq_residual_avg_threshold = 5 ##### to be changed

sel_fitted_forQC = amp_avg.loc[amp_avg['cond_num']==1,:]
ROI_goodFit1 = sel_fitted_forQC.query("r > @r_threshold").ROI.values
# ROI_goodFit2 = sel_fitted_forQC.query("sq_residual_avg < @sq_residual_avg_threshold").ROI.values
ROI_goodFit = ROI_goodFit1 #np.intersect1d(ROI_goodFit1, ROI_goodFit2)

amp_goodFit_slopeAllOnAveragedAmp = amp_avg.loc[amp_avg['ROI'].isin(ROI_goodFit)]

amp_goodFit_slopeAllOnAveragedAmp = amp_goodFit_slopeAllOnAveragedAmp.assign(
    exp_cond_ordered = amp_goodFit_slopeAllOnAveragedAmp['cond_num'].astype(str) + amp_goodFit_slopeAllOnAveragedAmp['exp_cond']
)

print(f'Ratio of ROIs pass regression QC: {len(ROI_goodFit)/len(sel_fitted_forQC.ROI.unique())}')

# %% plot

x_name='exp_cond_ordered'
y_name='slopeLowAng'
units = 'ROI'

p = plt_categorical_grid2(
    gridcol=None,
    gridrow=None,
    data=amp_goodFit_slopeAllOnAveragedAmp,
    x_name=x_name,
    y_name=y_name,
    units=units,
    height=3,
    aspect=1.3,
)








# %%
# %% plot good tuning ROI

slopeAll_threshold = np.percentile(amp_goodFit_slopeAllOnAveragedAmp.slopeAll_rawAmp, 80)
print(f"slopeAll thresholded {slopeAll_threshold}")
ROI_goodTuning = amp_goodFit_slopeAllOnAveragedAmp.query("cond_num == 1 & slopeAll_rawAmp > @slopeAll_threshold").ROI.values
amp_goodTuning = amp_goodFit_slopeAllOnAveragedAmp.loc[amp_goodFit_slopeAllOnAveragedAmp.ROI.isin(ROI_goodTuning)]

sns.lineplot(data=amp_goodTuning, 
             x='cond_num', 
             y='slopeAll_rawAmp')

x_name='exp_cond_ordered'
y_name='slopeAll_rawAmp'
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

# selROI = amp_goodFit_slopeAllOnAveragedAmp.sample(1).ROI.values[0]
# ppdf = amp_goodFit_slopeAllOnAveragedAmp.query("ROI == @selROI")
# sns.lineplot(data=ppdf, x='cond_num', y='slopeAll_rawAmp')
#%%























# %% look at raw data
sti_map = dict([(ii+1, sti) for ii, sti  in enumerate(STIMULUS)])

tuning_goodFit_longdf = amp_selected.loc[amp_selected.ROI.isin(amp_goodFit_slopeAllOnAveragedAmp.ROI.unique())]
tuning_goodFit_longdf = tuning_goodFit_longdf.assign(
    stimulus = tuning_goodFit_longdf['nsti'].map(sti_map),
    exp_cond = tuning_goodFit_longdf['area'].astype(str).map(area_cond_dict),

)

average_tuning_goodFit = tuning_goodFit_longdf.groupby(['exp_cond','ROI','stimulus']).median().reset_index()
sns.relplot(
    kind='line',
    data=average_tuning_goodFit,
    row='exp_cond',
    x='stimulus',
    y=DATA_FILE,
    estimator=None,
    units='ROI',
    alpha=0.2,
    palette=sns.color_palette("tab10"),
)

# %%
ROI_metadata = ROI_metadata.assign(
    if_goodFit=0,
    if_goodTuning=0, 
    slopeAll_area1=0,
)
ROI_metadata.loc[ROI_metadata['id'].isin(ROI_goodFit), 'if_goodFit'] = 1
ROI_metadata.loc[ROI_metadata['id'].isin(ROI_goodTuning), 'if_goodTuning'] = 1
ROI_metadata = ROI_metadata.merge(amp_avg.query("cond_num == 1")[['ROI', 'amp_4', 'slopeAll_avgAmp']], how='left',left_on='id', right_on='ROI')

slopeAll_chg = slopeAll_chg.merge(ROI_metadata, how='left', left_on='ROI', right_on='id')

# %%
# look at diff of slopeAlls between conditions

df = slopeAll_chg.query("if_goodFit == 1")

plt.figure()
sns.scatterplot(data=df,
                x='slopeAll_chg',
                y='slopeAll_chg_ratio')
plt.figure()
sns.histplot(data=df,
                y='slopeAll_chg',
                bins=4)
# %%
# look at location of good tuning neurons


df = slopeAll_chg.query("if_goodTuning == 1")
fig = px.scatter_3d(df, x='xCenter', y='yCenter', z='zPos',
              color='slopeAll_chg_ratio',
              color_continuous_scale=px.colors.diverging.balance
              )
fig.update_layout(
    scene=dict(
           aspectmode='data', #this string can be 'data', 'cube', 'auto', 'manual'
           )
)
fig.show()
# %%

dFF_long_roi_corrected.to_hdf(f'{root}/dFF_analyzed.h5', key='long_data', mode='w', format='table')
ROI_metadata.to_hdf(f'{root}/dFF_analyzed.h5', key='roi_metadata', format='table')



# %%
# 