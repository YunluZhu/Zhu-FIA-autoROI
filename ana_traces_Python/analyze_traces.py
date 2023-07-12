'''
nMLF pipeline visualization
prerequisite: fluorescent traces extracted by MATLAB script. 
HARD CODED for DURATION of STIMULUS 20s
for details of the whole pipeline, see: https://benchling.com/s/etr-GCcRV7cmhRg7u4JQuDfr?m=slm-n9tOEQ0aekst3JyASOkw

'''

# %%
import pandas as pd
import numpy as np
import seaborn as sns
import os
import scipy.io
import matplotlib.pyplot as plt
import itertools
from plot_functions.plt_functions import plt_categorical_grid


DATA_FILE = 'dFF_adj' # dFF_adj 
if_plot = False
sel_area = [1,2]
STIMULUS = [5, 10, 20, 30]
nsti = len(STIMULUS)

root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/_2analyze/230622_fish4 prox_bi"

# %%
raw_df = pd.read_csv(os.path.join(root, "rawF_df.csv"))
    
mat = scipy.io.loadmat(os.path.join(root, DATA_FILE+'.mat'))[DATA_FILE]

dFF = pd.DataFrame(data=mat)
dFF = dFF.assign(
    area = raw_df.area,
    repeat = raw_df.repeat,
)
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

parent_dir = os.path.abspath(os.path.join(root, os.pardir))
for file in os.listdir(parent_dir):
    if str.endswith(file, " metadata.csv"):
        fish_metaDATA_FILE = os.path.join(parent_dir, file)
        break
    
fish_metadata = pd.read_csv(fish_metaDATA_FILE)

frame_rate = fish_metadata.loc[
    (fish_metadata['exp_date']==exp_date) & (fish_metadata['fish_num']==fishNum), 'frame_rate'
    ]
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

if if_plot:

    sns.relplot(
        data = dFF_long_roi_corrected.loc[dFF_long_roi_corrected['area'].isin(sel_area)],
        kind = 'line',
        row = 'ROI',
        col = 'area',
        x = 'time',
        y = DATA_FILE,
        hue='area',
        facet_kws={'sharey': "row", 'sharex': True},
        height=2.5,
        errorbar=None,
    )
    filename = os.path.join(root, f"{DATA_FILE} {''.join(map(str, sel_area))}_roiConsolidated.pdf")
    plt.savefig(filename,format='PDF')

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

dFF_THRESHOLD = 0.5  # mean pean amp of dFF
BASE_THRESHOLD = 0.8  # std/mean

# 1. find ROIs with mean peak amp across all stimulus greater than threshold 

ROI_pass1 = amp_long[amp_long.groupby(sel_unique_roi_level)[DATA_FILE].transform('mean') > dFF_THRESHOLD].query("area == 1").ROI.unique()

# 2. find ROIs with relatively stable baseline after each sti
if str.endswith(DATA_FILE, 'adj'):
    baseline_file = 'baseline_adj'
else:
    baseline_file = 'baseline'
    
mat = scipy.io.loadmat(os.path.join(root, baseline_file+'.mat'))
baseline_file_name = list(mat.keys())[-1]

baseline_wide = pd.DataFrame(data=mat[baseline_file_name])
idx_for_head1EachSti = np.arange(0, len(baseline_wide), trial_frames)
idx_for_head1EachSti = idx_for_head1EachSti[0:len(dFF.query("area == 1").repeat.unique())] # only use baseline of first condition
bb = baseline_wide.loc[idx_for_head1EachSti]
bb = bb.std()/bb.mean()
ROI_pass2 = bb.loc[bb < BASE_THRESHOLD].index

ROI_passQC = np.intersect1d(ROI_pass1, ROI_pass2)

amp_selected = amp_long.loc[amp_long['ROI'].isin(ROI_passQC)]

# %%

amp_wide = pd.pivot(
    amp_long,
    index=sel_unique_roi_level,
    columns='nsti',
    values=DATA_FILE,
).reset_index()

amp_wide.columns.name = None
amp_wide.columns = sel_unique_roi_level + amp_cols
amp_wide = amp_wide.rename(columns={'area':'cond_num'})
amp_wide = amp_wide.merge(ROI_metadata, how='left', left_on='ROI', right_on='id')

amp_wide = amp_wide.assign(
    exp_cond = amp_wide['cond_num'].astype(str).map(area_cond_dict),
)
# %%
slope = []
sqErr = []

X = np.array(STIMULUS)
X = np.vstack((np.ones(len(X)), X)).T

for index, row in amp_wide.iterrows():
    Y = np.float64(row[amp_cols].values)
    b = (np.linalg.inv(X.T @ X) @ (X.T) @ Y.T)[1]
    slope.append(b)
    e = np.sum(np.square(Y - X.dot(b)[:,1]))
    sqErr.append(e)

amp_wide = amp_wide.assign(
    amp_slope = slope,
    sqErr = sqErr
)

# %%
# use squared error to find ROIs with tuning

tuning_qc_df = amp_wide.groupby(['cond_num', 'ROI'])[['amp_slope', 'sqErr']].median()
sqErr_threshold = np.percentile(tuning_qc_df.loc[1,:,:]['sqErr'], 95) 
ROI_goodFit = tuning_qc_df.loc[1,:,:].query("sqErr < @sqErr_threshold").index

amp_goodFit = amp_wide.loc[amp_wide['ROI'].isin(ROI_goodFit)]
# %%
x_name='exp_cond'
y_name='amp_slope'
units = 'ROI'

slope_avg = amp_goodFit.groupby([units, x_name])[y_name].median().reset_index()


p = plt_categorical_grid(
    gridcol=None,
    gridrow=None,
    data=slope_avg,
    x_name=x_name,
    y_name=y_name,
    units=units,
    height=3
)



# %%

dFF_long_roi_corrected.to_hdf(f'{root}/dFF_analyzed.h5', key='long_data', mode='w', format='table')
ROI_metadata.to_hdf(f'{root}/dFF_analyzed.h5', key='roi_metadata', format='table')


# %%



# %%
