'''
prerequisite: fluorescent traces extracted by MATLAB script. 
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
import plotly.express as px
import plotly.graph_objs as go

DATA_FILE = 'dFF_ksDensity' # dFF_adj dFF_ksDensity rawF
root = "/Volumes/LabDataPro/CB/111111_fish1 blahblah"
if_share_Y_axis_scale = False

# %%
fig_dir = os.path.join(os.getcwd(),'figures')
try:
    os.makedirs(fig_dir)
except:
    pass

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
time_stamp = np.arange(0,trial_frames/vol_rate+1/vol_rate,1/vol_rate)[1:]

# sti_frames = trial_frames / nsti

#%%
dFF_long = pd.wide_to_long(dFF.reset_index(), stubnames='roi', i='index', j='ROI', sep='_').reset_index()
dFF_long = dFF_long.rename(columns={'roi':DATA_FILE})
dFF_long = dFF_long.assign(
    exp_cond = dFF_long['area'].astype(str).map(area_cond_dict),
    fish_id = fish_info[0],
    fish_info = fish_info[1],
    time_stamp = dFF_long.frames.map({dFF_long.frames.unique()[i]: time_stamp[i] for i in range(len(time_stamp))})
)

ROI_metadata = ROI_metadata.assign(
    fish_id = fish_info[0],
    fish_info = fish_info[1],
)

# %%
f1 = sns.relplot(
    data=dFF_long,
    x='time_stamp',
    y=DATA_FILE,
    row='ROI',
    height=2,
    aspect=2.5,
    kind='line',
    facet_kws={'sharey': if_share_Y_axis_scale, 'sharex': True}
)

if if_share_Y_axis_scale:
    fig_suffix = 'shareY'
else:
    fig_suffix = ""
    
plt.savefig(fig_dir+f"/ROI_traces_{fig_suffix}.pdf",format='PDF')
# %%
