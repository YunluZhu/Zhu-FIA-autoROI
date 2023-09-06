#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.plt_functions import plt_categorical_grid2
from functions.doQC_getSlope import doQC_getSlope_wTimedAvgAmp, sel_exp
from functions.doCategorization import get_peakTimingCat
import matplotlib.pyplot as plt

#%%
sel_dir = 'lesion' # lesion or light
sel_qc = 'good_tuning'

#%%
fig_root = f"/Users/yunluzhu/Documents/Lab2/caiman/Volumetric_code/YZ_nMLF_speed/figures"
fig_folder_name = sel_dir
STIMULUS_EXT = [0,5,10,20,30]
fig_dir = os.path.join(fig_root, fig_folder_name)
try:
    os.makedirs(fig_dir)
except:
    pass

root, cond4qc = sel_exp(sel_dir)

# %%
_, amp_smval, slope = doQC_getSlope_wTimedAvgAmp(root, cond4qc)
slope = slope.reset_index(drop=True)

traces = pd.read_hdf(f'{root}/res_concatenated.h5', key='long_data')

#%%
ROI_sel = amp_smval.loc[amp_smval['good_tuning'], 'ROI_id'].unique()
traces_sel = traces.query("exp_cond == 'before'")
traces_sel_toplt = traces_sel.query("ROI_id in @ROI_sel")
# %%
traces_sel_toplt = traces_sel_toplt.assign(
    adj_frames = traces_sel_toplt['frames'] + traces_sel_toplt['nsti'] * 12.5
)

#%%
plt.figure(figsize=[5,3])
g = sns.lineplot(
    traces_sel_toplt,
    x='adj_frames',
    y='dFF',
    hue='nsti',
    palette='dark',
    legend=False,
    errorbar=None
)
g = sns.lineplot(
    traces_sel_toplt,
    x='adj_frames',
    y='dFF',
    units='ROI_id',
    estimator=None,
    alpha=0.01,
    hue='nsti',
    legend=False,

)
g.set(
    ylim=[-0.5, 11]
)
plt.savefig(f"{fig_dir}/traces_good_tuning.pdf", format='PDF')
# %%
