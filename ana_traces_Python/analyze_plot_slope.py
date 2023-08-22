#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from plot_functions.plt_functions import plt_categorical_grid2
from doQC_getSlope import doQC_getSlope
import matplotlib.pyplot as plt


#%%

root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT"


#%%

amp_QC, amp_goodFit, amp_goodTuning = doQC_getSlope(root)
fig_dir = f"{root}/figures"
try:
    os.makedirs(fig_dir)
except:
    pass
# %% plot

toplt = amp_goodFit.loc[amp_goodFit['cond_num'].isin([1,2])]

x_name='exp_cond_ordered'
y_name='slopeLowAng'
units = 'ROI_id'

p = plt_categorical_grid2(
    gridrow='which_exp',
    data=toplt,
    x_name=x_name,
    y_name=y_name,
    units=units,
    height=3,
    aspect=1,
)
plt.savefig(f"{fig_dir}/goodFit {y_name}X{x_name}.pdf", format='PDF')
# %%
toplt = amp_goodTuning.loc[amp_goodTuning['cond_num'].isin([1,2])]

x_name='exp_cond_ordered'
y_name='slopeLowAng'
units = 'ROI_id'

p = plt_categorical_grid2(
    gridrow='which_exp',
    # gridcol='peak_cat_cond1',
    data=toplt,
    x_name=x_name,
    y_name=y_name,
    units=units,
    height=3,
    aspect=1,
)
plt.savefig(f"{fig_dir}/goodTuning {y_name}X{x_name}.pdf", format='PDF')

# %%
# calculate percentage change
cond1_tuning = amp_goodTuning.query("cond_num == 1")
cond2_tuning = amp_goodTuning.query("cond_num == 2")
cond1_tuning = cond1_tuning.assign(
    slope_chg_ratio = (cond2_tuning['slopeLowAng'].values - cond1_tuning['slopeLowAng'].values) /  (cond2_tuning['slopeLowAng'].values + cond1_tuning['slopeLowAng'].values)
)
y_name='slope_chg_ratio'
x_name='which_exp'

plt_categorical_grid2(
    data=cond1_tuning,
    y_name=y_name,
    x_name=x_name,
    units='ROI_id'
)
plt.savefig(f"{fig_dir}/normSlopeChg {y_name}X{x_name}.pdf", format='PDF')

# %%
