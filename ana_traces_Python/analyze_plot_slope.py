#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.plt_functions import plt_categorical_grid2
from functions.doQC_getSlope import doQC_getSlope, sel_exp
import matplotlib.pyplot as plt


#%%
sel_dir = 'light' # lesion or light

#%%

if 'lesion' in sel_dir:
    which_slope = 'slopeAll_rawAmp'
else:
    which_slope = 'slopeLowAng'

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
_, slope = doQC_getSlope(root, cond4qc)

# %% plot
sel_qc = 'good_fit'
toplt = slope.loc[slope[sel_qc]]
toplt = toplt.loc[toplt['cond_num'].isin([1,2])]

x_name='exp_cond_ordered'
y_name=which_slope
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
sel_qc = 'good_tuning'
toplt = slope.loc[slope[sel_qc]]
toplt = toplt.loc[toplt['cond_num'].isin([1,2])]

x_name='exp_cond_ordered'
y_name=which_slope
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
# calculate  change
cond1_tuning = toplt.query("cond_num == 1")
cond2_tuning = toplt.query("cond_num == 2")
cond1_tuning = cond1_tuning.assign(
    slope_chg = cond2_tuning[which_slope].values - cond1_tuning[which_slope].values,
    slope_chg_ratio = (cond2_tuning[which_slope].values - cond1_tuning[which_slope].values) /  cond1_tuning[which_slope].values,
    slope_chg_norm = (cond2_tuning[which_slope].values - cond1_tuning[which_slope].values) /  (cond2_tuning[which_slope].values + cond1_tuning[which_slope].values)
)

exclude_for_plotting1 = cond1_tuning.loc[cond1_tuning['slope_chg_norm'].abs() > 1].ROI_id.unique()
exclude_for_plotting2 = cond1_tuning.loc[(cond1_tuning['amp_1'] < 0)].ROI_id.unique()
exclude_for_plotting = np.union1d(exclude_for_plotting1, exclude_for_plotting2)

df_change = cond1_tuning.loc[~cond1_tuning['ROI_id'].isin(exclude_for_plotting)]

for y_name in ['slope_chg', 'slope_chg_ratio', 'slope_chg_norm']:
    x_name='which_exp'

    g = plt_categorical_grid2(
        data=df_change,
        y_name=y_name,
        x_name=x_name,
        units='ROI_id'
    )
    if y_name != 'slope_chg_norm':
        g.set(ylim=[np.percentile(df_change[y_name], 1),np.percentile(df_change[y_name], 98.5)])
    else:
        g.set(ylim=[-1,1])
    plt.savefig(f"{fig_dir}/normSlopeChg {y_name}X{x_name}.pdf", format='PDF')

# %% ------- Peak time -------


# %%
x_name='exp_cond_ordered'
y_name=which_slope
units = 'ROI_id'

p = plt_categorical_grid2(
    gridrow='which_exp',
    # gridcol='peak_cat_cond1',
    data=toplt,
    x_name=x_name,
    y_name=y_name,
    gridcol='peak_cat',
    units=units,
    height=3,
    aspect=1,
)
plt.savefig(f"{fig_dir}/goodTuning_{y_name}X{x_name}_timing.pdf", format='PDF')

# %%
