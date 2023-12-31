#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.plt_functions import plt_categorical_grid2
from functions.doQC_getSlope import doQC_getSlope_4lesion
import matplotlib.pyplot as plt


#%%

root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_lesion"
which_slope = 'slopeAll_rawAmp'

#%%

amp_QC, amp_goodFit, amp_goodTuning = doQC_getSlope_4lesion(root)
STIMULUS_EXT = [0,5,10,20,30]
fig_root = f"/Users/yunluzhu/Documents/Lab2/caiman/Volumetric_code/YZ_nMLF_speed/figures"
fig_folder_name = "lesionSlope"
fig_dir = os.path.join(fig_root, fig_folder_name)

try:
    os.makedirs(fig_dir)
except:
    pass

#%%
df_toplt = amp_goodFit
df_toplt.loc[df_toplt['cond_num']==2, 'exp_cond_ordered'] = '2lesion'

x_name='exp_cond_ordered'
y_name=which_slope
units = 'ROI_id'

p = plt_categorical_grid2(
    gridrow='which_exp',
    data=df_toplt,
    x_name=x_name,
    y_name=y_name,
    units=units,
    height=3,
    aspect=1,
)
plt.savefig(f"{fig_dir}/goodFit {y_name}X{x_name}.pdf", format='PDF')
# %%
df_toplt = amp_goodTuning
df_toplt.loc[df_toplt['cond_num']==2, 'exp_cond_ordered'] = '2lesion'

x_name='exp_cond_ordered'
y_name=which_slope
units = 'ROI_id'

p = plt_categorical_grid2(
    gridrow='which_exp',
    # gridcol='peak_cat_cond1',
    data=df_toplt,
    x_name=x_name,
    y_name=y_name,
    units=units,
    height=3,
    aspect=1,
)
plt.savefig(f"{fig_dir}/goodTuning {y_name}X{x_name}.pdf", format='PDF')

# %%
# calculate  change
cond1_tuning = df_toplt.query("cond_num == 1")
cond2_tuning = df_toplt.query("cond_num == 2")
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
        g.set(ylim=[np.percentile(df_change[y_name], .1),np.percentile(df_change[y_name], 99)])
    else:
        g.set(ylim=[-1,1])
    plt.savefig(f"{fig_dir}/slopeChg {y_name}X{x_name}.pdf", format='PDF')

# %%
