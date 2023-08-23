#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.plt_functions import plt_categorical_grid2
from functions.doQC_getSlope import doQC_getSlope_4LD
import matplotlib.pyplot as plt


#%%

root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_light"
which_slope = 'slopeLowAng'
# which_slope = 'slopeAll_rawAmp'

#%%

amp_QC, amp_goodFit, amp_goodTuning = doQC_getSlope_4LD(root)
fig_dir = f"{root}/figures"
try:
    os.makedirs(fig_dir)
except:
    pass
# %% plot
toplt = amp_goodFit
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
toplt = amp_goodTuning.loc[amp_goodTuning['cond_num'].isin([1,2])]

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
cond1_tuning = amp_goodTuning.query("cond_num == 1")
cond2_tuning = amp_goodTuning.query("cond_num == 2")
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

df_toplt = toplt.assign(
    peak_cat_sti_peakTime = toplt.groupby(['ROI_id','cond_num'])['peak_time'].transform(
        'mean'
    )
)

df_toplt = df_toplt.assign(
    peak_cat_avg = pd.cut(df_toplt['peak_cat_sti_peakTime'], bins=[-1,0.4,3])
)

# %%
x_name='exp_cond_ordered'
y_name=which_slope
units = 'ROI_id'

p = plt_categorical_grid2(
    gridrow='which_exp',
    # gridcol='peak_cat_cond1',
    data=df_toplt,
    x_name=x_name,
    y_name=y_name,
    gridcol='peak_cat_avg',
    units=units,
    height=3,
    aspect=1,
)
plt.savefig(f"{fig_dir}/goodTuning_{y_name}X{x_name}_timing.pdf", format='PDF')