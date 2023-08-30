#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.plt_functions import plt_categorical_grid2
from functions.doQC_getSlope import doQC_getSlope, sel_exp
import matplotlib.pyplot as plt
from functions.plt_tools import set_font_type

#%%
sel_dir = 'light' # lesion or light
sel_qc = 'good_tuning'

#%%
set_font_type()
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
all_amp, slope = doQC_getSlope(root, cond4qc)

df = slope.loc[slope['good_tuning']]

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

# %%
# timing

df = all_amp#.loc[all_amp_long['ROI_id'].isin(amp_goodTuning.ROI_id.unique())]

df = df.sort_values(by=['ROI_id','exp_cond_ordered','nsti','repeat']).reset_index(drop=True)

df_oneTimePersti = df.groupby(['which_exp','ROI_id','cond_num','exp_cond_ordered','nsti'])[['peak_time_onTrialAvg','peak_time_onAvgTrial']].mean().reset_index()

df_oneTimePersti = df_oneTimePersti.assign(
    peak_time_adj = (df_oneTimePersti['peak_time_onTrialAvg']+df_oneTimePersti['peak_time_onAvgTrial'])/2
)

#%% use sti 10,20 for timing calculation

df_oneTimePersti_sel = df_oneTimePersti.loc[df_oneTimePersti['nsti'].isin([2,3])]

# determine peak timing category
df_oneTimePersti_sel_oneTimePerROI = df_oneTimePersti_sel.groupby(['ROI_id','cond_num', 'nsti'])['peak_time_adj'].mean().reset_index()
df_oneTimePersti_sel_oneTimePerROI = df_oneTimePersti_sel_oneTimePerROI.loc[df_oneTimePersti_sel_oneTimePerROI['cond_num'].isin([1,2])]
df_oneTimePersti_sel_oneTimePerROI = df_oneTimePersti_sel_oneTimePerROI.groupby(['ROI_id', 'nsti'])['peak_time_adj'].max().reset_index()
ROI_id = df_oneTimePersti_sel_oneTimePerROI['ROI_id'].values
ROI_cat = pd.cut(df_oneTimePersti_sel_oneTimePerROI['peak_time_adj'], bins=[-1,0.5,5], labels=['1slow','2fast']).values

ROI_cat_mapper = dict(zip(ROI_id, ROI_cat))
#%%
df_toplt = df_toplt.loc[df_toplt['cond_num'].isin([1,2])]

df_toplt = df_toplt.assign(
    peak_cat = df_toplt['ROI_id'].map(ROI_cat_mapper)
)

df_toplt = df_toplt.sort_values(by=['ROI_id','exp_cond_ordered','stimulus','peak_cat']).reset_index(drop=True)
# %%

peek_cat_map = df_toplt.set_index('ROI_id').to_dict()['peak_cat']

# %%
traces = pd.read_hdf(f'{root}/res_concatenated.h5', key='long_data')
traces = traces.loc[traces['ROI_id'].isin(df_toplt.ROI_id.unique())]
# %%
time_list = [list(np.arange(1/1.2987,20,1/1.2987)+(nsti-1)*20) for nsti in [1,2,3,4]]
time_list = np.array(time_list).flatten().tolist()
time_map = dict(zip(np.arange(0,100).tolist(), time_list))


traces = traces.assign(
    peak_cat = traces['ROI_id'].map(peek_cat_map),
    time = traces['frames'].map(time_map),
    exp_cond_nsti = traces['exp_cond'] + traces['nsti'].astype(str)
)

# %%
g = sns.relplot(
    data=traces.loc[traces['area'].isin([1,2])],
    row='which_exp',
    col='peak_cat',
    x='time',
    y='dFF',
    kind='line',
    hue='exp_cond_nsti',
    height=2.5,
    aspect=2
)
plt.savefig(f"{fig_dir}/traces by peak time.pdf", format='PDF')
# %%
