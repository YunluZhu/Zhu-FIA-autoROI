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

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap
from sklearn.cluster import DBSCAN
#%%
sel_dir = 'light' # lesion or light
sel_qc = 'good_endRes'

#%%
fig_root = f"/Users/yunluzhu/Documents/Lab2/caiman/Volumetric_code/YZ_nMLF_speed/figures/UMAP"
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
amp_cat, cat_col =  get_peakTimingCat(amp_smval)

# %%
amp_cat = amp_cat.loc[amp_cat[sel_qc]]
# amp_cat = amp_cat.loc[amp_cat['which_exp']=='nMLF']

value_col = ['amp_smval', 
            #  'half_decay_time', 
             'peak_time_smval']
ROI_wide = amp_cat.query("cond_num == 1").pivot(index='ROI_id', columns=['stimulus'], values=value_col)
ROI_wide = ROI_wide.reset_index()
ROI_wide.columns = ["_".join(map(str, tup)) for tup in ROI_wide.columns.to_flat_index()]
ROI_wide.rename(columns={'ROI_id_': 'ROI_id'}, inplace=True)

ROI_wide.drop(columns=['peak_time_smval_0'], inplace=True)
try:
    ROI_wide.drop(columns=['half_decay_time_0'], inplace=True)
    ROI_wide.drop(columns=['half_decay_time_1'], inplace=True)
except:
    pass

roi_list = ROI_wide['ROI_id'].values

df = ROI_wide.drop(columns=['ROI_id'])
df_std = StandardScaler().fit_transform(df)#.drop(index=bout_feature[bout_feature['to_bout'].isna()].index))



# %% UMAP

exp_map = dict(zip(amp_cat['ROI_id'], amp_cat['which_exp']))
timing_map = dict(zip(amp_cat['ROI_id'].values, amp_cat[cat_col].values))

standard_embedding = umap.UMAP().fit_transform(df_std)
umap_toplt = ROI_wide.assign(
    umap1 = standard_embedding[:, 0],
    umap2 = standard_embedding[:, 1],
    peak_cat = ROI_wide['ROI_id'].map(timing_map),
    which_exp = ROI_wide['ROI_id'].map(exp_map)
)
# #%%
# g = sns.scatterplot(umap_toplt, x='umap1', y='umap2', hue='peak_cat', alpha=0.5, size=0.1)
# p = sns.relplot(kind='scatter', hue='peak_cat', data= umap_toplt, x='umap1', y='umap2', alpha=0.5, size=0.1)

# %% re umap for clustering 
clusterable_embedding = umap.UMAP(
    n_neighbors=30,
    min_dist=0.1,
    n_components=2,
    random_state=40,
).fit_transform(df_std)


umap_toplt = umap_toplt.assign(
    umapC1 = clusterable_embedding[:, 0],
    umapC2 = clusterable_embedding[:, 1],
)

# g = sns.scatterplot(umap_toplt, x='umapC1', y='umapC2', hue='peak_cat', alpha=0.5, size=0.1)

# % cluste4ring
# cluster_labels = hdbscan.HDBSCAN(
#     min_samples=200,
#     min_cluster_size=5,
# ).fit_predict(clusterable_embedding)

get_clusters = DBSCAN(eps=0.5, 
                      min_samples=12
                      ).fit_predict(clusterable_embedding)

umap_toplt = umap_toplt.assign(
    cluster = get_clusters
)
p = sns.relplot(kind='scatter', hue='cluster', data=umap_toplt, x='umapC1', y='umapC2', alpha=0.5, linewidth=0, 
                palette = 'Set2',
                height=4,
                )
plt.savefig(f"{fig_dir}/UMAP scatter.pdf", format='PDF')

p = sns.relplot(kind='scatter', hue='which_exp', data=umap_toplt, x='umapC1', y='umapC2', alpha=0.5, linewidth=0, 
                palette = 'Set2',
                height=4,
                )
plt.savefig(f"{fig_dir}/UMAP scatter_byExp.pdf", format='PDF')
# %%


# %%
cluster_map = dict(zip(umap_toplt['ROI_id'], umap_toplt['cluster']))
df_toplt = amp_cat.assign(
    cluster = amp_cat['ROI_id'].map(cluster_map)
)
df_toplt = df_toplt#.loc[df_toplt['good_fit']]
# df_toplt.dropna(axis=0, inplace=True)
df_toplt = df_toplt.loc[df_toplt['cond_num'].isin([1,2])]

# g = sns.catplot(df_toplt.query("cluster != -1"), 
#                 x='stimulus', 
#                 y='peak_time_smval', 
#                 col = 'cond_num', 
#                 row = 'which_exp',
#                 hue='cluster', 
#                 kind='point',
#                 height=3)
# plt.savefig(f"{fig_dir}/UMAP cluster peakTiming.pdf", format='PDF')


# g = sns.catplot(df_toplt.query("cluster != -1"), 
#                 x='stimulus', 
#                 y='half_decay_time', 
#                 col = 'cond_num', 
#                 row = 'which_exp',
#                 hue='cluster', 
#                 kind='point',
#                 height=3)
# plt.savefig(f"{fig_dir}/UMAP cluster decay.pdf", format='PDF')


p = sns.catplot(df_toplt.query("cluster != -1"), 
                x='stimulus', 
                y='amp_smval', 
                hue = 'cond_num', 
                row = 'which_exp',
                col='cluster', 
                kind='point',
                height=3,
                sharey=True)
plt.savefig(f"{fig_dir}/UMAP cluster amp.pdf", format='PDF')

# q = sns.catplot(df_toplt.query("cluster != -1"), 
#                 x='stimulus', 
#                 y='amp_smval', 
#                 col = 'cond_num',                 
#                 row = 'which_exp',
#                 hue='cluster', 
#                 kind='point',
#                 height=3)
# %%
g = sns.pointplot(df_toplt.groupby(['ROI_id'])[['peak_time_smval','cluster']].mean().reset_index(),
                x='cluster', 
                y='peak_time_smval', 
                join=False
                )
plt.savefig(f"{fig_dir}/UMAP cluster timing.pdf", format='PDF')

# %%
