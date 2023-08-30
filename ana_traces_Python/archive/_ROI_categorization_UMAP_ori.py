#%%
import os,glob
import pandas as pd
import numpy as np 
import seaborn as sns
import math
from functions.plt_functions import plt_categorical_grid2
from functions.doQC_getSlope import (doQC_getSlope, sel_exp)
from functions.doCategorization import get_peakTime
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap
from sklearn.cluster import DBSCAN
#%%
sel_dir = 'light' # lesion or light
sel_qc = 'good_fit'
if_normalizeAmpByCell = False
filter_before_pca = False
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
df_peakTime, df_decay = get_peakTime(root, cond4qc)

# %%
ROI_timed = df_peakTime.copy()
amp_col = ['amp_1', 'amp_2', 'amp_3', 'amp_4', 'amp_0']
# value_col = ['amp_1', 'amp_2', 'amp_3', 'amp_4', 'amp_0', 'peak_time']

if if_normalizeAmpByCell:
    amp_matrix = ROI_timed[amp_col].T
    amp_norm = amp_matrix.apply(
        lambda x: (x - x.min())/(x - x.min()).max()
    )
    ROI_timed.loc[:, amp_col] = amp_norm.T.values
    ROI_timed = ROI_timed.assign(
        normed_amp = amp_matrix.apply(
            lambda x: x.max()-x.min()
        ).values
    )
    amp_col = ['amp_1', 'amp_2', 'amp_3', 'amp_4', 'amp_0', 'peak_time', 'normed_amp']

value_col = amp_col + ['halfDecT_1', 'halfDecT_2', 'halfDecT_3', 'halfDecT_4']

if filter_before_pca:
    ROI_timed = ROI_timed.loc[ROI_timed['ROI_id'].isin(df_peakTime.loc[df_peakTime[sel_qc],'ROI_id'].unique())]
# %%

reshaped = ROI_timed.loc[ROI_timed['cond_num'].isin([1])].pivot(index='ROI_id', columns='exp_cond_ordered', values=value_col).reset_index()
reshaped.columns = ["_".join(tup) for tup in reshaped.columns.to_flat_index()]
roi_list = reshaped['ROI_id_'].values
df = reshaped.drop(columns=['ROI_id_'])

df_std = StandardScaler().fit_transform(df)


# %% UMAP

standard_embedding = umap.UMAP().fit_transform(df_std)
umap_toplt = reshaped.assign(
    umap1 = standard_embedding[:, 0],
    umap2 = standard_embedding[:, 1],
    peak_cat = reshaped['ROI_id_'].map(dict(zip(ROI_timed['ROI_id'].values, ROI_timed['peak_cat'].values)))
)
# #%%
# g = sns.scatterplot(umap_toplt, x='umap1', y='umap2', hue='peak_cat', alpha=0.5, size=0.1)
# p = sns.relplot(kind='scatter', hue='peak_cat', data= umap_toplt, x='umap1', y='umap2', alpha=0.5, size=0.1)

# %% re umap for clustering 
clusterable_embedding = umap.UMAP(
    n_neighbors=15,
    min_dist=0.1,
    n_components=2,
    random_state=30,
).fit_transform(df_std)


umap_toplt = umap_toplt.assign(
    umapC1 = clusterable_embedding[:, 0],
    umapC2 = clusterable_embedding[:, 1],
)

g = sns.scatterplot(umap_toplt, x='umapC1', y='umapC2', hue='peak_cat', alpha=0.5, size=0.1)

# %% cluste4ring
# cluster_labels = hdbscan.HDBSCAN(
#     min_samples=200,
#     min_cluster_size=5,
# ).fit_predict(clusterable_embedding)

get_clusters = DBSCAN(eps=0.6, 
                      min_samples=15
                      ).fit_predict(clusterable_embedding)

umap_toplt = umap_toplt.assign(
    cluster = get_clusters
)
p = sns.relplot(kind='scatter', hue='cluster', data= umap_toplt, x='umapC1', y='umapC2', alpha=1, 
                # s=2,
                # palette = 'Spectral'
                )


# %%
cluster_map = dict(zip(umap_toplt['ROI_id_'], umap_toplt['cluster']))
df_decay = df_decay.assign(
    cluster = df_decay['ROI_id'].map(cluster_map)
)
g = sns.catplot(df_decay.query("cluster != -1"), 
                x='nsti', 
                y='half_decay_time', 
                col = 'cond_num', hue='cluster', 
                kind='point',
                height=3)

p = sns.catplot(df_decay.query("cluster != -1"), 
                x='nsti', 
                y='amp', 
                hue = 'cond_num', col='cluster', 
                kind='point',
                height=3,
                sharey=True)

q = sns.catplot(df_decay.query("cluster != -1"), 
                x='nsti', 
                y='amp', 
                col = 'cond_num', hue='cluster', 
                kind='point',
                height=3)
# %%
g = sns.pointplot(df_decay.groupby(['ROI_id'])[['peak_time','cluster']].mean().reset_index(),
                x='cluster', 
                y='peak_time', 
                join=False
                )

# %%
