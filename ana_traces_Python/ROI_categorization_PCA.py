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
from sklearn.manifold import TSNE
from sklearn.cluster import (KMeans, SpectralClustering)

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
_, amp_smval, slope = doQC_getSlope_wTimedAvgAmp(root, cond4qc)
amp_cat, cat_col =  get_peakTimingCat(amp_smval)

# %%
value_col = ['amp_smval', 'half_decay_time', 'peak_time_smval']
ROI_wide = amp_cat.query("cond_num == 1").pivot(index='ROI_id', columns=['stimulus'], values=value_col)
ROI_wide = ROI_wide.reset_index()
ROI_wide.columns = ["_".join(map(str, tup)) for tup in ROI_wide.columns.to_flat_index()]
ROI_wide.rename(columns={'ROI_id_': 'ROI_id'}, inplace=True)
ROI_wide.drop(columns=['half_decay_time_0', 'peak_time_smval_0'], inplace=True)
roi_list = ROI_wide['ROI_id'].values

df = ROI_wide.drop(columns=['ROI_id'])
df_std = StandardScaler().fit_transform(df)#.drop(index=bout_feature[bout_feature['to_bout'].isna()].index))
pca = PCA(n_components=10)
principalComponents = pca.fit_transform(df_std)
PCA_components = pd.DataFrame(principalComponents)

# %%

pca = PCA(n_components=4)
pca_result = pca.fit_transform(df_std)

df_toplt = ROI_wide.assign(
    pca1 = principalComponents[:,0],
    pca2 = principalComponents[:,1],
    pca3 = principalComponents[:,2],
    pca4 = principalComponents[:,3],
    peak_cat = ROI_wide['ROI_id'].map(dict(zip(amp_cat['ROI_id'].values, amp_cat[cat_col].values)))
)
print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))
df_toplt = df_toplt
# # plot 2pc
plt.figure(figsize=(5,5))
sns.scatterplot(
    x="pca1", 
    y="pca2",
    data=df_toplt,
    legend="full",
    hue='peak_cat',
    alpha=0.2
)
# %%

# %%
