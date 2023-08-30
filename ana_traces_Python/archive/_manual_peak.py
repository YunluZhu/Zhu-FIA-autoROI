# %%
import pandas as pd
import numpy as np
import seaborn as sns
 
 
# %%
 
path = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/230609 nMLF speed TY functional/f1/pre_and_post.csv"

manual_df = pd.read_csv(path, header=None)
manual_df.columns = ['id', 'area', 'mean', 'min','max', 'stimulus', 'condition']
manual_df = manual_df.assign(
    trial_id = sum([[ele] * 4 for ele in np.arange(len(manual_df)/4)], [])
)

sns.relplot(
    data = manual_df,
    kind = 'line',
    col = 'condition',
    x = 'stimulus',
    y = 'mean',
    col_order= ['pre', 'post'],
    units='trial_id',
    hue='trial_id',
    estimator=None,
    # sorted=False
)
# %%
path = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/230609 nMLF speed TY functional/f1/post_left_and_right.csv"

manual_df = pd.read_csv(path, header=None)
manual_df.columns = ['id', 'area', 'mean', 'min','max', 'stimulus', 'condition','trial_id']

sns.relplot(
    data = manual_df,
    kind = 'line',
    col = 'condition',
    x = 'stimulus',
    y = 'mean',
    col_order= ['l', 'r'],
    units='trial_id',
    hue='trial_id',
    estimator=None,
    # sorted=False
)
# %%

# %%

# %%

# %%

# %%

# %%

# %%
