# %%
import numpy as np
import pandas as pd
import seaborn as sns 
import scipy.io
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

def set_font_type():
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['savefig.bbox'] = 'tight' 
    
def plt_categorical_grid(
    data:pd.DataFrame, 
    x_name:str, 
    y_name:str, 
    units:str, 
    gridrow:str=None, 
    gridcol:str=None, 
    errorbar=None, 
    sharey=True, 
    sns_palette='colorblind', 
    markertype='d', 
    height=4, 
    aspect=0.8,
    ):
    
    """build upon sns.catplot(), plots mean y_name vs x_name with individual repeats. Repeats are plotted as stripplot if number of repeats are different among groups defined by x_name, otherwise, repeats are connected.
    the errorbar arg requires Seaborn v0.12 to work.

    Args:
        data (pd.DataFrame): Dataframe with one value per "unit". Can be achieved using: data.groupby([x_name, gridcol, gridrow, units]).mean().reset_index()
        x_name (str): Column to plot on the X axis.
        y_name (str): Column to plot on the Y axis. 
        gridrow (str): Categorical variables that will determine the faceting of the grid.
        gridcol (str): Categorical variables that will determine the faceting of the grid.
        units (str): Units for plotting individual repeats. see sns.lineplot() units arg.
        errorbar (optional): Defines type of error bars to plot. see seaborn.catplot errorbar arg. Defaults to None.
        sharey (bool, optional): Whether to share Y axis ticks. Defaults to True.
        sns_palette (str, optional): Color palettes. Defaults to 'colorblind'.
        height (int, optional): Height of the graph in inches. Defaults to 4.
        aspect (float, optional): aspect ratio of the graph. Defaults to 0.8.
    """
    data = data.sort_values(by=x_name)
    
    assert_repeats = len(set(data.groupby([x_name])[units].apply(lambda x: len(x.unique())).values))
        
    g = sns.catplot(
        data = data,
        col = gridcol,
        row = gridrow,
        hue = x_name,
        x = x_name,
        y = y_name,
        kind = 'point',
        sharey = sharey,
        palette = sns_palette,
        errorbar = errorbar,
        markers = [markertype]*len(set(data[x_name])),
        height = height,
        aspect = aspect,
        )
    if assert_repeats == 1:
        g.map(sns.lineplot,x_name,y_name,
            estimator=None,
            units=units,
            data = data,
            sort=False,
            color='grey',
            alpha=0.2,
            zorder=0,
            )
    else:
        g.map(sns.stripplot,x_name,y_name,
            data = data,
            color='lightgrey',
            zorder=0,
            order=data[x_name].unique().sort(),
            )
    g.add_legend()
    sns.despine(offset=10, trim=False)
    return g
# %%
# get the name of all folders under root

root = "/Volumes/LabDataPro/2P nMLF speed/analyzed/230609 fish1"
folder_paths = []
all_conditions = []
for folder in os.listdir(root):
    if folder[0] != '.':
        folder_paths.append(root+'/'+folder)
        all_conditions.append(folder)
        
#%%
dFF = pd.DataFrame()

for condition_idx, folder in enumerate(folder_paths):
    this_cond_dFF = pd.DataFrame()
    repeats = os.listdir(folder)
    repeats.sort()
    for repeat in repeats:
        if repeat[0] != '.':
            exp_path = os.path.join(folder, repeat)
            for sti in np.arange(4):
                sti += 1
                load_dFF_thisSti = scipy.io.loadmat(f'{exp_path}/dFF{sti}.mat')
                this_repres = load_dFF_thisSti[f'dff{sti}']
                t_num, ROI_num = this_repres.shape
                this_repres_format = this_repres.flatten(order="F") # for each ROI, t
                frame_col = list(np.arange(t_num)) * ROI_num
                ROI_col = np.repeat(np.arange(ROI_num), t_num)
                
                this_repres_output = pd.DataFrame(data={
                    "dFF": this_repres_format,
                    "frame": frame_col,
                    "roi": ROI_col,
                    "repeat": repeat,
                })
                this_repres_output = this_repres_output.assign(
                    stimulus = sti,
                    condition = all_conditions[condition_idx].split("_")[1]
                )
                dFF = pd.concat([dFF, this_repres_output])
              
# %%
dFF['unilateral'] = 'l'
dFF.loc[dFF['roi']>=8, 'unilateral'] = 'r'
volume_rate = 11.171 / 9
dFF['time'] = dFF['frame']/volume_rate
sti_map = {
    1:5,
    2:10,
    3:20,
    4:30,
}

dFF['pitch'] = dFF['stimulus'].map(sti_map)
set_font_type()

# %%
g = sns.relplot(
    kind='line',
    data=dFF,
    x='time',
    y='dFF',
    col='condition',
    hue='unilateral',
    row='pitch',
)

# %%

p = sns.relplot(
    kind='line',
    data=dFF.query("condition == 'pre' & stimulus > 1"),
    x='time',
    y='dFF',
    hue='pitch',
    aspect=1,
    height=3,
    errorbar=None
)
filename = "nMLF to nose-up pitch time series.pdf"
plt.savefig(filename,format='PDF')
# %%
# get peak
how_many_frames = 2
dFF_get_peak = dFF.query("frame < @how_many_frames")
dFF_amplitude = dFF_get_peak.groupby(['condition','pitch','repeat','unilateral']).mean().reset_index()

# % plot
g = sns.relplot(
    data=dFF_amplitude,
    x='pitch',
    y='dFF',
    hue='unilateral',
    col='condition',
    kind='line',
    units='repeat',
    col_order=['pre','post'],
    estimator=None,
    sort=False,
    height=3
)
filename = "nMLF Tunings post lesion unilateral raw.pdf"
plt.savefig(filename,format='PDF',)
# %
g = sns.relplot(
    data=dFF_amplitude,
    x='pitch',
    y='dFF',
    hue='condition',
    col='unilateral',
    kind='line',
    facet_kws={'sharey':True},
    estimator='mean',
    height=3,
)

filename = "nMLF Tunings post lesion unilateral.pdf"
plt.savefig(filename,format='PDF',)


# %%
