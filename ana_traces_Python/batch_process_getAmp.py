'''
1. process all fish, extract peak amp, calculate slope. See function: getAmp_fitDualSlope_wBaseTrials
2. for each fish, save data under fish folder
3. save concatenated data files under root folder
'''

#%%
import os,glob
import pandas as pd
import numpy as np 
import math
from getAmp_fitDualSlope_wBaseTrials import getAmp_fitDualSlope_wBaseTrials

# %%
root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT"
# if_reanalyze = 'n'
# %%

def batch_getAmp_fitDualSlope_wBaseTrials(root, if_reanalyze):
    STIMULUS = [5, 10, 20, 30]
    folder_paths = []
    all_fish = []
    for folder in os.listdir(root):
        if folder[0] != '.':
            folder_paths.append(root+'/'+folder)
            all_fish.append(folder)
    # %%
    traces_avg = pd.DataFrame()
    amp_long = pd.DataFrame() 
    slope = pd.DataFrame()
    ROI_metadata = pd.DataFrame()
    
    folder_paths.sort()
    if if_reanalyze == 'y':
        for fish_idx, fishfolder in enumerate(fishfolder):
            if (os.path.isdir(fishfolder)) and ('fish' in folder_paths):
                this_traces_avg, this_amp_long, this_slope, this_ROI_metadata, STIMULUS = getAmp_fitDualSlope_wBaseTrials(fishfolder)
                traces_avg = pd.concat([traces_avg, this_traces_avg])
                amp_long = pd.concat([amp_long, this_amp_long])
                slope = pd.concat([slope, this_slope])
                ROI_metadata = pd.concat([ROI_metadata, this_ROI_metadata])
    else:
        for fish_idx, fishfolder in enumerate(folder_paths):
            if (os.path.isdir(fishfolder)) and ('fish' in fishfolder):
                this_traces_avg = pd.read_hdf(f'{fishfolder}/dFF_analyzed.h5', key='traces')
                this_amp_long = pd.read_hdf(f'{fishfolder}/dFF_analyzed.h5', key='amp')
                this_slope = pd.read_hdf(f'{fishfolder}/dFF_analyzed.h5', key='slope')
                this_ROI_metadata = pd.read_hdf(f'{fishfolder}/dFF_analyzed.h5', key='roi_metadata')
                traces_avg = pd.concat([traces_avg, this_traces_avg])
                amp_long = pd.concat([amp_long, this_amp_long])
                slope = pd.concat([slope, this_slope])
                ROI_metadata = pd.concat([ROI_metadata, this_ROI_metadata])

#%%
    nsti = len(STIMULUS)
    sti_map = dict([(ii+1, sti) for ii, sti  in enumerate(STIMULUS)])

    traces_avg = traces_avg.assign(
        ROI_id = traces_avg['fish_id'] + '_' + traces_avg['ROI'].astype(str),
        which_exp = 'nMLF',
        stimulus = traces_avg['nsti'].map(sti_map),
    )
    traces_avg.loc[traces_avg['fish_info'].str.contains('S'), 'which_exp'] = 'TAN'

    slope = slope.assign(
        ROI_id = slope['fish_id'] + '_' + slope['ROI'].astype(str),
        which_exp = 'nMLF',
    )
    slope.loc[slope['fish_info'].str.contains('S'), 'which_exp'] = 'TAN'

    amp_long = amp_long.assign(
        ROI_id = amp_long['fish_id'] + '_' + amp_long['ROI'].astype(str),
        NS = 'nMLF',
        stimulus = amp_long['nsti'].map(sti_map),
    )
    amp_long.loc[amp_long['fish_info'].str.contains('S'), 'which_exp'] = 'TAN'

    # %%
    traces_avg.to_hdf(f'{root}/res_concatenated.h5', key='long_data', mode='w', format='table')
    ROI_metadata.to_hdf(f'{root}/res_concatenated.h5', key='roi_metadata', format='table')
    slope.to_hdf(f'{root}/res_concatenated.h5', key='slope', format='table')
    amp_long.to_hdf(f'{root}/res_concatenated.h5', key='amp', format='table')
    
#%%
if __name__ == "__main__":
    try:
        print(f'- Directory: {root}')
    except:
        root = input("- Where's the root folder?: ")

    try:
        if_reanalyze
    except NameError:
        if_reanalyze = input("- Reanalyze ROIs? (y/n): ")
        
    batch_getAmp_fitDualSlope_wBaseTrials(root, if_reanalyze)