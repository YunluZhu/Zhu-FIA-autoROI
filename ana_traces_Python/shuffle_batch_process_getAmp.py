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
from functions.getAmp_fitDualSlope import getAmp_fitDualSlope_kdeBaseCond1base


# %%
root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_light"
if_reanalyze = 'y'

# %%

def batch_getAmp_fitDualSlope_wBaseTrials_shuffle(root, if_reanalyze):
    print("Shuffle conditions by fish:")
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
    amp_avg = pd.DataFrame() 
    slope = pd.DataFrame()
    ROI_metadata = pd.DataFrame()
    
    folder_paths.sort()
    if if_reanalyze == 'y':
        for fish_idx, fishfolder in enumerate(folder_paths):
            if (os.path.isdir(fishfolder)) and ('fish' in fishfolder):
                this_traces_avg, this_amp_long, this_amp_avg, this_slope, this_ROI_metadata, STIMULUS = getAmp_fitDualSlope_kdeBaseCond1base(fishfolder, if_shuffle=True)
                traces_avg = pd.concat([traces_avg, this_traces_avg])
                amp_long = pd.concat([amp_long, this_amp_long])
                amp_avg = pd.concat([amp_avg, this_amp_avg])
                slope = pd.concat([slope, this_slope])
                ROI_metadata = pd.concat([ROI_metadata, this_ROI_metadata])
    else:
        for fish_idx, fishfolder in enumerate(folder_paths):
            if (os.path.isdir(fishfolder)) and ('fish' in fishfolder):
                this_traces_avg = pd.read_hdf(f'{fishfolder}/dFF_shuffled.h5', key='traces')
                this_amp_long = pd.read_hdf(f'{fishfolder}/dFF_shuffled.h5', key='amp')
                amp_avg = pd.read_hdf(f'{fishfolder}/dFF_analyzed.h5', key='amp_avg')
                this_slope = pd.read_hdf(f'{fishfolder}/dFF_shuffled.h5', key='slope')
                this_ROI_metadata = pd.read_hdf(f'{fishfolder}/dFF_shuffled.h5', key='roi_metadata')
                traces_avg = pd.concat([traces_avg, this_traces_avg])
                amp_long = pd.concat([amp_long, this_amp_long])
                slope = pd.concat([slope, this_slope])
                ROI_metadata = pd.concat([ROI_metadata, this_ROI_metadata])

#%%

    sti_map = dict([(ii+1, sti) for ii, sti  in enumerate(STIMULUS)])
    
    if 'lesion' in root:
        exp1_name = 'prox'
        exp2_name = 'dista'
        exp2id_char = 'distal'   
        
    elif 'light' in root:
        exp1_name = 'nMLF'
        exp2_name = 'TAN'
        exp2id_char = 'S'

    traces_avg = traces_avg.assign(
        ROI_id = traces_avg['fish_id'] + '_' + traces_avg['ROI'].astype(str),
        which_exp = exp1_name,
        stimulus = traces_avg['nsti'].map(sti_map),
    )
    traces_avg.loc[traces_avg['fish_info'].str.contains(exp2id_char), 'which_exp'] = exp2_name

    slope = slope.assign(
        ROI_id = slope['fish_id'] + '_' + slope['ROI'].astype(str),
        which_exp = exp1_name,
    )
    slope.loc[slope['fish_info'].str.contains(exp2id_char), 'which_exp'] = exp2_name

    amp_long = amp_long.assign(
        ROI_id = amp_long['fish_id'] + '_' + amp_long['ROI'].astype(str),
        which_exp = exp1_name,
        stimulus = amp_long['nsti'].map(sti_map),
    )
    amp_long.loc[amp_long['fish_info'].str.contains(exp2id_char), 'which_exp'] = exp2_name

    amp_avg = amp_avg.assign(
        ROI_id = amp_avg['fish_id'] + '_' + amp_avg['ROI'].astype(str),
        which_exp = exp1_name,
    )
    amp_avg.loc[amp_avg['fish_info'].str.contains(exp2id_char), 'which_exp'] = exp2_name

    # %%
    traces_avg.to_hdf(f'{root}/res_shuffled.h5', key='long_data', mode='w', format='table')
    ROI_metadata.to_hdf(f'{root}/res_shuffled.h5', key='roi_metadata', format='table')
    slope.to_hdf(f'{root}/res_shuffled.h5', key='slope', format='table')
    amp_long.to_hdf(f'{root}/res_shuffled.h5', key='amp', format='table')
    amp_avg.to_hdf(f'{root}/res_shuffled.h5', key='amp_avg', format='table')

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

    # try:
    #     if_contain_base_trials
    # except NameError:
    #     if_contain_base_trials = input("- Get baseline from trials or not? (y/n): ")
             
    batch_getAmp_fitDualSlope_wBaseTrials_shuffle(root, if_reanalyze)