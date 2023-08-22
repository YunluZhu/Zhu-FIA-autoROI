# %% 
import sys
import os,glob
from _func_analyze_traces import analyze_traces 
from tqdm import tqdm

# %%
root = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/analyzed_UP"
DATA_FILE = 'dFF_ksDensity'
STIMULUS = [5, 10, 20, 30]

fish_count = 0
for _, _, files in os.walk(root):
    fish_count = fish_count + len([ele for ele in files if "RoiSet.zip" in ele])

all_folders = os.walk(root)
with tqdm(total=fish_count) as pbar:  # Do tqdm
    # determine if FIJI ROI is under root folder
    filenames = glob.glob(os.path.join(root,"RoiSet.zip"))
    if filenames:  # if FIJI ROI under root, process them
        print(f"\n\n- In {root}")
        analyze_traces(root, DATA_FILE, STIMULUS)
        pbar.update(len(filenames)) # update progress bar after processing FIJI ROI in the current folder

    for path, dir_list, file_list in all_folders: # look for FIJI ROI in all subfolders
        # loop through each subfolder
        dir_list.sort()
        for folder_name in dir_list:
            # get the folder dir by joining path and subfolder name
            folder = os.path.join(path, folder_name)
            filenames = glob.glob(os.path.join(folder,"RoiSet.zip"))
            if filenames:
                analyze_traces(folder, DATA_FILE, STIMULUS)
# %%
# to concatenate processed files and save as a master outside
