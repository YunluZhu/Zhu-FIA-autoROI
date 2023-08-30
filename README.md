# Analyzing vol functional imaging data using auto ROI detection

## Introduction
- This document describes YZ's pipeline for analyzing volumetric functional imaging data generated using Schoppik 2P system. 
- Currently it is soft coded for Z stacks but hardcoded for the number of stimulus per trial (4) and the duration of each response (20s). 
- This pipeline is optimized for images with too many ROIs to circle by hand. Be noted that during capturing, make sure there's minimum z shift among all data for one fish. Namely, there should be little to no change in z from repeat to repeat, condition (area) to condition, etc, for a given fish. X and Y shift less than 20 pix will be corrected. 
- All incomplete/truncated/unwanted data should be removed before running the pipeline

## Outline
1. imagej: Read image sequence or OME Tiff from 2P dataset. Concatenate into hyperstack
2. ImageJ: Separate hyperstack into stacks of individual Z slices
3. ImageJ:Concatenate Z slices across repeats and then areas, generating one stack per Z per fish
4. Matlab: Run motion correction on all z stacks using EZcalcium
5. ImageJ: Detect ROIs by setting threshold and using particle analyzer
6. Matlab: Get traces for ROIs and calculate dFF
7. Python: Analyze dFF, merge ROIs with similar response and close approximity. Plot traces
8. Python: Further analysis such as Amplitude and Tuning

## Description
1. Generate Hyperstacks
	+  Data structure is rigid. Under the root dir there should be multiple fish folders in the format of "xxxxxx_fish# description" where xxxxxx is 6 digit date of experiment, # is fish number in any number of digits followed by a space and then description of any kind. These information will be used to write into fish metadata later in step 6
	+ Inside each fish folder there should be multiple area folders named as "area#\_description"
	+ each manipulation/condition is called an "area" (due to historic reasons)
	+ Inside each are folder, there should be multiple experiment repeat folders starting with "f"
	+ Here's a sample file structure after complete analysis
2. Split Hyperstacks
	+ splited hyperstacks will be saved inside each repeat folder
3. Concatenate across conditions
	+ Stacks of single z will be concatenated across repeats first, saved under area folders
	+ and across areas, saved under fish folders
4. Motion corrections
	+ Max Shift 20. Initial Batch Size 200. Bin Width 200. Upsampling Factor 50.
5. ROI detection
	+ Motion corrected Z stacks are first concatenated into a master hyperstack. This hyperstack contains ALL data for this fish and will be used for extraccting traces.
	+ An averaged image across T is generated for ROI detection
	+ Gamma is first set to 0.5 to reveal dim cells. Thresholds are determined on a stack by stack base, using Triangle. ROIs are analyzed using particle analylzer after watershed and saved as a zip under each fish folder `run("Gamma...", "value=0.5 stack"); setAutoThreshold("Triangle dark no-reset stack");`
6. Read ROIs and get traces
	+ The Matlab script saves 1 metadata for each analysis, under root dir with fish folders, named by the date of analysis, containing information of the number Here's an example of the fish metadata file
	+ It also saves an ROI metadata file under fish folders, containing names and positions of each ROI
	+ The script saves a few .mat data files, including baseline values, raw fluo, raw fluo dataframe containing additional information, dFF
	+ In addition, it saves a set of adjusted data which are normalized using the baseline period across different areas to account for manipulations that may change background noise. Can be turned off in the script.
7. ROI QC and dFF analysis
	+ First, very similar ROIs that likely come from the same neuron in different slices and show a correlation of response > 0.85 are removed
	+ Plot time series data. User can define which areas to plot.
	
## Analyze and Visualize Results Using Python
Refer to the README in folder `ana_traces_Python`