'''
nMLF pipeline visualization
prerequisite: fluorescent traces extracted by MATLAB script. 
HARD CODED for DURATION of STIMULUS 20s
for details of the whole pipeline, see: https://benchling.com/s/etr-GCcRV7cmhRg7u4JQuDfr?m=slm-n9tOEQ0aekst3JyASOkw

NOTE reads dF_ksDensity and uses Dark KS baseline (which includes DarkBase) to calculate dFF for all conditions
NOTE slopeAll on raw data with repeats
-----todo-----
NOTE under light, response to lower angles are reduced. Compare amplituds
'''

# %%
import pandas as pd
import numpy as np
import seaborn as sns
import os
import scipy.io
import matplotlib.pyplot as plt
import itertools
from functions.plt_functions import plt_categorical_grid2
import plotly.express as px
# import plotly.graph_objs as go


# if_plot = False

STIMULUS = [5, 10, 20, 30]
nsti = len(STIMULUS)

ori_base = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT 230822 copy with ori kde baseline/230808_fish1 TNY/baseline_KsDensity.mat"
mat = scipy.io.loadmat(ori_base)
baseline_file_name = list(mat.keys())[-1]
baseline_ori = pd.DataFrame(data=mat[baseline_file_name])

new_base = "/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_LT/230808_fish1 TNY/baseline_KsDensity.mat"
mat = scipy.io.loadmat(new_base)
baseline_file_name = list(mat.keys())[-1]
baseline_new = pd.DataFrame(data=mat[baseline_file_name])

