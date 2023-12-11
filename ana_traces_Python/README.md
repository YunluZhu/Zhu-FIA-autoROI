# README

## Introduction

## Usage

### Extract amplitudes

1. After extracting traces from new .tif files using Matlab code, move folders to the directory with the rest of the dataset.
2. Run `batch_process_getAmp.py`. This generates a `res_concatenated.h5` file under the root directory

### Plot results

1. Use `analyze_plot_amp.py` to visualize tuning
2. Use `analyze_plot_traces.py` to visualize traces
2. use `ROI_categorization_UMAP.py` to categorize neurons and view tuning by clusters
