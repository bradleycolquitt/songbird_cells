# Songbird single-cell/nucleus transcriptomics

This repository contains R code for analysis and figure generation of songbird single-cell RNA-sequencing.

## Directories

comparative - Comparisons between songbird, mouse, and turtle neurons. Scripts are divided three analysis approaches: correlations, scAlign dataset integration, and spatial analyses (Allen Brain Atlas). 

dataset_comparison - Script that compares the distribution of each cell type between HVC and RA

grn - Gene regulatory network inference using GRNBoost2. Included are the scripts to export data for analysis (export_to_numpy_glut.R), script to run the Arboreto framework (run_grnboost2_np.py), and post-processing analysis (grnboost2_glut.R)

marker_gene - Differential expression analysis

reduction_viz - Dimensionality reduction of different subsets of the data (all cells, glutamatergic neurons, GABAergic neurons) and plotting.

trees - Hierarchical clustering of cell clusters

utils - Several files with utility functions used in other scripts
