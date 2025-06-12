# Analysis of single-cell chromatin and transcriptional dynamics following drug exposure

Welcome! This GitHub repository is complementary to my PhD thesis and our manuscript (currently in preparation) exploring how a single exposure to cocaine alters chromatin accessibility and gene expression in dopamine neurons (DNs) over time. Using paired single-nucleus ATAC-seq and RNA-seq from the ventral tegmental area (VTA), we investigate the temporal dynamics of neuronal activation and long-term molecular memory.

All the figures in the manuscript and extended data can be reproduced using the scripts in this repository.

If you have any questions, comments, or suggestions for improvement, feel free to reach out by opening an issue.


## Overview of the repository

This repository contains the following folders:

- **`code`**: R scripts, bash scripts, and Jupyter Notebooks used for preprocessing, integration, and downstream analyses.
- **`data`**: Processed input files required for the analysis. Raw data should be downloaded separately from GEO.
- **`figures`**: Empty by default; populated when scripts are executed. Stores generated figures for the manuscript and supplementary materials.
- **`results`**: Contains non-figure outputs such as gene lists, differential analysis results, motif enrichments, and cluster assignments.
