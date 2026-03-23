# Analysis of single-cell chromatin and transcriptional dynamics following drug exposure

Welcome! This GitHub repository is complementary to my PhD thesis and our manuscript (currently in preparation) exploring how a single exposure to cocaine alters chromatin accessibility and gene expression in dopamine neurons (DNs) over time. Using paired single-nucleus ATAC-seq and RNA-seq from the ventral tegmental area (VTA), we investigate the temporal dynamics of neuronal activation and long-term molecular memory.

All the figures in the manuscript and extended data can be reproduced using the scripts in this repository.

This repository contains code used for the analyses in the thesis and manuscript. It is shared for transparency and reproducibility. Please note that the code was developed for research purposes and is not intended to be production-ready or optimized for general use.

If you have any questions, comments, or suggestions for improvement, feel free to reach out by opening an issue.


## Overview of the repository

Currently, only the code is available (remaining folders will be populated following publication).

This repository contains the following folders:

- **`code`**: R scripts, bash scripts, and Jupyter Notebooks used for preprocessing, integration, and downstream analyses. 
- **`data`**: Processed input files required for the analysis. Raw data should be downloaded separately from GEO.
- **`figures`**: Empty by default; populated when scripts are executed. Stores generated figures for the manuscript and supplementary materials.
- **`results`**: Contains non-figure outputs such as gene lists, differential analysis results, motif enrichments, and cluster assignments.


## Summary of the project:

### Overview  
Experience-driven stimuli can induce lasting cellular memory through molecular mechanisms, and in the context of addiction, these processes are thought to underlie persistent maladaptive changes long after drug exposure. However, the temporal molecular cascades triggered by a single drug exposure remain poorly understood, particularly at the level of specific neuronal populations. In this project, I set out to investigate genome-wide changes in transcription and chromatin accessibility in dopamine neurons (DNs) from the ventral tegmental area (VTA) following a single cocaine exposure, using paired single-nucleus ATAC-seq and RNA-seq across immediate (0–8h), short-term (24h), and long-term (14d) time points.

<p align="center">
  <img src="figures/summary/experimental_design.png" width="100%">
  <br>
  <em>Experimental design and sampling time points following single cocaine exposure</em>
</p>

**Transcriptional responses to cocaine exposure**  
The results show that even a single exposure to cocaine is sufficient to induce widespread and long-lasting transcriptional changes in VTA DNs. In particular, genes that become upregulated are enriched for addiction-associated and psychiatric disorder–related genes, suggesting that early transcriptional activation may already engage pathways linked to long-term vulnerability. Several neuropeptides, including *Cartpt, Ucn, Penk,* and *Npy*, remain upregulated for up to two weeks, pointing to sustained alterations in neuronal signaling. Overall, these findings indicate that transcriptional activation is a dominant feature of the response to cocaine.

<p align="center">
  <img src="figures/summary/gex_heatmap.png" width="100%">
  <br>
  <em>Gene expression patterns of DEGs in each vs. saline comparisons</em>
</p>

**Temporal gene expression dynamics**  
To better understand how these changes evolve over time, differential gene expression patterns were organized into distinct temporal trajectories, including Transient, Recovered, Memory, and Delayed responses. This framework reveals that early transcriptional changes are primarily associated with neuronal activation, plasticity, stress responses, and energy metabolism, whereas later changes reflect processes such as metabolic regulation, synaptic remodeling, and RNA processing. Together, these observations suggest a shift from rapid, activity-dependent responses toward more stable and long-lasting cellular reprogramming.

<p align="center">
  <img src="figures/summary/kmeans_workflow.png" width="100%">
  <br>
  <em>Workflow and classification of gene expression trajectories across time</em>
</p>

**Chromatin accessibility dynamics and transcription factor activity**  
At the level of chromatin accessibility, dynamic changes over time point to distinct regulatory programs underlying these transcriptional responses. Regions that gain accessibility early after cocaine exposure are strongly enriched for AP-1 motifs and pioneer transcription factors, consistent with canonical activity-dependent chromatin remodeling. In contrast, later changes in accessibility are associated with transcription factors involved in oxidative stress regulation and synaptic remodeling. Regions that lose accessibility display more diverse patterns, including transient repression of activity-related motifs and longer-lasting reductions in factors linked to neuronal identity and plasticity, highlighting the complexity of regulatory repression mechanisms.

<p align="center">
  <img src="figures/summary/TFenrichment.png" width="100%">
  <br>
  <em>TF motif enrichment across dynamic chromatin accessibility classes</em>
</p>

**Epigenetic features and regulatory mechanisms**  
Additional analyses provide insight into the epigenetic features associated with these changes. Regions gaining accessibility are enriched for GC-rich sequences, particularly in early and persistent response categories, suggesting a role for CpG-dependent regulatory mechanisms. Furthermore, integration with transcription factor target databases and published ChIP-seq datasets points to a disruption of Polycomb-mediated repression. Genes upregulated following cocaine exposure are enriched for PRC2 targets, and regions gaining accessibility overlap with H3K27me3-marked chromatin, supporting a model in which cocaine exposure leads to the derepression of previously silenced genomic regions.

<p align="center">
  <img src="figures/summary/PRC2_summary.png" width="100%">
  <br>
  <em>A single cocaine exposure induces long-lasting upregulation of Polycomb targets</em>
</p>

**Conclusion**  
Taken together, these findings demonstrate that a single cocaine exposure induces temporally structured transcriptional and epigenomic remodeling in VTA dopamine neurons. The response appears to transition from rapid neuronal activation to more persistent regulatory adaptations, involving both chromatin remodeling and the disruption of repressive mechanisms, which may contribute to long-term vulnerability to addiction.
