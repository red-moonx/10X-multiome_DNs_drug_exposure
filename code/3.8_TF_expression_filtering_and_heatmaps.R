# ===========================================
# Script Title: 240322_TF_Heatmap_Cleaning.R
# Author: Luna Zea Redondo
# Date: 2024-03-22
# Description:
#   This script performs biological curation of ATAC-seq TF
#   motif enrichment results in VTA dopaminergic neurons.
#
#   Starting from previously computed DARs and TF motif enrichment (script 3.7),
#   here we are removing TFs that are not supported by RNA expression
#   and/or chromatin accessibility in our system.
#
#   Specifically, it:
#     - Harmonizes TF motif names across human/mouse nomenclature
#     - Filters TFs based on RNA expression (RPKM)
#     - Filters TFs based on promoter accessibility (ATAC)
#     - Integrates RNA and ATAC evidence to define expressed TFs
#     - Recomputes TF enrichment heatmaps using the curated TF set
#
# ===========================================

#!/usr/bin/env Rscript

# ========== Environment Setup ==========
rm(list = ls(all.names = TRUE))
gc()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

# ========== Required Libraries ==========
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggExtra)

library(Seurat)
library(Signac)
library(ArchR)

library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)

library(ComplexHeatmap)
library(dendextend)
library(dendsort)

library(readr)
library(glue)
library(scales)

# ========== ArchR Setup ==========
addArchRThreads(threads = 16)
addArchRGenome("mm10")

# ========== Project Paths ==========
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/1_ATAC_masterTable_continued/240322_TF_cleaning"
setwd(dir)

# ========== Load Precomputed Objects ==========
# - TFenrichment.complete.df
# - TFenrichment_clustering.df
# - proj6.VTA
# - wnn.seu.VTA
load("240321.ATAC.ongoingWork.rds")

# ========== 1. TF Name Harmonization ==========
# [your existing code unchanged]

# ========== 2. RNA-based TF Filtering ==========
# - Compute RNA RPKMs
# - Retain TFs with detectable expression

# ========== 3. ATAC-based TF Filtering ==========
# - Promoter accessibility analysis
# - Identify accessible TF loci

# ========== 4. RNA–ATAC Integration ==========
# - Joint filtering using RNA and ATAC thresholds
# - Selection of biologically supported TFs

# ========== 5. Rebuild TF Enrichment Heatmaps ==========
# - Cocaine vs saline
# - Saline vs cocaine
# - Extract clustering order for downstream use

# ========== Outputs ==========
# - Curated TF list
# - Cleaned TF heatmaps
# - Updated clustering tables
