# ===========================================
# Script Title: 230601 ATAC Arrow Creation + RNA Integration
# Author: Luna Zea Redondo
# Date: 2023-06-01
# Description:
#   This script creates ArchR Arrow files from fragment files of 18 ATAC-seq samples 
#   (multiome), applying QC filters and adding gene score and tile matrices.
#   It also loads gene expression matrices from 10X and integrates them into the 
#   ArchR project. This sets up the base project for downstream peak calling and analysis.
#   Recomended to be run in the background. 
# ===========================================

#!/usr/bin/env Rscript

# ========== Set Environment ==========
rm(list = ls())
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_archR/site-library"))
set.seed(1)

# ========== Load Required Libraries ==========
library(Seurat)
library(Signac)
library(ArchR)
library(dplyr)
library(glue)
library(ggplot2)
library(tidyr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(readr)
library(SummarizedExperiment)

# ========== Set Paths ==========
data_dir <- "/data/pombo/Luna/MultiomeCocaineTreatments/data_dir/"
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR"
setwd(dir)

# ========== ArchR Setup ==========
addArchRThreads(threads = 16)
addArchRGenome("mm10")
print(glue("How many threads are we using: {getArchRThreads()}"))

# ========== 0. Create Arrow Files ==========
inputFiles <- c(
  glue("{data_dir}m30_cocaine_R1.atac_fragments.tsv.gz"),
  glue("{data_dir}h1_saline_R1.atac_fragments.tsv.gz"),
  glue("{data_dir}h1_cocaine_R1.atac_fragments.tsv.gz"),
  glue("{data_dir}h1_cocaine_R2.atac_fragments.tsv.gz"),
  glue("{data_dir}h4_saline_R1.atac_fragments.tsv.gz"),
  glue("{data_dir}h4_cocaine_R1.atac_fragments.tsv.gz"),
  glue("{data_dir}h4_cocaine_R2.atac_fragments.tsv.gz"),
  glue("{data_dir}h8_saline_R1.atac_fragments.tsv.gz"),
  glue("{data_dir}h8_cocaine_R1.atac_fragments.tsv.gz"),
  glue("{data_dir}h8_cocaine_R2.atac_fragments.tsv.gz"),
  glue("{data_dir}h24_saline_R1.atac_fragments.tsv.gz"),
  glue("{data_dir}h24_cocaine_R1.atac_fragments.tsv.gz"),
  glue("{data_dir}h24_cocaine_R2.atac_fragments.tsv.gz"),
  glue("{data_dir}h24_cocaine_R3.atac_fragments.tsv.gz"),
  glue("{data_dir}d14_saline_R1.atac_fragments.tsv.gz"),
  glue("{data_dir}d14_cocaine_R1.atac_fragments.tsv.gz"),
  glue("{data_dir}d14_cocaine_R2.atac_fragments.tsv.gz"),
  glue("{data_dir}d14_cocaine_R3.atac_fragments.tsv.gz")
)

names(inputFiles) <- c(
  "m30_cocaine_R1", "h1_saline_R1", "h1_cocaine_R1", "h1_cocaine_R2",
  "h4_saline_R1", "h4_cocaine_R1", "h4_cocaine_R2",
  "h8_saline_R1", "h8_cocaine_R1", "h8_cocaine_R2",
  "h24_saline_R1", "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3",
  "d14_saline_R1", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3"
)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  QCDir = "QualityControl",
  minTSS = 0,
  minFrags = 325,
  maxFrags = 1e6,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)

# Optional backup
# save.image("220930_arrowfiles_backup.rds")

# ========== 1. Create ArchR Project ==========
proj_raw <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "220930_all_MCT",
  copyArrows = TRUE
)

samples <- names(inputFiles)
save.image("230601.projb1.rds")

# ========== 2. Add Gene Expression Matrices ==========
# Note: using raw 10X gene expression matrices from multiome

seRNA <- import10xFeatureMatrix(
  input = c(
    glue("{data_dir}m30_cocaine_R1.raw_feature_bc_matrix.h5"),
    glue("{data_dir}h1_saline_R1.raw_feature_bc_matrix.h5"),
    glue("{data_dir}h1_cocaine_R1.raw_feature_bc_matrix.h5"),
    glue("{data_dir}h1_cocaine_R2.raw_feature_bc_matrix.h5"),
    glue("{data_dir}h4_saline_R1.raw_feature_bc_matrix.h5"),
    glue("{data_dir}h4_cocaine_R1.raw_feature_bc_matrix.h5"),
    glue("{data_dir}h4_cocaine_R2.raw_feature_bc_matrix.h5"),
    glue("{data_dir}h8_saline_R1.raw_feature_bc_matrix.h5"),
    glue("{data_dir}h8_cocaine_R1.raw_feature_bc_matrix.h5"),
    glue("{data_dir}h8_cocaine_R2.raw_feature_bc_matrix.h5"),
    glue("{data_dir}h24_saline_R1.raw_feature_bc_matrix.h5"),
    glue("{data_dir}h24_cocaine_R1.raw_feature_bc_matrix.h5"),
    glue("{data_dir}h24_cocaine_R2.raw_feature_bc_matrix.h5"),
    glue("{data_dir}h24_cocaine_R3.raw_feature_bc_matrix.h5"),
    glue("{data_dir}d14_saline_R1.raw_feature_bc_matrix.h5"),
    glue("{data_dir}d14_cocaine_R1.raw_feature_bc_matrix.h5"),
    glue("{data_dir}d14_cocaine_R2.raw_feature_bc_matrix.h5"),
    glue("{data_dir}d14_cocaine_R3.raw_feature_bc_matrix.h5")
  ),
  names = samples
)

# Combine all gene expression matrices
seRNAcombined <- do.call(cbind, lapply(seRNA, assay))
seRNA2 <- SummarizedExperiment(assays = list(counts = seRNAcombined), rowRanges = rowRanges(seRNA[[1]]))

# Add to ArchR project
proj_raw <- addGeneExpressionMatrix(input = proj_raw, seRNA = seRNA2, threads = 1)

save.image("230601.projb2.rds")
