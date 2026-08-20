#!/usr/bin/env Rscript
# ===========================================
# Script Title: 3.1. ATAC Arrow Creation + RNA Integration
# Author: Luna Zea Redondo
# Date: 2026-05-20
# Description:
#   Creates ArchR Arrow files from fragment files of 18 ATAC-seq samples
#   (multiome), applies QC filters, adds gene score and tile matrices,
#   imports raw 10X RNA matrices, and adds them to the ArchR project.
# ===========================================

# ========== Set Environment ==========
rm(list = ls())
options(stringsAsFactors = FALSE)
`%notin%` <- Negate(`%in%`)
set.seed(1)

# Do not set .libPaths() manually. Use the Guix R environment only.

# ========== Load Required Libraries ==========
suppressPackageStartupMessages({
  library(parallel)
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
})

# ========== Set Paths ==========
data_dir <- "/data/pombo/Luna/MultiomeCocaineTreatments/data_dir"
project_dir <- "/fast/AG_Pombo/luna/2026_rebuttal/8_ATAC-preprocessing"
output_dir <- file.path(project_dir, "260520_all_MCT")
qc_dir <- file.path(project_dir, "QualityControl")

setwd(project_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

tmp_dir <- "/tmp/lzeared_archr_tmp"
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(TMPDIR = tmp_dir)

# ========== ArchR Setup ==========
addArchRThreads(threads = 1)
addArchRGenome("mm10")
message(glue("How many threads are we using: {getArchRThreads()}"))

# ========== Sample IDs ==========
sample_ids <- c(
  "m30_cocaine_R1", "h1_saline_R1", "h1_cocaine_R1", "h1_cocaine_R2",
  "h4_saline_R1", "h4_cocaine_R1", "h4_cocaine_R2",
  "h8_saline_R1", "h8_cocaine_R1", "h8_cocaine_R2",
  "h24_saline_R1", "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3",
  "d14_saline_R1", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3"
)

fragment_files <- setNames(
  file.path(data_dir, paste0(sample_ids, ".atac_fragments.tsv.gz")),
  sample_ids
)

rna_files <- setNames(
  file.path(data_dir, paste0(sample_ids, ".raw_feature_bc_matrix.h5")),
  sample_ids
)

# ========== Sanity Checks ==========
cat("\nFragment file check:\n")
print(data.frame(
  sample = names(fragment_files),
  file = unname(fragment_files),
  exists = file.exists(unname(fragment_files)),
  stringsAsFactors = FALSE
))

cat("\nRNA file check:\n")
print(data.frame(
  sample = names(rna_files),
  file = unname(rna_files),
  exists = file.exists(unname(rna_files)),
  stringsAsFactors = FALSE
))

if (!all(file.exists(unname(fragment_files)))) {
  stop("Some ATAC fragment files are missing or inaccessible.")
}

if (!all(file.exists(unname(rna_files)))) {
  stop("Some RNA HDF5 files are missing or inaccessible.")
}




# ========== 1. Create Missing Arrow Files ==========

for (sid in names(fragment_files)) {
  
  arrow_path <- file.path(project_dir, paste0(sid, ".arrow"))
  
  # Skip existing Arrow files
  if (file.exists(arrow_path)) {
    message("Arrow already exists: ", sid)
    next
  }
  
  message("\n===== Processing: ", sid, " =====")
  
  one_file <- fragment_files[[sid]]
  
  createArrowFiles(
    inputFiles = one_file,
    sampleNames = sid,
    outputNames = sid,
    QCDir = qc_dir,
    minTSS = 0,
    minFrags = 325,
    maxFrags = 1e6,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    force = TRUE,
    threads = 1,
    subThreading = FALSE
  )
  
  message("Finished: ", sid)
}

# Collect ALL Arrow files (existing + newly created)
ArrowFiles <- file.path(
  project_dir,
  paste0(sample_ids, ".arrow")
)

# Check that all Arrow files exist
print(ArrowFiles)
print(file.exists(ArrowFiles))

if (!all(file.exists(ArrowFiles))) {
  stop("Some Arrow files are still missing!")
}

message("\nAll Arrow files are present.")


# ========== 2. Create ArchR Project ==========
proj_raw <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = output_dir,
  copyArrows = FALSE
)

# ========== 3. Add Gene Expression Matrices ==========
seRNA <- import10xFeatureMatrix(
  input = unname(rna_files),
  names = names(rna_files)
)

if (!is.list(seRNA)) {
  seRNA <- list(seRNA)
}

rna_counts <- do.call(cbind, lapply(seRNA, SummarizedExperiment::assay))
rna_ranges <- SummarizedExperiment::rowRanges(seRNA[[1]])

seRNA2 <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = rna_counts),
  rowRanges = rna_ranges
)

proj_raw <- addGeneExpressionMatrix(
  input = proj_raw,
  seRNA = seRNA2,
  threads = 1
)

# ========== 4. Save Project ==========
proj_raw <- saveArchRProject(
  ArchRProj = proj_raw,
  outputDirectory = output_dir,
  overwrite = TRUE,
  load = TRUE
)

saveRDS(proj_raw, "250520_VTA-DNs_proj_raw.rds")
#save.image("260520.proj-ATAC_VTA-DNs.rds")