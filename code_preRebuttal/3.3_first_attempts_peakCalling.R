# ===========================================
# Script Title: ATAC Peak Calling (first attempts)
# Author: Luna Zea Redondo
# Date: 2023-09-12
# Description:
#   This script performs ATAC-seq peak calling for DNs
#   using ArchR and MACS2. Multiple peak-calling strategies are explored
#   (cluster-based, condition-based, and threshold-based), but the final
#   implementation calls peaks on all DN ATAC cells as a single group.
#
#   The resulting peak set is used for downstream analyses including
#   differential accessibility, motif enrichment, and regulatory inference.
#
#   Earlier peak-calling attempts are retained as commented code for
#   methodological transparency and reproducibility.
#   
#   Intermediate files (scuh as proj3 and proj4) have been used in some other scripts abd
#   it has been documented accordingly.

# ===========================================

#!/usr/bin/env Rscript

# ========== Set Environment ==========
rm(list = ls())
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230310/site-library"))
set.seed(1)

# ========== Load Required Libraries ==========
library(ArchR)
library(dplyr)
library(glue)

# ========== Set Working Directory ==========
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR/peakCalling"
setwd(dir)

# ========== ArchR Setup ==========
addArchRThreads(threads = 2)
addArchRGenome("mm10")

# ========== Define MACS2 Path ==========
pathToMacs2 <- "/home/lzeared/profiles/r_multiome230310/bin/macs2"

##################################################
# Peak Calling – Exploratory Strategies (Archived)
##################################################

# Calling peaks: ATAC clusters based (no sample)
###################################################
# 
# # Load ATAC UMAP info
# load("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR/230602_ATAC_b2.rds") # proj2
# 
# proj2 <- addCellColData(
#   ArchRProj = proj2,
#   data = as.character(DNs_fromRNA$simpleIdent),
#   cells = as.character(rownames(DNs_fromRNA)),
#   name = "condition"
# )
# 
# proj3 <- addGroupCoverages(
#   ArchRProj = proj2,
#   groupBy = "Clusters",
#   force = TRUE
# )
# 
# proj3 <- addReproduciblePeakSet(
#   ArchRProj = proj3,
#   groupBy = "Clusters",
#   pathToMacs2 = pathToMacs2,
#   force = TRUE
# )
# 
# proj3 <- addPeakMatrix(proj3, force = TRUE)
# saveRDS(
#   proj3,
#   "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR/peakCalling/230602_ATAC_proj3.ATACclusters.rds"
# )

##################################################
# Condition-based peak calling (exploratory)
##################################################

# proj4 <- addGroupCoverages(
#   ArchRProj = proj2,
#   groupBy = "condition",
#   force = TRUE
# )
# 
# proj4 <- addReproduciblePeakSet(
#   ArchRProj = proj4,
#   groupBy = "condition",
#   pathToMacs2 = pathToMacs2,
#   force = TRUE
# )
# 
# proj4 <- addPeakMatrix(proj4)
# saveRDS(
#   proj4,
#   "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR/peakCalling/230602_ATAC_proj4.condition.rds"
# )

##################################################
# Increased peak number / relaxed thresholds
##################################################

# proj4_increased <- readRDS("230602_ATAC_proj4.condition.rds")
# proj4_increased <- addReproduciblePeakSet(
#   ArchRProj = proj4_increased,
#   groupBy = "condition",
#   pathToMacs2 = pathToMacs2,
#   force = TRUE,
#   maxPeaks = 200000
# )
# saveRDS(
#   proj4_increased,
#   "230824_ATAC_proj4_increased.condition.rds"
# )

##################################################
# Final Selected Strategy: All DN Cells
##################################################

# Load ArchR project
proj4_increased <- readRDS(
  "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR/peakCalling/230602_ATAC_proj4.condition.rds"
)

# Assign all DN cells to one group
proj4_increased$cellType <- "VTA-DNs"

# Create pseudo-bulk coverage
proj4_increased <- addGroupCoverages(
  ArchRProj = proj4_increased,
  groupBy = "cellType",
  force = TRUE
)

# Call reproducible peaks with relaxed maxPeaks
proj4_increased <- addReproduciblePeakSet(
  ArchRProj = proj4_increased,
  groupBy = "cellType",
  maxPeaks = 400000,
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

# Save final peak-called project
saveRDS(
  proj4_increased,
  "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR/peakCalling/230912_ATAC_proj4_increased.condition.rds"
)

# End of script

