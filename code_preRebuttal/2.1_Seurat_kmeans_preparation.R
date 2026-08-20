# ===========================================
# Script Title: 231003 RNA data processing for K-means
# Author: Luna Zea Redondo
# Date: 2025-07-22 (for github)
# Description:
#   This script processes the VTA-DNs scRNA-seq data to be suitable for
#   k-means analysis (script 1.6) The Seurat object is filteredfor VTA-DNs,
#   removes unwanted samples, updates metadata, 
#   and exports the results in H5AD format and as a gene list.
# ===========================================

# ========== Set Environment ==========
rm(list = ls(all.names = TRUE))
gc()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

# ========== Load Required Libraries ==========
library(muscat)
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(glue)
library(ggplot2)
library(readr)
library(SeuratDisk)

# ========== Set Directory + Load Data ==========
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/4_240209_GEX_kmeans"
setwd(dir)

load("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002.RNAonly.excluding_putative_SN_cells.rds")
DNs.RNA.seu  # 2411 cells; "region" separates VTA from SN cells

# ========== 1. Filter Samples ==========
# Keep selected salines and cocaine replicates
Idents(DNs.RNA.seu) <- "orig.ident"
samples_totest <- c("h4_saline_R1", "h8_saline_R1", "h24_saline_R1", "d14_saline_R1", #salines
                    "h1_cocaine_R1", "h1_cocaine_R2", "h4_cocaine_R1", "h4_cocaine_R2", "h8_cocaine_R1", "h8_cocaine_R2", #ETP
                    "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3") #LTP

DNs.RNA.seu.2reps <- subset(x = DNs.RNA.seu, idents = samples_totest)

#Keep only region == "VTA"
DNs.RNA.seu.VTA <- subset(DNs.RNA.seu.2reps, subset = region == "VTA")

#Harmonize metadata
selected_metadata <- DNs.RNA.seu.VTA@meta.data[, c("timepoint", "treatment", "replicate")]
selected_metadata_mod <- selected_metadata %>% 
  dplyr::mutate(replicate = ifelse(timepoint == "h4" & treatment == "saline", "R1", 
                                   ifelse(timepoint == "h8" & treatment == "saline", "R2",
                                          ifelse(timepoint == "h24" & treatment == "saline", "R3",
 
                                                                                                 ifelse(timepoint == "d14" & treatment == "saline", "R4", replicate)))))
#Update the Seurat object with the new metadata
DNs.RNA.seu.VTA@meta.data <- selected_metadata_mod

# ========== 2. Convert and Export ==========
# Use DietSeurat and SeuratDisk (keep normalized data; it can be different and could be tested for improvement)
# Delete the slot "scale.data", so Seuratdisk uses "data" as default. 
DNs.RNA.seu.VTA.data <- DietSeurat(DNs.RNA.seu.VTA, counts = TRUE, data = TRUE, scale.data = FALSE)

# Save as .h5Seurat and convert to .h5ad
SaveH5Seurat(DNs.RNA.seu.VTA.data, filename = "240210.DNs.VTA.RNA.seu.sub.data.1435.h5Seurat")
Convert("240210.DNs.VTA.RNA.seu.sub.data.1435.h5Seurat", dest = "h5ad")

# ========== 4. Extract DEG Genes ==========
vta.degs <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/231024_discrete_updated/231026_muscat/Analysis2/deg_results/231002_DEG_complete_results_Analysis2.tsv") %>% 
  filter(cluster_id == "VTA", significant != "No significant") %>% dplyr::select(gene) %>% distinct()
write_tsv(vta.degs, "231026_vta.degs.txt")
