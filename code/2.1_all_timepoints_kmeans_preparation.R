# ===========================================
# Script Title: 2.1. K-means analysis prep
# Author: Luna Zea Redondo
# Date: 2026-05-13
# Description:
#   This script processes the VTA-DNs scRNA-seq data to be suitable for
#   k-means analysis (script 2.2) The Seurat object is filtered for VTA-DNs,
#   removes unwanted samples, updates metadata, 
#   and exports the results in H5AD format and as a gene list.
# ===========================================


# ========== Set Environment ==========
rm(list = ls(all.names = TRUE))
gc()
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

source("/fast/AG_Pombo/luna/2026_rebuttal/config.R", local = FALSE)

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
dir <- "/fast/AG_Pombo/luna/2026_rebuttal/5_kmeans_GEX"
setwd(dir)

DEG_complete_results <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/3_muscat/deg_results/260513_DEG_complete_results_Analysis1.tsv")
seu.VTA_DNs <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/2_DN-evaluation/250507_seu.VTA_DNs.rds")
seu.VTA_DNs  # 1520 cells; all coming from VTA

seu.VTA_DNs@meta.data <- seu.VTA_DNs@meta.data %>%
  tidyr::separate(
    orig.ident,
    into = c("timepoint", "treatment", "replicate"),
    sep = "_",
    remove = FALSE
  )

# ========== 1. Filter Samples ==========
# Keep selected salines and cocaine replicates
Idents(seu.VTA_DNs) <- "orig.ident"
samples_totest <- c("h4_saline_R1", "h8_saline_R1", "h24_saline_R1", "d14_saline_R1", #salines
                    "h1_cocaine_R1", "h1_cocaine_R2", "h4_cocaine_R1", "h4_cocaine_R2", "h8_cocaine_R1", "h8_cocaine_R2", #ETP
                    "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3") #LTP

seu.VTA_DNs.2reps <- subset(x = seu.VTA_DNs, idents = samples_totest)



#Harmonize metadata
selected_metadata <- seu.VTA_DNs.2reps@meta.data[, c("timepoint", "treatment", "replicate")]
selected_metadata_mod <- selected_metadata %>% 
  dplyr::mutate(replicate = ifelse(timepoint == "h4" & treatment == "saline", "R1", 
                                   ifelse(timepoint == "h8" & treatment == "saline", "R2",
                                          ifelse(timepoint == "h24" & treatment == "saline", "R3",
                                                 ifelse(timepoint == "d14" & treatment == "saline", "R4", replicate)))))

#Update the Seurat object with the new metadata
seu.VTA_DNs.2reps@meta.data <- selected_metadata_mod

# ========== 2. Convert and Export ==========
# Use DietSeurat and SeuratDisk (keep normalized data; it can be different and could be tested for improvement)
# Delete the slot "scale.data", so Seuratdisk uses "data" as default. 
obj <- seu.VTA_DNs.2reps

# convert Assay5 -> Assay
obj[["RNA"]] <- as(obj[["RNA"]], Class = "Assay")

# keep only counts + data
obj <- DietSeurat(
  obj,
  assays = "RNA",
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE
)

SaveH5Seurat(obj, filename = "260513.DNs.VTA.RNA.seu.sub.data.1387.h5Seurat", overwrite = TRUE)
Convert("260513.DNs.VTA.RNA.seu.sub.data.1387.h5Seurat", dest = "h5ad", overwrite = TRUE)

# ========== 4. Extract DEG Genes ==========
vta.degs <- DEG_complete_results %>% 
  dplyr::filter(significant != "No significant") %>% dplyr::select(gene) %>% distinct()
write_tsv(vta.degs, "260513_vta.degs.txt")
