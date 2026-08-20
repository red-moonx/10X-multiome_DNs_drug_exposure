# ===========================================
# Script Title: 3.5 VTA-DN ATAC DAR k-means memory classification
# Author: Luna Zea Redondo
# Date: 2026-06-09
#
# Description:
#   This script prepares late-condition ATAC DARs for k-means clustering
#   and integrates manually assigned k-means memory classes with transient DARs.
# ===========================================

rm(list = ls(all.names = TRUE))
gc()

`%notin%` <- Negate(`%in%`)

# ========== Load Required Libraries ==========
library(Seurat)
library(Signac)
library(SeuratDisk)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(glue)

source("/fast/AG_Pombo/luna/2026_rebuttal/config.R", local = FALSE)



dir <- "/fast/AG_Pombo/luna/2026_rebuttal/12_DARs-TRaCE/"
setwd(dir)

# ========== 1. Subset ATAC object to selected samples / late conditions ==========
atac <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/9_ATAC-QCandIntegration/260608_signac_atac_peak_matrix_object_annotated.rds")
atac@meta.data <- atac@meta.data %>%
  tidyr::separate(
    sample,
    into = c("timepoint", "treatment", "replicate"),
    sep = "_",
    remove = FALSE
  )

Idents(atac) <- "sample"

#Extract only the "late" cells
atac.late <- subset(atac, subset = condition %in% c("saline", "h24_cocaine", "d14_cocaine") )

# ========== 2. Metadata cleanup and export ==========

selected_metadata <- atac.late@meta.data[, c("timepoint", "treatment", "replicate")]

selected_metadata_mod <- selected_metadata %>%
  dplyr::mutate(
    replicate = ifelse(timepoint == "h4" & treatment == "saline", "R1",
                       ifelse(timepoint == "h8" & treatment == "saline", "R2",
                              ifelse(timepoint == "h24" & treatment == "saline", "R3",
                                     ifelse(timepoint == "d14" & treatment == "saline", "R4", replicate))))
  )

atac.late@meta.data <- selected_metadata_mod

obj <- atac.late

DefaultAssay(obj) <- "ATAC"

# remove motifs to avoid HDF5 error
obj[["ATAC"]]@motifs <- NULL

# convert Assay5 -> Assay if needed
obj[["ATAC"]] <- as(obj[["ATAC"]], Class = "Assay")

obj <- DietSeurat(
  obj,
  assays = "ATAC",
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE
)

# SaveH5Seurat(obj, filename = "260609.DNs.VTA.ATAC.seu.sub.data.late.781cells.h5Seurat", overwrite = TRUE)
# Convert("260609.DNs.VTA.ATAC.seu.sub.data.late.781cells.h5Seurat", dest = "h5ad", overwrite = TRUE)


# =========================
# 3. Save late DARs used for k-means
# =========================
late_contrasts <- c("h24_cocaine-saline", "d14_cocaine-saline", "d14_cocaine-h24_cocaine")
#DARs_corrected <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/11_DARs_TFenrichment-contrast/260610.DARs_corrected.rds")
DARs_corrected <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/11_DARs_TFenrichment-contrast/260612.DARs_corrected.rds")

DARs_corrected_late_toKmeans <- DARs_corrected %>%
  dplyr::filter(diff == "Yes", contrast %in% late_contrasts) %>%
  dplyr::select(peakID) %>%
  distinct()

test <- DARs_corrected %>%
  dplyr::filter(diff == "Yes") %>%
  dplyr::select(peakID) %>%
  distinct()

#write_tsv(DARs_corrected_late_toKmeans, "260612_vta.late.DARs.txt")

# =========================
# 4. Manual k-means class annotation
# =========================

# kdars_recovered <- c(1, 3, 8, 10, 17, 18, 19, 22, 32, 33, 45, 46)
# kdars_memory    <- c(4, 9, 11, 12, 15, 23, 25, 30, 31, 34, 35, 37, 39, 42, 44, 47, 50)
# kdars_delayed   <- c(2, 5, 6, 7, 13, 14, 16, 20, 21, 24, 26, 27, 28, 29, 36, 38, 40, 41, 43, 48, 49)
# 
# kdars_up   <- c(1, 2, 4, 6, 10, 11, 14, 15, 16, 18, 19, 20, 21, 28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 39, 43, 45, 47, 49, 50)
# kdars_down <- c(3, 5, 7, 8, 9, 12, 13, 17, 22, 23, 24, 25, 26, 27, 37, 40, 41, 42, 44, 46, 48)

# kdars_recovered <- c(1, 4, 15, 16, 21, 25, 26, 28, 39, 40, 42, 48, 50)
# kdars_memory <- c(3, 5, 9, 10, 12, 14, 18, 19, 31, 32, 33, 34, 35, 37, 49)
# kdars_delayed <- c(2, 6, 7, 8, 11, 13, 17, 20, 22, 23, 24, 27, 29, 30, 36, 38, 41, 43, 44, 45, 46, 47)
# 
# kdars_up <- c(1, 2, 5, 7, 8, 10, 12, 13, 15, 16, 17, 19, 20, 22, 25, 26, 27, 28, 31, 34, 35, 36, 37, 39, 41, 44, 46, 47, 50)
# kdars_down <- c(3, 4, 6, 9, 11, 14, 18, 21, 23, 24, 29, 30, 32, 33, 38, 40, 42, 43, 45, 48, 49)

# kdars_recovered <- c(1,4,8,15,16,21,25,26,28,39,40,42,48,50)
# kdars_memory    <- c(3,5,9,10,14,18,19,31,32,33,34,35,37,49)
# kdars_delayed   <- c(2,6,7,11,12,13,17,20,22,23,24,27,29,30,36,38,41,43,44,45,46,47)
# 
# kdars_up   <- c(1,2,5,7,10,12,13,15,16,17,19,20,22,25,26,27,28,31,34,35,36,37,39,41,44,46,47,50)
# kdars_down <- c(3,4,6,8,9,11,14,18,21,23,24,29,30,32,33,38,40,42,43,45,48,49)


kdars_recovered <- c(1, 2, 11, 15, 18, 20, 27, 30, 35, 37, 40, 42, 46, 50)
kdars_memory <- c(6, 7, 12, 14, 16, 21, 22, 23, 25, 32, 36, 38, 39, 43, 47, 48)
kdars_delayed <- c(3, 4, 5, 8, 9, 10, 13, 17, 19, 24, 26, 28, 29, 31, 33, 34, 41, 44, 45, 49)

kdars_up <- c(2, 3, 5, 7, 8, 11, 13, 14, 15, 17, 19, 20, 23, 25, 28, 29, 34, 35, 36, 38, 39, 41, 42, 43, 50)
kdars_down <- c(1, 4, 6, 9, 10, 12, 16, 18, 21, 22, 24, 26, 27, 30, 31, 32, 33, 37, 40, 44, 45, 46, 47, 48, 49)

memory_DARs_kmeans <- read_csv("260612_DARs_memory_kmeans_version2_smooth.csv") %>%
  dplyr::rename(memory_DARs_kmeans = "level_0", peakID = "0") %>%
  dplyr::mutate(
    memory_DARs_kmeans = as.numeric(str_replace(memory_DARs_kmeans, "cluster", "")),
    memory_class = case_when(
      memory_DARs_kmeans %in% kdars_recovered ~ "recovered",
      memory_DARs_kmeans %in% kdars_memory ~ "memory",
      memory_DARs_kmeans %in% kdars_delayed ~ "delayed"
    ),
    direction = case_when(
      memory_DARs_kmeans %in% kdars_up ~ "up",
      memory_DARs_kmeans %in% kdars_down ~ "down"
    )
  ) %>%
  dplyr::select(memory_DARs_kmeans, peakID, memory_class, direction)

# =========================
# 5. Define transient DARs
# =========================

all_DARs <- DARs_corrected %>%
  dplyr::filter(diff == "Yes") %>%
  dplyr::select(peakID, contrast, logFC, diff)

all_DARs$contrast <- factor(all_DARs$contrast, levels = all_contrasts, ordered = TRUE)
all_DARs <- all_DARs[order(all_DARs$contrast), ] %>%
  dplyr::mutate(peakID = sub("^(chr[^-]+)-", "\\1:", peakID))

transient_DARs_first_occurrences <- all_DARs %>%
  dplyr::filter(peakID %notin% memory_DARs_kmeans$peakID) %>%
  mutate(direction = ifelse(logFC > 0, "up", "down")) %>%
  group_by(peakID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(
    memory_class = "transient",
    memory_DARs_kmeans = NA
  ) %>%
  dplyr::select(memory_DARs_kmeans, peakID, memory_class, direction)

# =========================
# 6. Combine and save
# =========================

DAR_complete_results_kmeans <- bind_rows(
  transient_DARs_first_occurrences,
  memory_DARs_kmeans
)

# write_tsv(DAR_complete_results_kmeans, "260612_DAR_complete_results_kmeans.tsv")

table(DAR_complete_results_kmeans$memory_class, useNA = "ifany")
table(DAR_complete_results_kmeans$direction, useNA = "ifany")
