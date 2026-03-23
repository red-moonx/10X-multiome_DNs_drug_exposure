# ===========================================
# Script Title: 4.3_PRCinATAC_files_generation.R
# Author: Luna Zea Redondo
# Date: 2024-11-13
# Description:
#   Integrates PRC peaks with ATAC-seq peak annotations and DAR tables.
#   Main steps:
#     - load PRC peak table and updated DAR files
#     - format peak coordinates
#     - build GRanges objects for all peaks and PRC peaks
#     - find overlaps between ATAC peaks and PRC peaks
#     - merge PRC annotation into ATAC peak tables
#     - merge PRC annotation with DARs per contrast and memory-class tables
# ===========================================

#!/usr/bin/env Rscript

# ========== Environment setup ==========
rm(list = ls(all.names = TRUE))
gc()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

# ========== Libraries ==========
library(dplyr)
library(tidyr)
library(readr)
library(glue)
library(ggplot2)

library(GenomicRanges)
library(IRanges)

# ========== Project setup ==========
source("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/7_240510_GEX_memory_only/240510_functions_upset.sh")
setwd("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/9_PRCwork_with_ATAC")



# ===========================================
# 1) Input data
# ===========================================

# 1A. Load PRC peaks table
PRC_peaks_table <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/7_midbrain_PRCwork_ATAC/input_files/DN.invivo.WT.K27me3.BCP_peaks.multicov.bed") %>% 
  dplyr::mutate(peakID = glue("{chrom}:{start}-{stop}")) %>% 
  dplyr::select(chrom, start, stop, peakID, lenght, WT_4mon_K27me3_normed, WT4mon_vs_Mut4mon_LFC, WT4mon_vs_Mut4mon_diff, WT4mon_vs_Mut8mon_LFC, WT4mon_vs_Mut8mon_diff)

# 1B. Load ATAC/DAR files and format DARs_memory_classes
#Updated files:
DARs_contrast <- read_tsv("../data/DARs_per_contrast.tsv")
DARs_memory_classes <- read_tsv("../data/DAR_complete_results_kmeans.tsv")
TFmotifs <- readRDS("motif_matrix_reduced_collapsed_DARs_8201x122.rds")

#Format updated files:
DARs_memory_classes <- DARs_memory_classes %>%
  tidyr::separate(peakID, into = c("chr", "range"), sep = ":", remove = FALSE) %>%
  tidyr::separate(range, into = c("start", "end"), sep = "-") %>%
  dplyr::mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  dplyr::select(chr, start, end, peakID = peakID, memory_DARs_kmeans, memory_class, direction)



# ===========================================
# 2) Annotate all peaks with PRC information
# ===========================================
all_peaks <- DARs_contrast %>% dplyr::select(seqnames, start, end, peakID, length, distToGeneStart, nearestGene, peakType, distToTSS, nearestTSS, GC) %>% distinct()

all_peaks_gr <- GRanges(seqnames = all_peaks$seqnames,
                        ranges = IRanges(start = all_peaks$start, end = all_peaks$end),
                        peakID = all_peaks$peakID)

mcols(all_peaks_gr) <- all_peaks %>%
  dplyr::select(-seqnames, -start, -end)  # Exclude seqnames, start, and end as they are already used for GRanges creation

# Convert PRC_peaks_table to GRanges
PRC_gr <- GRanges(seqnames = PRC_peaks_table$chrom,
                  ranges = IRanges(start = PRC_peaks_table$start, end = PRC_peaks_table$stop),
                  peakID = paste0(PRC_peaks_table$chrom, ":", PRC_peaks_table$start, "-", PRC_peaks_table$stop),
                  lenght = PRC_peaks_table$lenght,
                  WT_4mon_K27me3_normed = PRC_peaks_table$WT_4mon_K27me3_normed,
                  WT4mon_vs_Mut4mon_LFC = PRC_peaks_table$WT4mon_vs_Mut4mon_LFC,
                  WT4mon_vs_Mut4mon_diff = PRC_peaks_table$WT4mon_vs_Mut4mon_diff,
                  WT4mon_vs_Mut8mon_LFC = PRC_peaks_table$WT4mon_vs_Mut8mon_LFC,
                  WT4mon_vs_Mut8mon_diff = PRC_peaks_table$WT4mon_vs_Mut8mon_diff)

# Find overlaps between all_peaks and PRC peaks
overlaps <- findOverlaps(all_peaks_gr, PRC_gr)

all_peaks_hits <- queryHits(overlaps)
PRC_hits <- subjectHits(overlaps)

# Create a dataframe with the overlapping PRC peak information
overlap_df <- data.frame(all_peaks_peakID = all_peaks_gr$peakID[all_peaks_hits],
                         PRC_peakID = PRC_gr$peakID[PRC_hits],
                         length = PRC_gr$lenght[PRC_hits],
                         WT_4mon_K27me3_normed = PRC_gr$WT_4mon_K27me3_normed[PRC_hits],
                         WT4mon_vs_Mut4mon_LFC = PRC_gr$WT4mon_vs_Mut4mon_LFC[PRC_hits],
                         WT4mon_vs_Mut4mon_diff = PRC_gr$WT4mon_vs_Mut4mon_diff[PRC_hits],
                         WT4mon_vs_Mut8mon_LFC = PRC_gr$WT4mon_vs_Mut8mon_LFC[PRC_hits],
                         WT4mon_vs_Mut8mon_diff = PRC_gr$WT4mon_vs_Mut8mon_diff[PRC_hits])

# Convert all_peaks_gr back to a data frame to merge with the PRC overlap information
all_peaks_df <- as.data.frame(all_peaks_gr)

# Merge the overlap data with all_peaks, keeping all original columns
all_peaks_withPRC <- all_peaks_df %>%
  left_join(overlap_df, by = c("peakID" = "all_peaks_peakID"))

# View the resulting dataframe
all_peaks_withPRC
#write_tsv(all_peaks_withPRC, "241215_all_peaks_withPRC.tsv")

# ===========================================
# 3) Combine all peaks with DARs per contrast
# ===========================================
DARs_contrast_toMerge <- DARs_contrast %>% dplyr::select(peakID, contrast, logFC, unshrunk.logFC, logCPM, p_val, p_adj.loc, p_adj.glb, diff)
all_peaks_withPRC_withDARs_contrast <- left_join(all_peaks_withPRC, DARs_contrast_toMerge, by = "peakID")
#write_tsv(all_peaks_withPRC_withDARs_contrast, "241215_all_peaks_withPRC_withDARs_contrast.tsv")



# ===========================================
# 4) Combine all peaks with DARs memory classes
# ===========================================
DARs_memory_classes_toMerge <- DARs_memory_classes %>% dplyr::select(peakID, memory_DARs_kmeans, memory_class, direction)
all_peaks_withPRC_withDARs_memory <- left_join(all_peaks_withPRC, DARs_memory_classes_toMerge, by = "peakID")
#write_tsv(all_peaks_withPRC_withDARs_memory, "241215_all_peaks_withPRC_withDARs_memory.tsv")
