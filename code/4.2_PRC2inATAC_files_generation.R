# ===========================================
# Script Title: 4.2_PRC2inATAC_files_generation.R
# Author: Luna Zea Redondo
# Date: 2026-05-21
#
# Description:
#   Generates PRC2-annotated ATAC peak tables.
#   Main steps:
#     - load PRC peak table, DAR contrast table, and memory-class DAR table
#     - overlap ATAC peaks with PRC2 peaks
#     - add PRC2 annotation to all ATAC peaks
#     - merge PRC2 annotation with DAR contrast and memory-class tables
# ===========================================

#!/usr/bin/env Rscript

# ========== Environment setup ==========
rm(list = ls(all.names = TRUE)); gc()
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))

# ========== Libraries ==========
library(dplyr)
library(readr)
library(glue)
library(GenomicRanges)
library(IRanges)

# ========== Project setup ==========
setwd("/fast/AG_Pombo/luna/2026_rebuttal/16_DARs_vs_PRC2")

# ===========================================
# 1) Input data
# ===========================================

PRC_peaks_table <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/DN.invivo.WT.K27me3.BCP_peaks.multicov.bed") %>%
  mutate(PRC_peakID = glue("{chrom}:{start}-{stop}")) %>%
  dplyr::select(chrom, start, stop, PRC_peakID, lenght, WT_4mon_K27me3_normed,
         WT4mon_vs_Mut4mon_LFC, WT4mon_vs_Mut4mon_diff,
         WT4mon_vs_Mut8mon_LFC, WT4mon_vs_Mut8mon_diff)

DARs_contrast <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/10_DARs_TFenrichment-contrast/260610.DARs_corrected.rds")
DARs_memory_classes <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/11_DARs-TRaCE/260612_DAR_complete_results_kmeans.tsv") %>%
  mutate(chr = sub(":.*", "", peakID),
         start = as.numeric(sub(".*:(\\d+)-.*", "\\1", peakID)),
         end = as.numeric(sub(".*-", "", peakID))) %>%
  dplyr::select(chr, start, end, peakID, memory_DARs_kmeans, memory_class, direction)

# ===========================================
# 2) Annotate all ATAC peaks with PRC2 information
# ===========================================

all_peaks <- DARs_contrast %>%
  dplyr::mutate(
    seqnames = sub("-.*", "", peakID),
    start = as.numeric(sub("^[^-]+-(\\d+)-.*", "\\1", peakID)),
    end = as.numeric(sub(".*-", "", peakID)),
    length = end - start
  ) %>%
  dplyr::select(seqnames, start, end, peakID, length,
                dplyr::any_of(c("distToGeneStart", "nearestGene", "peakType", "distToTSS", "nearestTSS", "GC"))) %>%
  dplyr::distinct(peakID, .keep_all = TRUE)

all_peaks_gr <- GRanges(seqnames = all_peaks$seqnames,
                        ranges = IRanges(start = all_peaks$start, end = all_peaks$end))
mcols(all_peaks_gr) <- all_peaks %>% dplyr::select(-seqnames, -start, -end)

PRC_gr <- GRanges(seqnames = PRC_peaks_table$chrom,
                  ranges = IRanges(start = PRC_peaks_table$start, end = PRC_peaks_table$stop))
mcols(PRC_gr) <- PRC_peaks_table %>% dplyr::select(-chrom, -start, -stop)

overlaps <- findOverlaps(all_peaks_gr, PRC_gr)

overlap_df <- data.frame(
  peakID = all_peaks_gr$peakID[queryHits(overlaps)],
  PRC_peakID = PRC_gr$PRC_peakID[subjectHits(overlaps)],
  PRC_length = PRC_gr$lenght[subjectHits(overlaps)],
  WT_4mon_K27me3_normed = PRC_gr$WT_4mon_K27me3_normed[subjectHits(overlaps)],
  WT4mon_vs_Mut4mon_LFC = PRC_gr$WT4mon_vs_Mut4mon_LFC[subjectHits(overlaps)],
  WT4mon_vs_Mut4mon_diff = PRC_gr$WT4mon_vs_Mut4mon_diff[subjectHits(overlaps)],
  WT4mon_vs_Mut8mon_LFC = PRC_gr$WT4mon_vs_Mut8mon_LFC[subjectHits(overlaps)],
  WT4mon_vs_Mut8mon_diff = PRC_gr$WT4mon_vs_Mut8mon_diff[subjectHits(overlaps)]
)

all_peaks_withPRC <- as.data.frame(all_peaks_gr) %>%
  left_join(overlap_df, by = "peakID")

# write_tsv(all_peaks_withPRC, "260615_all_peaks_withPRC.tsv")

# ===========================================
# 3) Combine PRC2 annotation with DARs per contrast
# ===========================================

DARs_contrast_toMerge <- DARs_contrast %>%
  dplyr::select(peakID, contrast, logFC, unshrunk.logFC, logCPM,
                p_val, p_adj.loc, p_adj.glb, diff)

all_peaks_withPRC_withDARs_contrast <- all_peaks_withPRC %>%
  dplyr::select(-dplyr::any_of(c("cluster_id", "logFC", "unshrunk.logFC", "logCPM",
                                 "p_val", "p_adj.loc", "p_adj.glb",
                                 "contrast", "diff", "peakID_v2"))) %>%
  dplyr::left_join(DARs_contrast_toMerge, by = "peakID")


# write_tsv(all_peaks_withPRC_withDARs_contrast, "260615_all_peaks_withPRC_withDARs_contrast.tsv")

# ===========================================
# 4) Combine PRC2 annotation with DAR memory classes
# ===========================================

all_peaks_withPRC <- all_peaks_withPRC %>%
  dplyr::mutate(peakID = sub("^(chr[^-]+)-", "\\1:", peakID))

DARs_memory_classes_toMerge <- DARs_memory_classes %>%
  dplyr::select(peakID, memory_DARs_kmeans, memory_class, direction)

all_peaks_withPRC_withDARs_memory <- all_peaks_withPRC %>%
  left_join(DARs_memory_classes_toMerge, by = "peakID")

# write_tsv(all_peaks_withPRC_withDARs_memory, "260615_all_peaks_withPRC_withDARs_memory.tsv")