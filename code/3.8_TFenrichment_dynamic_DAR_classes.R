# ===========================================
# Script Title: 3.7 VTA-DN TF enrichment per DAR dynamic class
# Author: Luna Zea Redondo
# Date: 2026-06-09
#
# Description:
#   This script performs TF motif enrichment analysis across dynamic DAR classes
#   derived from k-means clustering. Enrichment is tested separately for DARs
#   classified as transient, recovered, memory or delayed, and split by direction
#   of accessibility change.
# ===========================================

rm(list = ls(all.names = TRUE))
gc()

`%notin%` <- Negate(`%in%`)

# ========== Load Required Libraries ==========
library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(glue)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(data.table)
library(BSgenome.Mmusculus.UCSC.mm10)

set.seed(1)

# ========== Set Working Directory ==========
setwd("/fast/AG_Pombo/luna/2026_rebuttal/12_TFenrichment-memory/")

# ========== Load Inputs ==========
motifs_for_heatmap <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/12_TFenrichment-memory/260610_motifs_for_heatmap.rds")
DAR_complete_results_kmeans <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/11_DARs-TRaCE/260612_DAR_complete_results_kmeans.tsv")

atac <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/10_DARs_TFenrichment-contrast/signac_atac_with_motifs.rds")
DARs_corrected <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/10_DARs_TFenrichment-contrast/260610.DARs_corrected.rds")

# ========== Add Motifs if Missing ==========
DefaultAssay(atac) <- "ATAC"

if (is.null(atac[["ATAC"]]@motifs)) {
  pfm <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/231217.pfm.rds")
  atac <- AddMotifs(object = atac, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)
}

# ========== Prepare Peak IDs and Backgrounds ==========
valid_peaks <- intersect(rownames(atac[["ATAC"]]), rownames(atac[["ATAC"]]@motifs@data))

bdg_allPeaks <- valid_peaks


bdg_allDARs <- DARs_corrected %>%
  dplyr::filter(diff == "Yes") %>%
  pull(peakID) %>%
  gsub(":", "-", .) %>%
  intersect(valid_peaks)

DARs_simp <- DAR_complete_results_kmeans %>%
  mutate(
    peakID = gsub(":", "-", peakID),
    category_name = glue("{memory_class}_{direction}")
  ) %>%
  dplyr::select(category_name, peakID)

unique_categories <- unique(DARs_simp$category_name)

# ========== Run TF Enrichment per Dynamic Class ==========
memory_TFenrichment.complete.df <- data.frame()

for (category in unique_categories) {
  message(glue("Analyzing: {category}"))
  
  selected_dars <- DARs_simp %>%
    dplyr::filter(category_name == category) %>%
    pull(peakID) %>%
    gsub(":", "-", .) %>%
    intersect(valid_peaks)
  
  if (length(selected_dars) == 0) next
  
  selected_dars_allPeaks.tfenrich <- FindMotifs(
    object = atac,
    features = selected_dars,
    background = bdg_allPeaks
  ) %>%
    mutate(category_name = category, bdg = "allPeaks")
  
  selected_dars_allDARs.tfenrich <- FindMotifs(
    object = atac,
    features = selected_dars,
    background = bdg_allDARs
  ) %>%
    mutate(category_name = category, bdg = "allDARs")
  
  memory_TFenrichment.complete.df <- bind_rows(
    memory_TFenrichment.complete.df,
    selected_dars_allPeaks.tfenrich,
    selected_dars_allDARs.tfenrich
  )
}

write_tsv(
  memory_TFenrichment.complete.df,
  "260609_memory_TFenrichment.complete.original_table.tsv"
)

# ========== Filter TF Enrichment Results ==========
memory_TFenrichment.complete.df <- memory_TFenrichment.complete.df %>%
  dplyr::filter(bdg == "allDARs", motif.name %in% motifs_for_heatmap) %>%
  mutate(mlog10_pval = -log10(pvalue))

write_tsv(
  memory_TFenrichment.complete.df,
  "260609_memory_TFenrichment.complete.allDARs.filtered.tsv"
)

# ========== Select Significant TFs ==========
up_categories <- c("transient_up", "recovered_up", "memory_up", "delayed_up")
down_categories <- c("transient_down", "recovered_down", "memory_down", "delayed_down")

upreg_incocaine <- memory_TFenrichment.complete.df %>%
  dplyr::filter(category_name %in% up_categories, pvalue < 0.05) %>%
  distinct(motif.name) %>%
  pull(motif.name)

upreg_insaline <- memory_TFenrichment.complete.df %>%
  dplyr::filter(category_name %in% down_categories, pvalue < 0.05) %>%
  distinct(motif.name) %>%
  pull(motif.name)

# ========== Cocaine-Up Dynamic Classes Heatmap ==========
coc_vs_sal_heatmap_TFs.df <- memory_TFenrichment.complete.df %>%
  dplyr::filter(category_name %in% up_categories, motif.name %in% upreg_incocaine) %>%
  dplyr::select(motif.name, category_name, mlog10_pval)

coc_vs_sal_heatmap_TFs.cast <- data.table::dcast(
  as.data.table(coc_vs_sal_heatmap_TFs.df),
  motif.name ~ category_name,
  value.var = "mlog10_pval",
  fill = 0
) %>%
  as.data.frame() %>%
  column_to_rownames("motif.name")

coc_vs_sal_heatmap_TFs.matrix <- as.matrix(coc_vs_sal_heatmap_TFs.cast)

coc_vs_sal_heatmap_TFs.matrix_trans <- t(coc_vs_sal_heatmap_TFs.matrix)[
  c("transient_up", "recovered_up", "memory_up", "delayed_up"),
  ,
  drop = FALSE
]

dend <- as.dendrogram(hclust(dist(coc_vs_sal_heatmap_TFs.matrix)))
dend <- color_branches(dend, k = min(30, nrow(coc_vs_sal_heatmap_TFs.matrix)))
col_fun <- colorRamp2(c(0, 6), c("white", "red"))

coc_vs_sal.TFenrichment <- Heatmap(
  coc_vs_sal_heatmap_TFs.matrix_trans,
  name = "-log10(pvalue)",
  cluster_rows = FALSE,
  cluster_columns = dend,
  column_split = 10,
  col = col_fun,
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  column_names_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 10)
)

coc_vs_sal.TFenrichment

# ========== Reorder Cocaine-Up Heatmap ==========
col_order <- column_order(coc_vs_sal.TFenrichment)
new_cluster_order <- c(1,2,3,4,6,10,5,9,7,8) 

new_col_order <- col_order[new_cluster_order]
flattened_order <- unlist(new_col_order)
column_groups_cat <- factor(rep(seq_along(new_col_order), lengths(new_col_order)), levels = 1:10)

reordered_matrix <- coc_vs_sal_heatmap_TFs.matrix_trans[, flattened_order, drop = FALSE]

coc_vs_sal_reordered_HM <- Heatmap(
  reordered_matrix,
  column_split = column_groups_cat,
  name = "-log10(pvalue)",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  column_names_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 10)
)

coc_vs_sal_reordered_HM

coc_vs_sal_reordered_HM_vertical <- Heatmap(
  t(reordered_matrix),
  row_split = column_groups_cat,
  name = "-log10(pvalue)",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 10)
)

coc_vs_sal_reordered_HM_vertical

# dev.size()
# width_original = 6
# height_original= 8.4
# 
# plot_name <- "chapter3_atac_TFs_memory_upregulated.pdf"
# 
# pdf(plot_name, width = width_original, height =height_original)
# coc_vs_sal_reordered_HM_vertical
# dev.off()

# ========== Saline-Up Dynamic Classes Heatmap ==========
sal_vs_coc_heatmap_TFs.df <- memory_TFenrichment.complete.df %>%
  dplyr::filter(category_name %in% down_categories, motif.name %in% upreg_insaline) %>%
  dplyr::select(motif.name, category_name, mlog10_pval)

sal_vs_coc_heatmap_TFs.cast <- data.table::dcast(
  as.data.table(sal_vs_coc_heatmap_TFs.df),
  motif.name ~ category_name,
  value.var = "mlog10_pval",
  fill = 0
) %>%
  as.data.frame() %>%
  column_to_rownames("motif.name")

sal_vs_coc_heatmap_TFs.matrix <- as.matrix(sal_vs_coc_heatmap_TFs.cast)

sal_vs_coc_heatmap_TFs.matrix_trans <- t(sal_vs_coc_heatmap_TFs.matrix)[
  c("transient_down", "recovered_down", "memory_down", "delayed_down"),
  ,
  drop = FALSE
]

dend <- as.dendrogram(hclust(dist(sal_vs_coc_heatmap_TFs.matrix)))
dend <- color_branches(dend, k = min(30, nrow(sal_vs_coc_heatmap_TFs.matrix)))

sal_vs_coc.TFenrichment <- Heatmap(
  sal_vs_coc_heatmap_TFs.matrix_trans,
  name = "-log10(pvalue)",
  cluster_rows = FALSE,
  cluster_columns = dend,
  column_split = 10,
  col = col_fun,
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  column_names_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 10)
)

sal_vs_coc.TFenrichment

# ========== Reorder Saline-Up Heatmap ==========
col_order <- column_order(sal_vs_coc.TFenrichment)
new_cluster_order <- c(5,6,7,8,2,9,10,3,4,1)

new_col_order <- col_order[new_cluster_order]
flattened_order <- unlist(new_col_order)
column_groups_cat <- factor(rep(seq_along(new_col_order), lengths(new_col_order)), levels = 1:15)

reordered_matrix <- sal_vs_coc_heatmap_TFs.matrix_trans[, flattened_order, drop = FALSE]

sal_vs_coc_reordered_HM <- Heatmap(
  reordered_matrix,
  column_split = column_groups_cat,
  name = "-log10(pvalue)",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  column_names_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 10)
)

sal_vs_coc_reordered_HM

sal_vs_coc_reordered_HM_vertical <- Heatmap(
  t(reordered_matrix),
  row_split = column_groups_cat,
  name = "-log10(pvalue)",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 10)
)

sal_vs_coc_reordered_HM_vertical

# dev.size()
# width_original = 6
# height_original= 8.4
# 
# plot_name <- "chapter3_atac_TFs_memory_downregulated.pdf"
# 
# pdf(plot_name, width = width_original, height =height_original)
# sal_vs_coc_reordered_HM_vertical
# dev.off()



# ========== Save Workspace ==========
# save.image("260609_TFenrichment_dynamic_class_final_complete_analysis.rds")