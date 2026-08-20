# ===========================================
# Script Title: 3.11_TFenrichment_dynamic_DAR_class.R
# Author: Luna Zea Redondo
# Date: 2024-11-21
# Description:
#   This script performs transcription factor (TF) motif enrichment
#   analysis across dynamic chromatin accessibility classes defined
#   from ATAC-seq differential accessible regions (DARs).
#
#   Using DARs grouped into dynamic categories (e.g. transient,
#   recovered, memory, delayed; up- and down-regulated), the script:
#
#     - Runs motif enrichment analysis per dynamic class
#     - Filters enrichment results to biologically supported TF motifs
#     - Generates heatmaps of -log10(p-value) enrichment scores
#     - Computes z-score–normalized heatmaps for comparative visualization
#     - Reorders motif clusters for improved interpretability
#     - Extracts motif subsets of interest (e.g. AP1 variants)
#
#   The goal is to characterize dynamic TF regulatory programs
#   associated with memory-related chromatin accessibility changes
#   in VTA-DNs following cocaine exposure.
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
library(stringr)
library(readr)
library(glue)
library(reshape2)
library(ggplot2)

library(Seurat)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)

library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dendsort)
library(grid)

# ========== ArchR Setup (if needed) ==========
# addArchRThreads(threads = 16)
# addArchRGenome("mm10")

# ========== Paths ==========

dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/4_241113_TFenrichment_dynamic"
setwd(dir)

# ========== Input Files ==========

wnn.seu.VTA <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/6_231107_updated_analyses/2_ATAC_alternative_processing/3_231205_ArchR_DARsandMotifs/231219_rpkm_based_DARs/240119_thresholds_setC/240122.wnn.seu.VTA.withMotifs.rds")
wnn.seu.VTA.1291 <- subset(wnn.seu.VTA, subset = orig.ident %notin% c("m30_cocaine_R1", "h1_saline_R1"))

motifs_for_heatmap <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/2_241113_TFenrichment_contrast/241121_motifs_for_heatmap.rds")
VTA_RPKM_formatted_metadata_withDARs <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/241121_VTA_RPKM_formatted_metadata_withDARs.rds")

# Add motif information to seurat object
pfm <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/1_ATAC_masterTable_continued/231217.pfm.rds")
wnn.seu.VTA.1291 <- AddMotifs(
  object = wnn.seu.VTA.1291,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
# ========== Analysis Parameters ==========

pval_threshold <- 0.05
logfc_threshold <- 0.25
logcpm_threshold <- 3

all_contrasts <- c(
  "h1_cocaine-saline","h4_cocaine-saline","h8_cocaine-saline",
  "h24_cocaine-saline","d14_cocaine-saline",
  "h4_cocaine-h1_cocaine","h8_cocaine-h4_cocaine",
  "h24_cocaine-h8_cocaine","d14_cocaine-h24_cocaine"
)


# ========== 1. Prepare Background Peak Sets ==========

#All backgrounds: In this iteration, fitlered peaks  == all peaks
bdg_allPeaks <- VTA_RPKM_formatted_metadata_withDARs %>% dplyr::select(peakID) %>%  mutate(peakID = str_replace(peakID, ":", "-")) %>% distinct() %>% pull()
bdg_filteredPeaks <- VTA_RPKM_formatted_metadata_withDARs %>% dplyr::filter(complete.cases(.)) %>% dplyr::select(peakID_v2) %>% distinct() %>% pull()
bdg_allDARs <- VTA_RPKM_formatted_metadata_withDARs %>% dplyr::filter(complete.cases(.), p_val < pval, abs(logFC) >=logfc, logCPM >= logcpm) %>% dplyr::select(peakID_v2) %>% distinct() %>% pull()

DAR_complete_results_kmeans <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/3_241113_kmeans/241114_DAR_complete_results_kmeans.tsv")
memory_TFenrichment.complete.df <- data.frame()
DARs_simp <- DAR_complete_results_kmeans %>% 
  dplyr::mutate(category_name = glue("{memory_class}_{direction}")) %>% 
  dplyr::select(category_name, peakID) 

# ========== 2. TF Enrichment per Dynamic Class ==========
unique_categories <- unique(DARs_simp$category_name)
for (category in unique_categories) {
  message(glue("Analyzing: {category}"))
  
  #A. Extract DARs: 
  selected_dars <- DARs_simp %>% dplyr::filter(category_name == category) %>% dplyr::select(peakID) %>% dplyr::mutate(peakID = str_replace(peakID, ":", "-")) %>% pull()
  
  #B.1. TFenrichment; bdg: allPeaks
  selected_dars_allPeaks.tfenrich <- FindMotifs(
    object = wnn.seu.VTA.1291,
    features = selected_dars, 
    background = bdg_allPeaks) %>% 
    dplyr::mutate(category_name = category, bdg = "allPeaks")
  
  #B.2. TFenrichment; bdg: filteredPeaks
  selected_dars_filteredPeaks.tfenrich <- FindMotifs(
    object = wnn.seu.VTA.1291,
    features = selected_dars, 
    background = bdg_filteredPeaks) %>% 
    dplyr::mutate(category_name = category, bdg = "filteredPeaks")
  
  #B.3. TFenrichment; bdg: allDARs
  selected_dars_allDARs.tfenrich <- FindMotifs(
    object = wnn.seu.VTA.1291,
    features = selected_dars, 
    background = bdg_allDARs) %>% 
    dplyr::mutate(category_name = category, bdg = "allDARs")
  
  memory_TFenrichment.complete.df <-rbind(memory_TFenrichment.complete.df,
                                          selected_dars_allPeaks.tfenrich,
                                          selected_dars_filteredPeaks.tfenrich,
                                          selected_dars_allDARs.tfenrich)
}
#write_tsv(memory_TFenrichment.complete.df, "241121.TFenrichment.complete.df.original_table.tsv")

#Filter for bdg = allDARs and "expressed" TFs
memory_TFenrichment.complete.df <- memory_TFenrichment.complete.df %>%
  dplyr::filter(bdg == "allDARs", motif.name %in% motifs_for_heatmap)  %>% #motifs_for_heatmap = 134 in this run (enriched and expressed)
  dplyr::mutate(mlog10_pval = -log10(pvalue))
write_tsv(memory_TFenrichment.complete.df, glue("241121_memory_TFenrichment.complete.allDARs.tsv"))

upreg_incocaine <- memory_TFenrichment.complete.df %>% filter(category_name %in% c("transient_up", "recovered_up", "memory_up", "delayed_up"), pvalue < 0.05) %>% distinct(motif.name) %>% pull()
upreg_insaline <- memory_TFenrichment.complete.df %>% filter(category_name %in% c("transient_down", "recovered_down", "memory_down", "delayed_down"), pvalue < 0.05) %>% distinct(motif.name) %>% pull()


# ========== 3. Heatmap Construction ==========
#A] TF enrichment (mlogpval)
coc_vs_sal_heatmap_TFs.df <- memory_TFenrichment.complete.df %>% 
  dplyr::filter(category_name %in% c("transient_up", "recovered_up", "memory_up", "delayed_up"), 
                motif.name %in% upreg_incocaine) %>%
  dplyr::select(motif.name,category_name, mlog10_pval)

#Cast dataframe
coc_vs_sal_heatmap_TFs.cast <- dcast(coc_vs_sal_heatmap_TFs.df, motif.name ~ category_name) %>% column_to_rownames("motif.name") #samples X TFs

#Convert to matrix:
coc_vs_sal_heatmap_TFs.matrix <- as.matrix(coc_vs_sal_heatmap_TFs.cast)

dend = as.dendrogram(hclust(dist(coc_vs_sal_heatmap_TFs.matrix)))
dend = color_branches(dend, k = 30)
col_runif = colorRamp2(c(0, 6), c("white", "red"))
coc_vs_sal_heatmap_TFs.matrix_trans <- t(coc_vs_sal_heatmap_TFs.matrix)[c("transient_up", "recovered_up", "memory_up", "delayed_up"), ]
coc_vs_sal.TFenrichment <- Heatmap(coc_vs_sal_heatmap_TFs.matrix_trans, name = "-log10(pvalue)", cluster_rows = FALSE, cluster_columns = dend, column_split = 10, col =col_runif,
                                   row_gap = unit(2, "mm"), column_gap = unit(2, "mm"), column_names_gp = gpar(fontsize = 9), 
                                   row_names_gp = gpar(fontsize = 10))

coc_vs_sal.TFenrichment


#B] corrected zscores
coc_vs_sal_heatmap_TFs.matrix_trans_zscores <- scale(coc_vs_sal_heatmap_TFs.matrix_trans)

dend = as.dendrogram(hclust(dist(coc_vs_sal_heatmap_TFs.matrix_trans)))
dend = color_branches(dend, k = 3)

coc_vs_sal.Zscores<- Heatmap(coc_vs_sal_heatmap_TFs.matrix_trans_zscores, name = "Z-score", cluster_rows = FALSE, cluster_columns = TRUE, column_split = 15,
                             row_gap = unit(2, "mm"), column_gap = unit(2, "mm"), column_names_gp = gpar(fontsize = 9),
                             row_names_gp = gpar(fontsize = 10))
coc_vs_sal.Zscores




#3.2. Saline vs cocaine: 
#########################
#A] TF enrichment (mlogpval)
sal_vs_coc_heatmap_TFs.df <- memory_TFenrichment.complete.df %>% 
  dplyr::filter(category_name %in% c("transient_down", "recovered_down", "memory_down", "delayed_down"), 
                motif.name %in% upreg_insaline) %>%
  dplyr::select(motif.name,category_name, mlog10_pval)

#Cast dataframe
sal_vs_coc_heatmap_TFs.cast <- dcast(sal_vs_coc_heatmap_TFs.df, motif.name ~ category_name) %>% column_to_rownames("motif.name") #samples X TFs

#Convert to matrix:
sal_vs_coc_heatmap_TFs.matrix <- as.matrix(sal_vs_coc_heatmap_TFs.cast)

dend = as.dendrogram(hclust(dist(sal_vs_coc_heatmap_TFs.matrix)))
dend = color_branches(dend, k = 30)
col_runif = colorRamp2(c(0, 6), c("white", "red"))
sal_vs_coc_heatmap_TFs.matrix_trans <- t(sal_vs_coc_heatmap_TFs.matrix)[c("transient_down", "recovered_down", "memory_down", "delayed_down"), ]
sal_vs_coc.TFenrichment <- Heatmap(sal_vs_coc_heatmap_TFs.matrix_trans, name = "-log10(pvalue)", cluster_rows = FALSE, cluster_columns = dend, column_split = 15, col =col_runif,
                                   row_gap = unit(2, "mm"), column_gap = unit(2, "mm"), column_names_gp = gpar(fontsize = 9), 
                                   row_names_gp = gpar(fontsize = 10))

sal_vs_coc.TFenrichment

#B] corrected zscores: note zscores are calculated from -log10pval
sal_vs_coc_heatmap_TFs.matrix_trans_zscores <- scale(sal_vs_coc_heatmap_TFs.matrix_trans)

dend = as.dendrogram(hclust(dist(sal_vs_coc_heatmap_TFs.matrix_trans)))
dend = color_branches(dend, k = 3)

sal_vs_coc_Zscores<- Heatmap(sal_vs_coc_heatmap_TFs.matrix_trans_zscores, name = "Z-score", cluster_rows = FALSE, cluster_columns = TRUE, column_split = 15,
                             row_gap = unit(2, "mm"), column_gap = unit(2, "mm"), column_names_gp = gpar(fontsize = 9),
                             row_names_gp = gpar(fontsize = 10))
sal_vs_coc_Zscores
#save.image("241121_TFenrichment_dynamic_class_final_complete_analysis.rds")


# ========== 4. Cluster reordering for easy visualization ==========
#4.1. Cocaine vs saline
##########################
coc_vs_sal_heatmap_TFs.df <- memory_TFenrichment.complete.df %>% 
  dplyr::filter(category_name %in% c("transient_up", "recovered_up", "memory_up", "delayed_up"), 
                motif.name %in% upreg_incocaine) %>%
  dplyr::select(motif.name,category_name, mlog10_pval)

#Cast dataframe
coc_vs_sal_heatmap_TFs.cast <- dcast(coc_vs_sal_heatmap_TFs.df, motif.name ~ category_name) %>% column_to_rownames("motif.name") #samples X TFs

#Convert to matrix:
coc_vs_sal_heatmap_TFs.matrix <- as.matrix(coc_vs_sal_heatmap_TFs.cast)

dend = as.dendrogram(hclust(dist(coc_vs_sal_heatmap_TFs.matrix)))
dend = color_branches(dend, k = 30)
col_runif = colorRamp2(c(0, 8), c("white", "red"))
coc_vs_sal_heatmap_TFs.matrix_trans <- t(coc_vs_sal_heatmap_TFs.matrix)[c("transient_up", "recovered_up", "memory_up", "delayed_up"), ]
coc_vs_sal.TFenrichment <- Heatmap(coc_vs_sal_heatmap_TFs.matrix_trans, name = "-log10(pvalue)", cluster_rows = FALSE, cluster_columns = dend, column_split = 10, col =col_runif,
                                   row_gap = unit(2, "mm"), column_gap = unit(2, "mm"), column_names_gp = gpar(fontsize = 9), 
                                   row_names_gp = gpar(fontsize = 10))

coc_vs_sal.TFenrichment

#a) get the vectors for reordering and split
col_order <- column_order(coc_vs_sal.TFenrichment)
new_cluster_order <- c(5, 2, 10, 6, 7, 9, 4, 8, 3, 1)  # Desired cluster order
new_col_order <- col_order[new_cluster_order]  # Reorder col_order

# Flatten the list into a single vector
flattened_order <- unlist(new_col_order)

# Create a factor for column splits (one group per sublist in new_col_order)
column_groups_cat <- factor(rep(seq_along(new_col_order), lengths(new_col_order)), levels = 1:10)

# Reorder the matrix columns
reordered_matrix <- coc_vs_sal_heatmap_TFs.matrix_trans[, flattened_order]

# Create the heatmap
col_runif = colorRamp2(c(0, 8), c("white", "red"))

coc_vs_sal_reordered_HM <- Heatmap(
  reordered_matrix,
  column_split = column_groups_cat,   # Split into groups
  name = "-log10(pvalue)",
  cluster_rows = FALSE,               # Rows are not clustered
  cluster_columns = FALSE,            # Turn off column clustering
  col = col_runif,                    # Your specified color palette
  row_gap = unit(2, "mm"),            # Gap between row groups
  column_gap = unit(2, "mm"),         # Gap between column groups
  column_names_gp = gpar(fontsize = 9),  # Font size for column names
  row_names_gp = gpar(fontsize = 10)
)
coc_vs_sal_reordered_HM

#In vertical (24/12):
# Transpose the reordered matrix
reordered_matrix_vertical <- t(reordered_matrix)

# Create the heatmap in vertical form
coc_vs_sal_reordered_HM_vertical <- Heatmap(
  reordered_matrix_vertical,
  row_split = column_groups_cat,       # Split rows into groups (previously columns)
  name = "-log10(pvalue)",
  cluster_rows = FALSE,                # Turn off row clustering (was column clustering before)
  cluster_columns = FALSE,             # Columns (TFs) are not clustered
  col = col_runif,                     # Your specified color palette
  row_gap = unit(2, "mm"),             # Gap between row groups
  column_gap = unit(2, "mm"),          # Gap between column groups
  row_names_gp = gpar(fontsize = 9),   # Font size for row names (was column names)
  column_names_gp = gpar(fontsize = 10) # Font size for column names (was row names)
)

# Display the vertical heatmap
coc_vs_sal_reordered_HM_vertical




#4.2. Saline vs cocaine
##########################
sal_vs_coc_heatmap_TFs.df <- memory_TFenrichment.complete.df %>% 
  dplyr::filter(category_name %in% c("transient_down", "recovered_down", "memory_down", "delayed_down"), 
                motif.name %in% upreg_insaline) %>%
  dplyr::select(motif.name,category_name, mlog10_pval)

#Cast dataframe
sal_vs_coc_heatmap_TFs.cast <- dcast(sal_vs_coc_heatmap_TFs.df, motif.name ~ category_name) %>% column_to_rownames("motif.name") #samples X TFs

#Convert to matrix:
sal_vs_coc_heatmap_TFs.matrix <- as.matrix(sal_vs_coc_heatmap_TFs.cast)

dend = as.dendrogram(hclust(dist(sal_vs_coc_heatmap_TFs.matrix)))
dend = color_branches(dend, k = 30)
col_runif = colorRamp2(c(0, 8), c("white", "red"))
sal_vs_coc_heatmap_TFs.matrix_trans <- t(sal_vs_coc_heatmap_TFs.matrix)[c("transient_down", "recovered_down", "memory_down", "delayed_down"), ]
sal_vs_coc.TFenrichment <- Heatmap(sal_vs_coc_heatmap_TFs.matrix_trans, name = "-log10(pvalue)", cluster_rows = FALSE, cluster_columns = dend, column_split = 15, col =col_runif,
                                   row_gap = unit(2, "mm"), column_gap = unit(2, "mm"), column_names_gp = gpar(fontsize = 9), 
                                   row_names_gp = gpar(fontsize = 10))

sal_vs_coc.TFenrichment


#a) get the vectors for reordering and split
col_order <- column_order(sal_vs_coc.TFenrichment)
new_cluster_order <- c(5,7,8,6,15,2,3,4,13,14,9,10,11,12,1)  # Desired cluster order

new_col_order <- col_order[new_cluster_order]  # Reorder col_order

# Flatten the list into a single vector
flattened_order <- unlist(new_col_order)

# Create a factor for column splits (one group per sublist in new_col_order)
column_groups_cat <- factor(rep(seq_along(new_col_order), lengths(new_col_order)), levels = 1:15)

# Reorder the matrix columns
reordered_matrix <- sal_vs_coc_heatmap_TFs.matrix_trans[, flattened_order]

# Create the heatmap
sal_vs_coc_reordered_HM <- Heatmap(
  reordered_matrix,
  column_split = column_groups_cat,   # Split into groups
  name = "-log10(pvalue)",
  cluster_rows = FALSE,               # Rows are not clustered
  cluster_columns = FALSE,            # Turn off column clustering
  col = col_runif,                    # Your specified color palette
  row_gap = unit(2, "mm"),            # Gap between row groups
  column_gap = unit(2, "mm"),         # Gap between column groups
  column_names_gp = gpar(fontsize = 9),  # Font size for column names
  row_names_gp = gpar(fontsize = 10)
)
sal_vs_coc_reordered_HM

#in vertical:
# Transpose the reordered matrix
reordered_matrix_vertical <- t(reordered_matrix)

# Create the heatmap in vertical form
sal_vs_coc_reordered_HM_vertical <- Heatmap(
  reordered_matrix_vertical,
  row_split = column_groups_cat,       # Split rows into groups (previously columns)
  name = "-log10(pvalue)",
  cluster_rows = FALSE,                # Turn off row clustering (was column clustering before)
  cluster_columns = FALSE,             # Columns (TFs) are not clustered
  col = col_runif,                     # Your specified color palette
  row_gap = unit(2, "mm"),             # Gap between row groups
  column_gap = unit(2, "mm"),          # Gap between column groups
  row_names_gp = gpar(fontsize = 9),   # Font size for row names (was column names)
  column_names_gp = gpar(fontsize = 10) # Font size for column names (was row names)
)

# Display the vertical heatmap
sal_vs_coc_reordered_HM_vertical

#save.image("241203_TFenrichment_dynamic_class_final_complete_analysis.rds")

# ========== 5. Exploration of motifs  ==========
#Which dynamic classes have significant AP1_var1 binding?
AP1_var1 <- c("FOSL2::JUN", "FOSL2::JUND", "FOSL2::JUNB", "FOS::JUN", "FOSL1::JUN", "FOS::JUND", "FOSL1::JUND", "FOSB::JUNB", "FOSL1::JUNB", "FOS::JUNB", "BATF::JUN", "JUN::JUNB")
AP1_var2 <- c("FOSB::JUNB(var.2)", "FOSL1::JUN(var.2)") 

TFenrich_ap1_var1 <- memory_TFenrichment.complete.df %>% 
  dplyr::filter(motif.name %in% AP1_var1) %>% 
  dplyr::mutate(is_significant = ifelse(pvalue <= 0.05, "yes", "no")) %>% 
  dplyr::filter(is_significant == "yes") %>% 
  dplyr::select(category_name, bdg, mlog10_pval, motif.name)

TFenrich_ap1_var2 <- memory_TFenrichment.complete.df %>% 
  dplyr::filter(motif.name %in% AP1_var2) %>% 
  dplyr::mutate(is_significant = ifelse(pvalue <= 0.05, "yes", "no")) %>% 
  dplyr::filter(is_significant == "yes") %>% 
  dplyr::select(category_name, bdg, mlog10_pval, motif.name)

# ========== Outputs ==========
# - Filtered TF enrichment tables
# - Heatmap objects and figures
# - Cluster ordering tables

