# ===========================================
# Script Title: 3.4 VTA-DN DARs and TF enrichment (partI)
# Author: Luna Zea Redondo
# Date: 2026-06-09

# Description:
#   This script explores DARs (all pairwise analyses)
#   and TF enrichment analysis per contrast
# ===========================================


# ========== Set Environment ==========
rm(list = ls())
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230310/site-library"))
set.seed(1)

# ========== Load Required Libraries ==========
library(ArchR)
library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(glue)
library(ComplexHeatmap)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(patchwork)
library(readr)
library(forcats)
library(GenomeInfoDb)
library(ggpubr)
library(PupillometryR)
library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library(Matrix)

source("/fast/AG_Pombo/luna/2026_rebuttal/config.R", local = FALSE)

# ========== Set Working Directory ==========
dir <- "/fast/AG_Pombo/luna/2026_rebuttal/10_DARs_TFenrichment-contrast"
setwd(dir)

# ========== ArchR Setup ==========
addArchRThreads(threads = 16)
addArchRGenome("mm10")


# ===========================================
# 1) Pre-process DARs and create master table
# ===========================================

# Load data:
atac <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/9_ATAC-QCandIntegration/260608_signac_atac_peak_matrix_object_annotated.rds")

#Calculate RPKM per peak per condition
conditions_use <- c("saline", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")

counts <- GetAssayData(atac, assay = "ATAC", layer = "counts")
peaks <- granges(atac[["ATAC"]])
names(peaks) <- rownames(counts)
peak_width_kb <- width(peaks) / 1000

missing_conditions <- setdiff(conditions_use, unique(atac$condition))
if (length(missing_conditions) > 0) stop("Missing conditions in atac$condition: ", paste(missing_conditions, collapse = ", "))

pseudo_bulk_counts <- sapply(conditions_use, function(cond) {
  cells_use <- colnames(atac)[atac$condition == cond]
  Matrix::rowSums(counts[, cells_use, drop = FALSE])
})

rownames(pseudo_bulk_counts) <- rownames(counts)

lib_size_million <- Matrix::colSums(pseudo_bulk_counts) / 1e6

rpkm_per_condition <- sweep(pseudo_bulk_counts, 1, peak_width_kb, "/")
rpkm_per_condition <- sweep(rpkm_per_condition, 2, lib_size_million, "/")

rpkm_per_condition_df <- as.data.frame(rpkm_per_condition) %>%
  tibble::rownames_to_column("peakID")

tail(rpkm_per_condition_df)

pval_threshold = 0.05
logfc_threshold = 0.5 
logcpm_threshold = 3.5 


#Load DARs dataframe:
DARs_corrected <- data.table::fread("260604.peaks-edgeR.tsv.gz") %>% 
  dplyr::filter(contrast %in% all_contrasts) %>% 
  dplyr::mutate(peakID_v2 = str_replace(gene_id, "(\\w+)-(\\d+)-(\\d+)", "\\1:\\2-\\3"), 
                diff = ifelse(abs(logFC) >= logfc_threshold & p_val < pval_threshold & logCPM >= logcpm_threshold, "Yes", "No")) %>% 
  dplyr::rename(peakID = "gene_id")

dplyr::n_distinct(DARs_corrected$peakID[DARs_corrected$diff == "Yes"])

# Create master table
# Peak metadata from ATAC assay
peak_metadata <- atac[["ATAC"]]@meta.features %>%
  tibble::rownames_to_column("peakID")

# DARs table: make sure peak column is called peakID
DARs_info <- DARs_corrected 

# RPKM table: add rpkm_ prefix to condition columns
rpkm_info <- rpkm_per_condition_df %>%
  dplyr::rename_with(
    ~ paste0("rpkm_", .x),
    .cols = -peakID
  )

# Merge everything
complete_peaks_table <- DARs_info %>%
  dplyr::left_join(peak_metadata, by = "peakID") %>%
  dplyr::left_join(rpkm_info, by = "peakID")

# Check
# dim(complete_peaks_table)
# head(complete_peaks_table)
# colnames(complete_peaks_table)

#Add two more metrics and save an unaltered copy of all VTA-DN peaks
# Use peakID_v2 if it contains coordinates like chr1:1000-1500
peak_coords <- complete_peaks_table$peakID_v2

# Convert chr:start-end to GRanges
gr <- GRanges(
  seqnames = str_extract(peak_coords, "^chr[^:]+"),
  ranges = IRanges(
    start = as.integer(str_extract(peak_coords, "(?<=:)\\d+")),
    end   = as.integer(str_extract(peak_coords, "(?<=-)\\d+"))
  )
)

# Peak length
complete_peaks_table$peak_length <- width(gr)

# Extract sequences
peak_seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr)

# Count CG dinucleotides
complete_peaks_table$CG_count <- vcountPattern("CG", peak_seqs)
complete_peaks_table$CG_percent <- 
  complete_peaks_table$CG_count * 2 / complete_peaks_table$peak_length * 100

complete_peaks_table <- complete_peaks_table %>%
  dplyr::rename(
    GC_count = CG_count,
    GC_percent = CG_percent
  )
#saveRDS(complete_peaks_table, "260610_VTA-DNs_complete_peakset.rds")

# ===========================================
# 2) Visualize DARs: MA plots per contrast
# ===========================================

condition_colors <- c("black", "#9F2365", "#617641", "#C48208", "#326186", "#AE430A", "#564686")
condition_names <- c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")
names(condition_colors) <- condition_names

MAplots <- list()
for (contrast in all_contrasts) {
  message(glue("Analyzing: {contrast}"))
  group1=str_split(contrast, "-")[[1]][1]
  group2=str_split(contrast, "-")[[1]][2]
  contrast_label=contrast
  
  volcano_colors <- c(condition_colors[group1], condition_colors[group2], "gray80")
  names(volcano_colors) <- c("Upregulated", "Downregulated", "No significant")
  
  DARS_to_plot <- DARs_corrected %>% #load
    dplyr::filter(contrast == contrast_label) %>% 
    dplyr::mutate(is_DAR = ifelse(p_val < pval_threshold & logFC >= logfc_threshold & logCPM >= logcpm_threshold, "Upregulated", 
                                  ifelse(p_val < pval_threshold & logFC <= -logfc_threshold & logCPM >= logcpm_threshold, "Downregulated", "No significant")))
  
  downreg_num <- table(DARS_to_plot$is_DAR)["Downregulated"] 
  upreg_num <- table(DARS_to_plot$is_DAR)["Upregulated"] 
  
  MA_plot <- ggplot(DARS_to_plot, aes(x = logCPM, y = logFC)) +
    geom_point(alpha=0.6, aes(color = is_DAR)) +
    geom_point(alpha=0.6,aes(color = is_DAR), data=DARS_to_plot[DARS_to_plot$is_DAR != "No significant",]) +
    scale_color_manual(values = volcano_colors) +
    theme_classic(base_size = 12) + theme(legend.position = "bottom") +
    labs(title = glue("{contrast_label}"),
         subtitle = glue("{upreg_num} up- and {downreg_num} down-regulated"), 
         x= "logCPM", 
         y="LogFC") +
    guides(fill=guide_legend(title="is DARs")) +
    geom_vline(xintercept = logcpm_threshold, linetype="dotted", color = "darkred", linewidth=1.5) + 
    theme(legend.position = "none")
  MA_plot
  #ggsave(filename = glue("{dirplots}/{contrast_label}_peakANDmotifs_results_MAplot.png"), MA_plot, units = "px", device = "png", width=6000,height=2000,dpi = 300)  
  MAplots[[paste("MAplot", contrast, sep = "_")]] <- MA_plot
}

MA_row1 <- grid.arrange(grobs = MAplots[1:5], nrow = 1)
MA_row2 <- grid.arrange(grobs = MAplots[6:9], nrow = 1)


blank <- ggplot() + theme_void() 

# Arrange the plots
combined_layout <- (MAplots[[1]] + MAplots[[2]] + MAplots[[3]] + MAplots[[4]] + MAplots[[5]] +
                      MAplots[[6]] + MAplots[[7]] + MAplots[[8]] + MAplots[[9]] + blank) + plot_layout(ncol = 5)
combined_layout
#Save
# width_original = 10.30
# height_original= 7.31
# 
# plot_name <- "chapter3_atac_DARs_MAplots.pdf"
# 
# pdf(plot_name, width = width_original, height =height_original)
# combined_layout
# dev.off()

# Save as TIFF
# width_original = 12
# height_original= 7.31
# plot_name <- "chapter3_atac_DARs_MAplots.tiff"
# tiff(filename = plot_name, width = width_original, height = height_original, units = "in", res = 300, compression = "lzw")
# combined_layout
# dev.off()



#Plot: rpkm, peak length and GC content panel for sup figure
# This QC is done in a different script: 260612_QCpanels.R



# ================================
# 3) TF enrichment per contrast
# ================================
#All backgrounds: In this iteration, fitlered peaks  == all peaks
bdg_allPeaks <- complete_peaks_table %>% dplyr::select(peakID) %>% distinct() %>% pull()
bdg_allDARs <- complete_peaks_table %>% dplyr::filter(complete.cases(.), diff == "Yes") %>%  dplyr::select(peakID) %>% distinct() %>% pull()

#Save DARs in table (for kmeans):
DARs.table <- complete_peaks_table %>% dplyr::filter(complete.cases(.), diff == "Yes") %>%
  dplyr::distinct(peakID)
#write_tsv(DARs.table, "260610_vta.edgeR.DARs.txt")

DARs.table.late <- complete_peaks_table %>% dplyr::filter(complete.cases(.), diff == "Yes", contrast %in% c("d14_cocaine-saline", "h24_cocaine-saline", "d14_cocaine-h24_cocaine")) %>%
  dplyr::distinct(peakID)
#write_tsv(DARs.table.late, "260610_vta.edgeR.onlyLATE_DARs.txt")

# Add motif information to seurat object
pfm <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/231217.pfm.rds")
atac <- AddMotifs(
  object = atac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

#saveRDS(atac, "/fast/AG_Pombo/luna/2026_rebuttal/11_DARs_TFenrichment-contrast/signac_atac_with_motifs.rds")

TFenrichment.complete.df <- data.frame()

#Extract enriched TFs
for (contrast_label in all_contrasts) {
  message(glue("Analyzing: {contrast_label}"))
  group1=str_split(contrast_label, "-")[[1]][1]
  group2=str_split(contrast_label, "-")[[1]][2]
  
  #A. Extract DARs and divide by contrast and direction
  contrast.group1.dars <- DARs_corrected %>% dplyr::filter(contrast == contrast_label, logFC >= logfc_threshold, p_val < pval_threshold, logCPM >= logcpm_threshold) %>% dplyr::select(peakID) %>% pull()
  contrast.group2.dars <- DARs_corrected %>% dplyr::filter(contrast == contrast_label, logFC <= -logfc_threshold, p_val < pval_threshold, logCPM >= logcpm_threshold) %>% dplyr::select(peakID) %>% pull()
  
  #B.1. TFenrichment; bdg: allPeaks
  contrast.group1vsallPeaks.tfenrich <- FindMotifs(
    object = atac,
    features = contrast.group1.dars, 
    background = bdg_allPeaks) %>% 
    dplyr::mutate(contrast = contrast_label, upreg_in = group1, bdg = "allPeaks")
  
  contrast.group2vsallPeaks.tfenrich <- FindMotifs(
    object = atac,
    features = contrast.group2.dars, 
    background = bdg_allPeaks) %>% 
    dplyr::mutate(contrast = contrast_label, upreg_in = group2, bdg = "allPeaks")

  
  #B.2. TFenrichment; bdg: allDARs
  contrast.group1vsallDARs.tfenrich <- FindMotifs(
    object = atac,
    features = contrast.group1.dars, 
    background = bdg_allDARs) %>% 
    dplyr::mutate(contrast = contrast_label, upreg_in = group1, bdg = "allDARs")
  
  contrast.group2vsallDARs.tfenrich <- FindMotifs(
    object = atac,
    features = contrast.group2.dars, 
    background = bdg_allDARs) %>% 
    dplyr::mutate(contrast = contrast_label, upreg_in = group2, bdg = "allDARs")
  
  TFenrichment.complete.df <- rbind(TFenrichment.complete.df,
                                    contrast.group1vsallPeaks.tfenrich, contrast.group2vsallPeaks.tfenrich,
                                    contrast.group1vsallDARs.tfenrich, contrast.group2vsallDARs.tfenrich)
}
#write_tsv(TFenrichment.complete.df, "260610.TFenrichment.complete.df.tsv")


# ================================
# 4) TF–DAR Matrix Construction
# ================================
#Candidate TFs to study (only comparisons vs saline):
preselected_TFs.df <- TFenrichment.complete.df %>%
  group_by(contrast, upreg_in, bdg) %>% 
  dplyr::select(pvalue, motif.name, contrast, upreg_in, bdg) %>% 
  dplyr::filter(pvalue < 0.05, bdg == "allDARs") %>% 
  dplyr::filter(contrast %in% all_contrasts[1:5]) 

preselected_TFs_unique_vec <- preselected_TFs.df %>%
  dplyr::ungroup() %>% 
  dplyr::select(motif.name) %>% distinct() %>% pull()
#saveRDS(preselected_TFs_unique_vec, "280610_preselected_TFs_unique_vec_allcontrast.rds")

#Extract motif matrix and convert motifID to motif names. 
motif_motifName_conversion <- TFenrichment.complete.df %>% dplyr::select(motif, motif.name)
motif_matrix <- as.matrix(atac[["ATAC"]]@motifs@data)

motif_names <- motif_motifName_conversion$motif.name
colnames(motif_matrix) <- motif_names[match(colnames(motif_matrix), motif_motifName_conversion$motif)]

#Reduce matrix to candidate TFs:
motif_matrix_DARs <- motif_matrix[bdg_allDARs ,] #4983 X 746


# Important: save the DARs_corrected file, will come back to it often. 
# load("260610_DARs_and_TFenrichment_contrast.rds")
# saveRDS(DARs_corrected, "260610.DARs_corrected.rds")


#Plots (heatmaps)
TFenrichment.toPlot.df <- TFenrichment.complete.df %>%
  filter(bdg == "allDARs") %>%
  mutate(mlog10_pval = -log10(pvalue))

saline_contrasts <- c("h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline", "h24_cocaine-saline", "d14_cocaine-saline")
cocaine_conditions <- c("h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")
mlog10_threshold <- 1.75
col_fun <- colorRamp2(c(0, 10), c("white", "red"))

# Cocaine-up
TFs_coc <- TFenrichment.toPlot.df %>%
  filter(contrast %in% saline_contrasts, upreg_in != "saline") %>%
  group_by(motif.name) %>%
  summarise(max_mlog10 = max(mlog10_pval), .groups = "drop") %>%
  filter(max_mlog10 >= mlog10_threshold) %>%
  pull(motif.name)

mat_coc <- TFenrichment.toPlot.df %>%
  filter(contrast %in% saline_contrasts, upreg_in != "saline", motif.name %in% TFs_coc) %>%
  select(motif.name, upreg_in, mlog10_pval) %>%
  as.data.table() %>%
  dcast(motif.name ~ upreg_in, value.var = "mlog10_pval", fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("motif.name") %>%
  as.matrix()

mat_coc <- mat_coc[, cocaine_conditions, drop = FALSE]
mat_coc_t <- t(mat_coc)

dend_coc <- as.dendrogram(hclust(dist(mat_coc)))
dend_coc <- color_branches(dend_coc, k = min(10, nrow(mat_coc)))


HM_coc_up <- Heatmap(mat_coc_t, name = "-log10(pvalue)", cluster_rows = FALSE, cluster_columns = dend_coc,
                     column_split = min(10, nrow(mat_coc)), col = col_fun,
                     row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
                     column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 10))

HM_coc_up


# Saline-up
TFs_sal <- TFenrichment.toPlot.df %>%
  filter(contrast %in% saline_contrasts, upreg_in == "saline") %>%
  group_by(motif.name) %>%
  summarise(max_mlog10 = max(mlog10_pval), .groups = "drop") %>%
  filter(max_mlog10 >= mlog10_threshold) %>%
  pull(motif.name)

mat_sal <- TFenrichment.toPlot.df %>%
  filter(contrast %in% saline_contrasts, upreg_in == "saline", motif.name %in% TFs_sal) %>%
  select(motif.name, contrast, mlog10_pval) %>%
  as.data.table() %>%
  dcast(motif.name ~ contrast, value.var = "mlog10_pval", fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("motif.name") %>%
  as.matrix()

mat_sal <- mat_sal[, saline_contrasts, drop = FALSE]
mat_sal_t <- t(mat_sal)

dend_sal <- as.dendrogram(hclust(dist(mat_sal)))
dend_sal <- color_branches(dend_sal, k = min(10, nrow(mat_sal)))

HM_sal_up <- Heatmap(mat_sal_t, name = "-log10(pvalue)", cluster_rows = FALSE, cluster_columns = dend_sal,
                     column_split = min(10, nrow(mat_sal)), col = col_fun,
                     row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
                     column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 10))

HM_sal_up


condition_order <- c("h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")
contrast_labels <- c("h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline", "h24_cocaine-saline", "d14_cocaine-saline")

# Cocaine-up vertical
mat_coc_vert <- mat_coc[, condition_order, drop = FALSE]

HM_coc_up_vertical <- Heatmap(mat_coc_vert, name = "-log10(pvalue)", cluster_rows = dend_coc, cluster_columns = FALSE,
                              row_split = min(10, nrow(mat_coc)), col = col_fun,
                              row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
                              row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 10))

HM_coc_up_vertical

#Save
# dev.size()
# width_original = 8.17
# height_original= 9
# 
# plot_name <- "chapter3_atac_DARs_heatmap_subpanel.pdf"
# 
# pdf(plot_name, width = width_original, height =height_original)
# HM_coc_up_vertical
# dev.off()



# Saline-up vertical
mat_sal_vert <- mat_sal[, contrast_labels, drop = FALSE]
colnames(mat_sal_vert) <- condition_order

HM_sal_up_vertical <- Heatmap(mat_sal_vert, name = "-log10(pvalue)", cluster_rows = dend_sal, cluster_columns = FALSE,
                              row_split = min(10, nrow(mat_sal)), col = col_fun,
                              row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
                              row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 10))

HM_sal_up_vertical
