# ===========================================
# Script Title: 3.5 VTA-DN TF filtering by RNA expression and promoter accessibility
# Author: Luna Zea Redondo
# Date: 2026-06-09
#
# Description:
#   This script filters preselected TF motifs based on RNA expression
#   and promoter ATAC accessibility across conditions.
#   TFs are retained if the corresponding gene passes the RNA or ATAC
#   signal threshold in at least one treatment.
# ===========================================

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() 
#.rs.restartR()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))

# ========== Load Required Libraries ==========
library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(glue)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(ggrepel)
library(ggExtra)

dir <- "/fast/AG_Pombo/luna/2026_rebuttal/13_TFenrichment-memory"
setwd(dir)

# =========================
# Parameters
# =========================

conditions_use <- c("saline", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")

threshold_quantile <- 0.25

#and data:
gene_info_file <- "/fast/AG_Pombo/luna/2026_rebuttal/gene_info.txt"
promoters_file <- "/fast/AG_Pombo/luna/2026_rebuttal/promoters.bed"

TFenrichment.complete.df <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/10_DARs_TFenrichment-contrast/260610.TFenrichment.complete.df.tsv")
preselected_TFs_unique_vec <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/10_DARs_TFenrichment-contrast/280610_preselected_TFs_unique_vec_allcontrast.rds")


mouse_human_conversion <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/human_mouse_conversion.txt") %>% 
  dplyr::select(human_gene_name, mouse_gene_name)

HUGO <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/HUGO_aliases.txt") %>% 
  dplyr::select(`Approved symbol` , `Previous symbols` , `Alias symbols`)  
names(HUGO) <- c("human_gene_name", "human_previous_symbols", "human_gene_aliases")

Jaspar_mouse_equivalence <- split_names %>%
  mutate(equivalent = case_when(
    motif.name %in% mouse_human_conversion$human_gene_name ~ mouse_human_conversion$mouse_gene_name[match(motif.name, mouse_human_conversion$human_gene_name)],
    motif.name %in% mouse_human_conversion$mouse_gene_name ~ motif.name,
    TRUE ~ NA_character_),
    pre_selected = ifelse(motif.name %in% preselected_TFs_unique_vec, "Yes", "No")
  )

Jaspar_mouse_curated <- Jaspar_mouse_equivalence %>% 
  dplyr::mutate(gene_symbol = case_when(
    motif.name == "EWSR1-FLI1" ~ "Ewsr1, Fli1",
    motif.name == "Znf281" ~ "Zfp281", 
    motif.name == "PITX2" ~ "Pitx2", 
    motif.name == "Znf423" ~ "Zfp423", 
    motif.name == "ZNF135" ~ "Zfp61", 
    motif.name == "YY2" ~ "Yy2", 
    TRUE ~ equivalent
  ))

# =========================
# Helper function: rpkm
# =========================

calculate_rpkm <- function(counts, length_bp) {
  total_reads <- colSums(counts)
  sweep(sweep(counts, 1, length_bp / 1000, "/"), 2, total_reads / 1e6, "/")
}

# =========================
# 1. Motif names -> TF gene symbols
# =========================

split_names <- TFenrichment.complete.df %>% 
  dplyr::mutate(motif.name = ifelse(motif.name == "EWSR1-FLI1", "EWSR1::FLI1", motif.name)) %>% 
  dplyr::select(motif.name) %>% 
  tidyr::separate_rows(motif.name, sep = "::|\t") %>%
  dplyr::filter(motif.name != "") %>%
  dplyr::mutate(motif.name = gsub("\\(.*?\\)", "", motif.name)) %>% 
  distinct()

TF_genes <- unique(Jaspar_mouse_curated$gene_symbol)

TF_genes_preselected <- Jaspar_mouse_curated %>% 
  dplyr::filter(pre_selected == "Yes") %>% 
  distinct(gene_symbol) %>% 
  pull(gene_symbol)

# =========================
# 2. RNA RPKM per condition
# =========================

gene_lengths <- read_table(
  gene_info_file,
  col_names = c("gene_id", "gene_symbol", "length_bp")
)

DefaultAssay(seu.VTA_DNs) <- "RNA"

GEX_raw_counts <- AverageExpression(object = seu.VTA_DNs, assay = "RNA", slot = "counts", group.by = "simpleIdent")$RNA
colnames(GEX_raw_counts) <- gsub("-", "_", colnames(GEX_raw_counts))

GEX_raw_counts <- GEX_raw_counts[, conditions_use, drop = FALSE]

rna_lengths <- gene_lengths$length_bp[match(rownames(GEX_raw_counts), gene_lengths$gene_symbol)]
keep_genes <- !is.na(rna_lengths)

GEX_raw_counts <- GEX_raw_counts[keep_genes, , drop = FALSE]
rna_lengths <- rna_lengths[keep_genes]

GEX_rpkms_all <- calculate_rpkm(as.matrix(GEX_raw_counts), rna_lengths)

GEX_rpkms_all_withTFinfo <- as.data.frame(GEX_rpkms_all) %>%
  rownames_to_column("gene_symbol") %>%
  dplyr::mutate(TF_class = case_when(
    gene_symbol %in% TF_genes_preselected ~ "pre_selected",
    gene_symbol %in% TF_genes ~ "other_TF",
    TRUE ~ "no_TF"
  ))

GEX_rpkms_all_withTFinfo_tomerge <- GEX_rpkms_all_withTFinfo %>%
  pivot_longer(
    cols = all_of(conditions_use),
    names_to = "condition",
    values_to = "rna_value"
  )

# =========================
# 3. Promoter ATAC RPKM per condition
# =========================

DefaultAssay(atac) <- "ATAC"

promoters <- read_table(
  promoters_file,
  col_names = c("chr", "gene_symbol", "gene_id", "transcript_id", "start", "end", "strand")) %>%
  dplyr::select(chr, gene_symbol, start, end) %>%
  dplyr::filter(chr %in% c(paste0("chr", 1:19), "chrX")) %>%
  distinct() %>%
  dplyr::mutate(peakID = glue("{chr}-{start}-{end}"))

promoters.gr <- GRanges(
  seqnames = promoters$chr,
  ranges = IRanges(start = promoters$start, end = promoters$end)
)

names(promoters.gr) <- promoters$peakID

promoter_counts <- FeatureMatrix(
  fragments = Fragments(atac),
  features = promoters.gr,
  cells = colnames(atac)
)

promoter_counts_condition <- sapply(conditions_use, function(cond) {
  cells_use <- colnames(atac)[atac$condition == cond]
  Matrix::rowSums(promoter_counts[, cells_use, drop = FALSE])
})

rownames(promoter_counts_condition) <- rownames(promoter_counts)

promoter_lengths <- width(promoters.gr)[match(rownames(promoter_counts_condition), names(promoters.gr))]

promoter_rpkm <- calculate_rpkm(as.matrix(promoter_counts_condition), promoter_lengths)
rownames(promoter_rpkm) <- gsub("^([^.]*)\\.([0-9]+)\\.([0-9]+)$", "\\1:\\2-\\3", rownames(promoter_rpkm))

promoters_averaged_annotated.df <- as.data.frame(promoter_rpkm) %>%
  tibble::rownames_to_column("peakID") %>%
  dplyr::mutate(peakID = gsub("\\.", "-", peakID)) %>%
  dplyr::left_join(
    promoters %>% dplyr::select(peakID, gene_symbol),
    by = "peakID"
  )

most_accessible_promoters <- promoters_averaged_annotated.df %>%
  tidyr::pivot_longer(
    cols = all_of(conditions_use),
    names_to = "condition",
    values_to = "atac_value"
  ) %>%
  dplyr::filter(!is.na(gene_symbol)) %>%
  dplyr::group_by(gene_symbol, condition) %>%
  dplyr::slice_max(atac_value, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# =========================
# 4. Combine RNA and ATAC
# =========================

RNA_and_ATAC_toPlot <- left_join(
  GEX_rpkms_all_withTFinfo_tomerge,
  most_accessible_promoters,
  by = c("gene_symbol", "condition")
) %>%
  mutate(
    rna_value = replace_na(rna_value, 0),
    atac_value = replace_na(atac_value, 0)
  )

rna_threshold <- as.numeric(quantile(
  RNA_and_ATAC_toPlot$rna_value[RNA_and_ATAC_toPlot$rna_value > 0],
  probs = threshold_quantile,
  na.rm = TRUE
))

atac_threshold <- as.numeric(quantile(
  RNA_and_ATAC_toPlot$atac_value[RNA_and_ATAC_toPlot$atac_value > 0],
  probs = threshold_quantile,
  na.rm = TRUE
))

message("RNA threshold: ", round(rna_threshold, 4))
message("ATAC threshold: ", round(atac_threshold, 4))
# =========================
# 5. Extract preselected TF 
# =========================
#Genes  passing RNA or ATAC signal
genes_passing <- RNA_and_ATAC_toPlot %>%
  dplyr::filter(
    TF_class == "pre_selected",
    atac_value >= atac_threshold | rna_value >= rna_threshold
  ) %>%
  distinct(gene_symbol) %>%
  pull(gene_symbol)

genes_passing_original <- Jaspar_mouse_curated %>%
  dplyr::filter(gene_symbol %in% genes_passing) %>%
  distinct(motif.name) %>%
  pull(motif.name)

selected_TFs <- as.data.frame(preselected_TFs_unique_vec) %>%
  mutate(motifs = str_extract_all(preselected_TFs_unique_vec, "[^::\\(]+")) %>%
  dplyr::filter(sapply(motifs, function(x) any(x %in% genes_passing_original)))

motifs_for_heatmap <- unique(c(selected_TFs$preselected_TFs_unique_vec, "CTCFL"))
#saveRDS(motifs_for_heatmap, "260610_motifs_for_heatmap.rds")

# =========================
# 6. Check output
# =========================

length(genes_passing)
length(motifs_for_heatmap)

genes_passing
motifs_for_heatmap


# =========================
# 7. Visualization
# =========================

# Saline only: all genes / TF classes
RNA_ATAC_saline_only <- RNA_and_ATAC_toPlot %>%
  dplyr::filter(condition == "saline")

RNA_ATAC_saline_only.plot <- ggplot(
  RNA_ATAC_saline_only,
  aes(x = atac_value, y = rna_value, color = TF_class)
) +
  geom_point(size = 2.5, alpha = 0.3) +
  scale_color_manual(values = c("no_TF" = "gray90", "other_TF" = "royalblue3", "pre_selected" = "darkred")) +
  scale_x_log10(limits = c(0.003, NA)) +
  scale_y_log10(limits = c(NA, 100)) +
  geom_vline(xintercept = atac_threshold, linetype = "dashed", color = "goldenrod") +
  geom_hline(yintercept = rna_threshold, linetype = "dashed", color = "goldenrod") +
  # geom_text_repel(
  #   data = subset(RNA_ATAC_saline_only, TF_class == "pre_selected"),
  #   aes(label = gene_symbol),
  #   size = 3,
  #   color = "darkred",
  #   max.overlaps = 100
  # ) +
  labs(title = "Saline", x = "Promoter ATAC RPKM", y = "RNA RPKM") +
  theme_minimal() +
  theme(legend.position = "bottom")

RNA_ATAC_saline_only.plot.complete <- ggExtra::ggMarginal(
  RNA_ATAC_saline_only.plot,
  groupColour = TRUE,
  groupFill = FALSE
)

RNA_ATAC_saline_only.plot.complete


# All conditions: pre-selected TFs only
RNA_ATAC_ALL.plot <- RNA_and_ATAC_toPlot %>%
  dplyr::filter(TF_class == "pre_selected") %>%
  ggplot(aes(x = atac_value, y = rna_value)) +
  geom_point(size = 2.5, alpha = 0.2, color = "darkred") +
  scale_x_log10() +
  scale_y_log10() +
  geom_vline(xintercept = atac_threshold, linetype = "dashed", color = "goldenrod") +
  geom_hline(yintercept = rna_threshold, linetype = "dashed", color = "goldenrod") +
  geom_text_repel(
    aes(label = gene_symbol),
    size = 3,
    color = "darkred",
    max.overlaps = 100
  ) +
  facet_wrap(~ condition, nrow = 3) +
  labs(x = "Promoter ATAC RPKM", y = "RNA RPKM") +
  theme_minimal()

RNA_ATAC_ALL.plot



#save.image("260610_TFcleaning.rds")

