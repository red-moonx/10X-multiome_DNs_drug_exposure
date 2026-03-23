# ===========================================
# Script Title: 1.6_gene-category_classifiers.R
# Author: Luna Zea Redondo
# Date: 2024-11-04
# Description:
#   Integrates cocaine DEGs with CAGs and several gene-category classifiers
#   (psychiatric, addiction, neuropeptide, neurodegenerative sets).
#   Main steps:
#     - load DEG tables and external classifier gene lists
#     - annotate DEGs with gene-category flags
#     - generate MA plots across contrasts
#     - build heatmaps for cocaine vs saline expression
#     - summarize classifier proportions by direction and timepoint
#     - create UpSet plots for classifier overlap
#     - run Fisher / chi-square tests for enrichment
#     - compare intersections between CAGs and psychiatric genes
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
library(stringr)
library(tibble)
library(glue)

library(ggplot2)
library(ggrepel)
library(patchwork)

library(Seurat)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)

library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(RColorBrewer)

library(PupillometryR)
library(scales)
library(grid)

# ========== Project setup ==========
source("/fast/AG_Pombo/luna/2023_vta_multiome/2024/5_figures/pixels_to_cm.R")
setwd("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/9_241104_CAGs_psychiatric_association")
load("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002.RNAonly.excluding_putative_SN_cells.rds")

# ========== 1. Load DEG tables and classifier gene sets ==========
DEG_complete_results <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002_muscat_VTAvsSN_all_subtypes/Analysis2/deg_results/231002_DEG_complete_results_Analysis2.tsv") %>% 
  dplyr::filter(cluster_id == "VTA")
genes <- DEG_complete_results %>% dplyr::select(gene) %>% distinct()

#Complete DEG info for label selection
#DEGs_complete_info <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/240606_DEG_kmeans_info_coords_GAM_final_table.tsv")
DEGs_complete_info <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/250205_DEG_complete_results_kmeans.corrected.4243.rds")

# Updated Gene lists
####################
#Add mouse human conversion
mouse_human_conversion <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/5_240219_dynamic_cluster_classifiers/human_mouse_conversion.txt") %>% 
  dplyr::select(human_gene_name, mouse_gene_name) 
mouse_human_conversion$human_gene_name[mouse_human_conversion$mouse_gene_name == "Hnrnpa1"] <- "HNRNPA1"

# ========== 2. Reorder and save DEG annotation table ==========
#1. Dominik list of CAGs:
cags_all <- read_tsv("dominik_CAGs_harmonized_n584.tsv") %>% dplyr::select(gene_name) %>% dplyr::pull(gene_name)

#2. Psychiatric disorders:
psychiatric_all <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/5_figures/240812_DEGs_MAandHeatmap/240826_psychotic_disorders_final.tsv") %>%
  distinct() %>%
  mutate(PsychiatricDisorder_simplified = case_when(
    PsychiatricDisorder == "Alcohol use disorders" ~ "AD",
    PsychiatricDisorder == "Cocaine use disorders" ~ "CD",
    PsychiatricDisorder == "Schizophrenia spectrum and other psychotic disorders" ~ "SCH",
    PsychiatricDisorder == "Depressive disorders" ~ "DD",
    PsychiatricDisorder == "Substance/drug induced depressive disorder" ~ "SIDD",
    PsychiatricDisorder == "Bipolar disorders and related disorders" ~ "BD",
    PsychiatricDisorder == "Cannabis use disorders" ~ "CAD",
    PsychiatricDisorder == "Drug-induced psychosis" ~ "DIP",
    TRUE ~ "No"  # Default case for unmatched disorders
  )) %>% 
  dplyr::select(mouse_gene_name, PsychiatricDisorder_simplified) %>% 
  group_by(mouse_gene_name) %>% 
  dplyr::summarise(PsychiatricDisorder_collapsed = paste(unique(PsychiatricDisorder_simplified), collapse = "|"),
                   .groups = "drop") %>% dplyr::pull(mouse_gene_name) 

#3. Other addiction:
iza_anco_updated <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/5_figures/240812_DEGs_MAandHeatmap/240823_gene_categories/ANCO_genes_curated_to_mouse_genes.txt") %>% 
  dplyr::select(gene_name, addiction_phenotype) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(addiction_category = if_else(
    addiction_phenotype == "not", "not", if_else(
      addiction_phenotype == "CD", "only_cocaine", "other")))
general_addiction_all <- iza_anco_updated %>% dplyr::filter(addiction_category == "other") %>% dplyr::pull(gene_name)

#4. Alzheimer:
alzheimer_all <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/5_figures/240812_DEGs_MAandHeatmap/240823_gene_categories/alzheimer_genes.tsv") %>% 
  left_join(mouse_human_conversion, by = c("human_gene" = "human_gene_name")) %>% 
  dplyr:: filter(!is.na(mouse_gene_name))
#5. Parkinson:
parkinson_all <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/5_figures/240812_DEGs_MAandHeatmap/240823_gene_categories/parkinson_genes.tsv") %>% 
  left_join(mouse_human_conversion, by = c("human_gene" = "human_gene_name")) %>% 
  dplyr:: filter(!is.na(mouse_gene_name))
#6. ALS
als_all <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/5_figures/240812_DEGs_MAandHeatmap/240823_gene_categories/als_genes.tsv") %>% 
  left_join(mouse_human_conversion, by = c("human_gene" = "human_gene_name")) %>% 
  dplyr:: filter(!is.na(mouse_gene_name))

#7. Neuropeptides
neuropeptides <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/5_figures/240812_DEGs_MAandHeatmap/neuropeptides.tsv") %>% 
  left_join(mouse_human_conversion, by = c("Approved_symbol" = "human_gene_name")) %>% 
  dplyr:: filter(!is.na(mouse_gene_name))


#Reorder and save the DEGs complete info dataframe
#####################################################
DEGs_complete_info_final <- DEGs_complete_info[, c(1,2,3,4,5,6,7,8,14,15,16,17,18,19)] %>%
  dplyr::mutate(
    is_CAG = ifelse(gene %in% cags_all, "yes", "no"), 
    is_addiction = ifelse(gene %in% general_addiction_all, "yes", "no"), 
    is_psychiatric = ifelse(gene %in% psychiatric_all, "yes", "no"), 
    is_neuropeptide = ifelse(gene %in% neuropeptides$mouse_gene_name, "yes", "no"),
    is_alzheimer = ifelse(gene %in% alzheimer_all$mouse_gene_name, "yes", "no"),
    is_parkinson = ifelse(gene %in% parkinson_all$mouse_gene_name, "yes", "no"),
    is_als = ifelse(gene %in% als_all$mouse_gene_name, "yes", "no")) 
#write_tsv(DEGs_complete_info_final, "/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/9_241104_CAGs_psychiatric_association/250127_DEGs_complete_info_final")

#####################################################
# The DEGs table has been corrected (05/02/2025)
#####################################################
DEGs_complete_info_final <- readRDS("../250205_DEG_complete_results_kmeans.corrected.4243.rds")


# ========== 3. MA plots across contrasts ==========
# select_genes_df <- left_join(DEG_complete_results, DEGs_complete_info, by = "gene") %>%
#   dplyr::filter(logCPM > 5, significant_and_5per != "No significant") %>% 
#   dplyr::select(gene, contrast, significant, logFC, memory_class, direction, neurodegenerative_disease, psychotic_disorders, CAGs_n583, is_neuropeptide, other_adddiction) %>% 
#   dplyr::filter(contrast == all_contrasts[1])


downreg_genes <- c("Fancg", "Map2k1", "Ntsr1", "Drd2", "Grik1")
upreg_genes <- c("Dtx3", "Bdnf", "Cartpt", "Vip", "Penk", "Ucn", "Ptprt")
selected_DEGs <- c(downreg_genes, upreg_genes)

#Establish all possible contrasts:
all_contrasts <- c("h1_cocaine-saline",
                   "h4_cocaine-saline",
                   "h8_cocaine-saline",
                   "h24_cocaine-saline",
                   "d14_cocaine-saline", 
                   "h4_cocaine-h1_cocaine",
                   "h8_cocaine-h4_cocaine", 
                   "h24_cocaine-h8_cocaine", 
                   "d14_cocaine-h24_cocaine")
condition_colors <- c("black", "#9F2365", "#617641", "#C48208", "#326186", "#AE430A", "#564686")
condition_names <- c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")
names(condition_colors) <- condition_names

# Create a list to store all plots
MAplots <- list()

# Loop through each contrast
for (i in seq_along(all_contrasts)) {
  comparison <- all_contrasts[i]
  plot.data <- DEG_complete_results %>%
    filter(contrast == comparison) %>%
    mutate(to_label = ifelse(gene %in% selected_DEGs & significant != "No significant", "yes", "no"))
  
  query <- str_split(comparison, "-")[[1]][1]
  control <- str_split(comparison, "-")[[1]][2]
  
  # Define color values for the plot based on the current contrast
  volcano_colors <- c(condition_colors[query], condition_colors[control], "gray")
  names(volcano_colors) <- c("Upregulated (pval < 0.05 and logFC > 0.5)", "Downregulated (pval < 0.05 and logFC < -0.5)", "No significant")
  
  # Generate the plot
  MA_plot <- ggplot(plot.data, aes(x = logCPM, y = logFC)) +
    geom_point(aes(color = significant), size = 1.5) +
    scale_color_manual(values = volcano_colors) +
    theme_classic(base_size = 12) +
    geom_hline(yintercept = c(-0.5, 0.5), colour = "goldenrod", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "goldenrod") +
    ylim(-6, 8) +
    labs(title = glue("{query} vs {control}"),
         x = "Log (CPM)", y = "Log (FC)") +
    theme(legend.position = "none", legend.title = element_blank(), text = element_text(size = 15)) +
    geom_label_repel(data = subset(plot.data, plot.data$to_label != "no"), aes(label = gene), col = "black", box.padding = 1)
  
  # Store the plot in the list
  MAplots[[paste("MAplot", i, sep = "_")]] <- MA_plot
}

# Example to print a specific plot
MA_row1 <- MAplots$MAplot_1 + MAplots$MAplot_2 + MAplots$MAplot_3 + MAplots$MAplot_4 + MAplots$MAplot_5 +
  plot_layout(nrow = 1)
MA_row2 <- MAplots$MAplot_6 + MAplots$MAplot_7 + MAplots$MAplot_8 + MAplots$MAplot_9+ MAplots$MAplot_1 +
  plot_layout(nrow = 1)

MA_row1/MA_row2

# #Save
# dev.size()
# width_original = 20.572917
# height_original= 7.145833
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter1_DEGs_MAplots"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# MA_row1/MA_row2
# dev.off()



# ========== 4. Heatmap of cocaine vs saline DEGs ==========
saline_contrasts <-c("h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline", "h24_cocaine-saline", "d14_cocaine-saline")
time_contrasts <- c("h1_cocaine-saline", "h4_cocaine-h1_cocaine", "h8_cocaine-h4_cocaine", "h24_cocaine-h8_cocaine", "d14_cocaine-h24_cocaine")

samples_totest <- c("h4_saline_R1", "h8_saline_R1", "h24_saline_R1", "d14_saline_R1", #salines
                    "h1_cocaine_R1", "h1_cocaine_R2", "h4_cocaine_R1", "h4_cocaine_R2", "h8_cocaine_R1", "h8_cocaine_R2", #ETP
                    "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3") #LTP
condition <- c("saline", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")


clusterID = "VTA"

DNs.RNA.seu #2411 cells; "region" separates VTA from SN cells
seu.temporal <- DNs.RNA.seu[, DNs.RNA.seu$region == clusterID]

#Calculate average expression and format columns
Idents(seu.temporal) <- "orig.ident"
seu.temporal.avgexp = as.data.frame(AverageExpression(seu.temporal, group.by = 'orig.ident')) 
names(seu.temporal.avgexp) <- substring(names(seu.temporal.avgexp), 5)

# DEG_perGroup_results <- DEG_complete_results %>% filter(cluster_id == "VTA", contrast %in% saline_contrasts)
DEG_perGroup_results <- DEG_complete_results %>% filter(cluster_id == clusterID)

heatmap_genes <- DEG_complete_results %>%
  dplyr::filter(contrast %in% all_contrasts,
                significant != "No significant",  #Changed this bit to include also < 5%
                cluster_id == clusterID) %>% dplyr::select("gene") %>% pull()



seu.temporal <- DNs.RNA.seu[, DNs.RNA.seu$region == clusterID]
# seu.temporal@assays[["ATAC"]] <- NULL
# seu.temporal@assays[["SCT"]] <- NULL
# seu.temporal@assays[["peaks"]] <- NULL

#Calculate average expression and format columns
Idents(seu.temporal) <- "orig.ident"
seu.temporal.avgexp = as.data.frame(AverageExpression(seu.temporal, group.by = 'orig.ident')) 
names(seu.temporal.avgexp) <- substring(names(seu.temporal.avgexp), 5)

DEG_perGroup_results <- DEG_complete_results %>% filter(cluster_id == clusterID, contrast %in% saline_contrasts) 

#Extract DEGs from the big table and create metadata
DEG_withIDs <- DEG_perGroup_results %>% 
  dplyr::filter(significant != "No significant") %>%#%>% distinct(gene, .keep_all = TRUE) %>% 
  dplyr::mutate(groupID = ifelse(logFC < 0, glue("{control}_{query}"), query)) %>% 
  dplyr::group_by(gene) %>%
  dplyr::filter(abs(logFC) == max(abs(logFC)))  %>% 
  dplyr::ungroup() %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  dplyr::select(gene, contrast, p_val, logFC, query, control, significant, groupID) %>% 
  dplyr::group_by(gene) %>%
  dplyr::filter(abs(logFC) == max(abs(logFC))) %>% 
  dplyr::mutate(groupID = factor(groupID, levels = c(paste0("saline_", condition_names[3:7]), condition_names[3:7]))) %>% 
  dplyr::ungroup()

desired_order <- c(condition[2:6], paste0("saline_", condition[2:6]))

#Get average expression    
seu.DEGs.avgexp.heatmap <- seu.temporal.avgexp %>% 
  rownames_to_column("gene") %>% 
  left_join(DEG_withIDs[, c("gene", "groupID")], by = "gene") %>% 
  filter(gene %in% DEG_withIDs$gene) %>% 
  rename_with(~ gsub("\\.", "_", .)) %>% 
  dplyr::select(all_of(samples_totest), gene, groupID) %>%
  dplyr::mutate(groupID = factor(groupID, levels = desired_order)) %>%
  dplyr::arrange(groupID) %>% 
  column_to_rownames("gene")

sample_order <- colnames(seu.DEGs.avgexp.heatmap)[1:16]
sample_order <- gsub("\\.", "_", sample_order)
cluster_cols_vector <- c(rep("saline", 4), rep("h1_cocaine", 2), rep("h4_cocaine", 2),  rep("h8_cocaine",2),  rep("h24_cocaine", 3), rep("d14_cocaine", 3))
cluster_cols_vector <- factor(cluster_cols_vector, levels = c("saline", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine"))
cluster_rows_vector <- seu.DEGs.avgexp.heatmap$groupID


#Prepare matrix for heatmap (in the order I want to be plotted)
mat <- as.matrix(seu.DEGs.avgexp.heatmap[, 1:16])
mat_scaled = t(scale(t(mat)))[, sample_order]


sample_colors <- condition_colors[cluster_cols_vector]

gene_groupings <- seu.DEGs.avgexp.heatmap$groupID
groupings_cols <- c("saline_h1_cocaine" = lighten("#617641", 0.5), "saline_h4_cocaine" = lighten("#C48208", 0.5), "saline_h8_cocaine" = lighten("#326186", 0.5), "saline_h24_cocaine" = lighten("#AE430A", 0.5), "saline_d14_cocaine" = lighten("#564686", 0.5),
                    "h1_cocaine" = "#617641", "h4_cocaine" = "#C48208", "h8_cocaine" = "#326186", "h24_cocaine" = "#AE430A", "d14_cocaine" = "#564686")


ha = HeatmapAnnotation(Condition = cluster_cols_vector, show_legend = FALSE, 
                       col = list(Condition = condition_colors), show_annotation_name = FALSE,labels = NULL)

row_ha = rowAnnotation(Group = gene_groupings, show_legend = FALSE, 
                       col = list(Group = groupings_cols, show_annotation_name = FALSE, simple_anno_size = unit(8, "mm")))
heatmap_color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)  

#Now, add some extra features: CAGs and psychotic genes
# 1. get order of genes in the heatmap:
cocaine_vs_saline_heatmap <- 
  Heatmap(mat_scaled, split = cluster_rows_vector, cluster_rows = FALSE, column_split = cluster_cols_vector, 
          cluster_row_slices = FALSE, cluster_column_slices = FALSE,  
          top_annotation = ha, 
          right_annotation = row_ha,
          name = "Z-score", row_title = NULL, show_row_names = FALSE,  use_raster = FALSE, 
          col = heatmap_color, 
          show_row_dend = FALSE, show_column_dend = FALSE, 
          heatmap_legend_param = list(title = "Z-score", direction = "horizontal"), 
          show_column_names = FALSE)
ht_row_order <- do.call(rbind, lapply(names(row_order(cocaine_vs_saline_heatmap)), function(group) {
  data.frame(index = row_order(cocaine_vs_saline_heatmap)[[group]], groupID = group)
}))
ht_row_order <- as.data.frame(ht_row_order)
mat_scaled_genes_sorted <- rownames(mat_scaled)[ht_row_order$index]

selected_genes_to_plot <- c("Ncor2", "Dtx3", "Ace", "Cdh9", #h1 
                            "Akt1","Bdnf", "Ptprt", # h4
                            "Pde7b", "Vip", "Cartpt", "Ucn", #h8
                            "Penk", "Cers2", #h24
                            "Npy", "Gal", "Il7", "Gbrg1", "Gabra5",  #d14
                            "Tacr3", "Xrcc5", "Polr2a", "Fancg", "Abca1",  #sal vs h1
                            "Map2k1", "Nrp1", "Csnk2b", "Ctdp1", "Grk6", #sal vs h4
                            "Ldlr", "Slc27a1", "Npy1r", "Ptgds", "Chrna6", #sal vs h8
                            "Egf", "Src", "Ctnna1", "Drd2", #sal vs h24
                            "Coq3", "Grik1")  #sal vs d14
indices <- which(mat_scaled_genes_sorted %in% selected_genes_to_plot)
#DEGs_complete_info[DEGs_complete_info$gene %in% selected_genes_to_plot, c(4, 15:21)] %>% filter(is_CAG == "yes")
DEGs_complete_info <- DEGs_complete_info_final


#2. Get the TRUE/FALSE for each gene in each category:
CAGs_degs <-  DEGs_complete_info %>% filter(is_CAG != "no") %>% pull(gene)
other_addiction_degs <-  DEGs_complete_info %>% filter(is_addiction != "no") %>% pull(gene)
psycho_degs <- DEGs_complete_info %>% filter(is_psychiatric != "no") %>% pull(gene)
neuropeptides_degs <- DEGs_complete_info %>% filter(is_neuropeptide != "no") %>% pull(gene) 
alzheimer_degs <- DEGs_complete_info %>% filter(is_alzheimer != "no") %>% pull(gene)
parkinson_degs <- DEGs_complete_info %>% filter(is_parkinson != "no") %>% pull(gene)
als_degs <- DEGs_complete_info %>% filter(is_als != "no") %>% pull(gene)

is_CAG <- mat_scaled_genes_sorted %in% CAGs_degs
is_other_addiction <- mat_scaled_genes_sorted %in% other_addiction_degs
is_psycho <- mat_scaled_genes_sorted %in% psycho_degs 
is_neuropeptide <- mat_scaled_genes_sorted %in% neuropeptides_degs
is_alzheimer <- mat_scaled_genes_sorted %in% alzheimer_degs
is_parkinson <- mat_scaled_genes_sorted %in% parkinson_degs
is_als <- mat_scaled_genes_sorted %in% als_degs


#Plot everything
ht_opt$ROW_ANNO_PADDING = unit(0.5, "cm")

row_ha = rowAnnotation(Group = gene_groupings, show_legend = FALSE, 
                       col = list(Group = groupings_cols, show_annotation_name = FALSE, simple_anno_size = unit(8, "mm")), 
                       foo = anno_mark(at = indices, 
                                       labels = mat_scaled_genes_sorted[indices]))
cocaine_vs_saline_heatmap <- 
  Heatmap(mat_scaled, split = cluster_rows_vector, cluster_rows = FALSE, column_split = cluster_cols_vector, 
          cluster_row_slices = FALSE, cluster_column_slices = FALSE,  
          top_annotation = ha, 
          right_annotation = row_ha,
          name = "Z-score", row_title = NULL, show_row_names = FALSE,  use_raster = FALSE, 
          col = heatmap_color, 
          show_row_dend = FALSE, show_column_dend = FALSE, 
          heatmap_legend_param = list(title = "Z-score", direction = "horizontal"), 
          show_column_names = FALSE) +
  Heatmap(is_CAG + 0, name = "CAGs", col = c("0" = "white", "1" = "black"), 
          show_heatmap_legend = FALSE, width = unit(10, "mm")) + 
  Heatmap(is_other_addiction + 0, name = "other addiction", col = c("0" = "white", "1" = "black"), 
          show_heatmap_legend = FALSE, width = unit(10, "mm")) +
  Heatmap(is_psycho + 0, name = "psychiatric", col = c("0" = "white", "1" = "black"), 
          show_heatmap_legend = FALSE, width = unit(10, "mm")) +
  Heatmap(is_neuropeptide + 0, name = "neuropeptide", col = c("0" = "white", "1" = "black"), 
          show_heatmap_legend = FALSE, width = unit(10, "mm")) +
  Heatmap(is_alzheimer + 0, name = "alzheimer", col = c("0" = "white", "1" = "black"), 
          show_heatmap_legend = FALSE, width = unit(10, "mm")) +
  Heatmap(is_parkinson + 0, name = "parkinson", col = c("0" = "white", "1" = "black"), 
          show_heatmap_legend = FALSE, width = unit(10, "mm")) +
  Heatmap(is_als + 0, name = "als", col = c("0" = "white", "1" = "black"), 
          show_heatmap_legend = FALSE, width = unit(10, "mm"))
cocaine_vs_saline_heatmap


# dev.size()
# 
# current_width <- 8
# current_height <- 11.84375

# # Set the resolution in DPI
# dpi <- 300
# 
# # Calculate pixel dimensions based on current device size and DPI
# width_pixels <- current_width * dpi
# height_pixels <- current_height * dpi
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter2_GEX_cocaine_vs_saline_heatmaps_v2"
# file_name <- glue("{pdf_dir}/{plot_name}.tiff" )


# Save the heatmap as a TIFF file with calculated dimensions
# tiff(file_name, width = width_pixels, height = height_pixels, res = dpi)
# draw(cocaine_vs_saline_heatmap) # Replace with your heatmap drawing function
# dev.off()
#  
# width_original <- 9.5
# height_original <- 11.84375
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter2_GEX_cocaine_vs_saline_heatmaps_v2"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# draw(cocaine_vs_saline_heatmap)
# dev.off()


# ========== 5. Classifier statistics and barplots ==========
#Bar plots:
ht_row_order_withFeatures <- cbind(ht_row_order, gene = mat_scaled_genes_sorted, isCAG = is_CAG, isPyscho = is_psycho) %>%
  dplyr::mutate(
    direction = ifelse(startsWith(as.character(groupID), "saline"), "down", "up"),
    timepoint = ifelse(startsWith(as.character(groupID), "saline"), sub("^saline_", "", as.character(groupID)), as.character(groupID)))

#Save files for enrichR (only for figure)
pairwise_genes_enrichR <- seu.DEGs.avgexp.heatmap %>% 
  rownames_to_column("gene") %>% 
  dplyr::select(gene, groupID)
grouped_genes <- split(pairwise_genes_enrichR$gene, pairwise_genes_enrichR$groupID)

# for(group in names(grouped_genes)) {
#   filename <- paste0(group, "_genes.txt")
#   writeLines(grouped_genes[[group]], con = filename)
# }
# writeLines(pairwise_genes_enrichR$gene, con = "background.txt")

#Define the desired reverse order for time points
timepoint_order <- c("d14_cocaine", "h24_cocaine", "h8_cocaine", "h4_cocaine", "h1_cocaine")

#Calculate the percentage of "isCAG" genes for each timepoint and direction
data_summary <- ht_row_order_withFeatures %>%
  group_by(timepoint, direction) %>%
  dplyr::summarise(total_genes = n(),
                   isCAG_genes = sum(isCAG == TRUE)) %>%
  mutate(percentage = round((isCAG_genes / total_genes) * 100, 1))

#Ensure the data includes all relevant time points and is ordered correctly
data_summary$timepoint <- factor(data_summary$timepoint, levels = timepoint_order)

#Create the bidirectional bar plot with percentage labels
plota <- ggplot(data_summary, aes(x = factor(timepoint, levels = timepoint_order), y = ifelse(direction == "up", percentage, -percentage), fill = direction)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  geom_text(aes(label = abs(percentage)), 
            position = position_stack(vjust = 0.5), 
            color = "black") +
  scale_y_continuous(labels = abs, limits = c(-5, 5)) +  # Set x-axis limits to -5 and +5
  coord_flip() +
  labs(title = "",
       x = "",
       y = "Percentage of Cocaine Assiciated Genes (CAGs)") +
  scale_fill_manual(values = c("up" = "#E75252", "down" = "#5290C1")) +no_legend()


#Calculate the percentage of "isPyscho" genes for each timepoint and direction
data_summary_psycho <- dim(ht_row_order_withFeatures) %>%
  group_by(timepoint, direction) %>%
  dplyr::summarise(total_genes = n(),
                   isPyscho_genes = sum(isPyscho == TRUE)) %>%
  mutate(percentage = round((isPyscho_genes / total_genes) * 100, 1))

#Ensure the data includes all relevant time points and is ordered correctly
data_summary_psycho$timepoint <- factor(data_summary_psycho$timepoint, levels = timepoint_order)

#Create the bidirectional bar plot with percentage labels for isPyscho genes
plotb <- ggplot(data_summary_psycho, aes(x = factor(timepoint, levels = timepoint_order), y = ifelse(direction == "up", percentage, -percentage), fill = direction)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  geom_text(aes(label = abs(percentage)), 
            position = position_stack(vjust = 0.5), 
            color = "black") +
  scale_y_continuous(labels = abs, limits = c(-12, 12)) +  # Set x-axis limits to -5 and +5
  coord_flip() +
  labs(title = "",
       x = "",
       y = "Percentage of psychiatric disorder genes") +
  scale_fill_manual(values = c("up" = "#E75252", "down" = "#5290C1")) +
  theme(
    axis.text.y = element_blank(),   # Remove Y-axis text
    axis.ticks.y = element_blank()   # Remove Y-axis ticks
  )
plota + plotb

# dev.size()
# width_original <- 9.714286
# height_original <- 2.771429
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter2_classifiers_upanddown_contrast"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# plota + plotb
# dev.off()

#Plots:
#########
DEGs_complete_info_Classifiers <- ht_row_order_withFeatures[, c("index", "groupID", "gene", "direction", "timepoint")] %>% 
  left_join(DEGs_complete_info[, c("gene", "is_CAG", "is_addiction", "is_psychiatric", "is_neuropeptide", "is_alzheimer", "is_parkinson", "is_als")])

data <- DEGs_complete_info_Classifiers
data$direction <- factor(data$direction, levels = c("up", "down"))

#Calculate percentages for each category and direction
#First, calculate the total number of upregulated and downregulated genes
total_DEGs <- data %>%
  group_by(direction) %>%
  dplyr::summarise(total = n(), .groups = "drop")

#Then calculate the count of 'yes' statuses for each classifier within each direction
percentages <- data %>%
  pivot_longer(cols = c("gene", "is_CAG", "is_addiction", "is_psychiatric", "is_neuropeptide", "is_alzheimer", "is_parkinson", "is_als"),
               names_to = "classifier", values_to = "status") %>%
  filter(status != "no") %>%
  group_by(direction, classifier) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  left_join(total_DEGs, by = "direction") %>%
  mutate(percentage = count / total * 100) %>% 
  filter(!is.na(percentage)) %>% 
  complete(direction, classifier, fill = list(count = 0, percentage = 0)) %>%   # Fill missing combinations
  mutate(total = if_else(is.na(total) & direction == "up" & classifier == "in_als", 
                         sum(total_DEGs$total[total_DEGs$direction == "up"]),  # Assuming total_DEGs has total for "up"
                         total),
         percentage = count / total * 100) %>% 
  filter(!is.na(percentage), classifier !="gene")  

classifier_order <- c("gene", "is_CAG", "is_addiction", "is_psychiatric", "is_neuropeptide", "is_alzheimer", "is_parkinson", "is_als")
percentages$classifier <- factor(percentages$classifier, levels = classifier_order)

max_percentage <- max(percentages$percentage)

#Percentage of DEGs Associated with Neurological and Addiction-Related Classifiers
percentages <- percentages %>% filter(classifier != "is_neuropeptide")
plot1 <- ggplot(percentages, aes(x = direction, y = percentage, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_text(aes(label = count), position = position_dodge(width = 0.7), vjust = -0.25) +
  facet_wrap(~ classifier, scales = "free_y", nrow=1) +
  labs(title = "",
       x = "Direction",
       y = "Percentage (%)") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
  scale_fill_manual(values = c("up" = "#EC5D5B", "down" = "#89B8D2")) +
  theme(legend.position = "none") 
plot1

# dev.size()
# width_original <- 11.979167
# height_original <- 4.666667
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter1_classifiers_upanddown_over_all_v3"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# plot1
# dev.off()

#Plot 2:
percentages <- percentages %>%
  dplyr::mutate(total_classifier = case_when(
    classifier == "is_CAG" ~  length(cags_all),
    classifier == "is_addiction" ~ length(unique(general_addiction_all)),
    classifier == "is_psychiatric" ~ length(unique(psychiatric_all)),
    classifier == "is_neuropeptide" ~ length(unique(neuropeptides$mouse_gene_name)),
    classifier == "is_alzheimer" ~ length(unique(alzheimer_all$mouse_gene_name)),
    classifier == "is_parkinson" ~ length(unique(parkinson_all$mouse_gene_name)),
    classifier == "is_als" ~ length(unique(als_all$mouse_gene_name)),
    TRUE ~ NA_real_
  )) %>%
  dplyr::mutate(proportion_over_all = (count / total_classifier) * 100)
max_percentage <- max(percentages$proportion_over_all) + 1

plot2 <- ggplot(percentages, aes(x = direction, y = proportion_over_all, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", proportion_over_all)), position = position_dodge(width = 0.7), vjust = -0.25) +
  facet_wrap(~ classifier, scales = "free_y", nrow=1) +
  labs(title = "Proportion of known addiction/neurological disorder genes mis-regulated after one dose of cocaine",
       x = "Direction",
       y = "Percentage (%)") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(0, max_percentage), breaks = seq(0, max_percentage, 5)) +
  scale_fill_manual(values = c("up" = "#EC5D5B", "down" = "#89B8D2")) +
  theme(legend.position = "right") 
plot2

# dev.size()
# width_original <- 11.752381
# height_original <- 3.257143
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter2_classifiers_upanddown_over_all"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# plot2
# dev.off()



# ========== 6. UpSet plots ==========
##Upset plot
##############
classifiers <-  c("is_CAG", "is_addiction", "is_psychiatric", "is_neuropeptide", "is_alzheimer", "is_parkinson", "is_als")

# Initialize lists
upregulated_upset <- list()
downregulated_upset <- list()

# Loop through classifiers to create lists
for (classifier in classifiers) {
  upregulated_upset[[classifier]] <- DEGs_complete_info_Classifiers %>%
    filter(direction == "up", .data[[classifier]] != "no") %>%
    pull(gene)
  
  downregulated_upset[[classifier]] <- DEGs_complete_info_Classifiers %>%
    filter(direction == "down", .data[[classifier]] != "no") %>%
    pull(gene)
}

m_up = make_comb_mat(upregulated_upset)
m_down = make_comb_mat(downregulated_upset)

upregulated_upset_plot <- UpSet(t(m_up), set_order =classifiers, comb_order = rev(order(comb_size(m_up))), 
                                top_annotation = upset_top_annotation(t(m_up), add_numbers = TRUE),
                                right_annotation = upset_right_annotation(t(m_up), add_numbers = TRUE)
)

downregulated_upset_plot <- UpSet(t(m_down), set_order =classifiers, comb_order = rev(order(comb_size(m_down))), 
                                  top_annotation = upset_top_annotation(t(m_down), add_numbers = TRUE),
                                  right_annotation = upset_right_annotation(t(m_down), add_numbers = TRUE))


# ========== 7. Significance tests ==========
no_degs <- genes[genes$gene %notin% data$gene, ] %>% pull(gene) 
upreg_genes <- data %>% filter(direction == "up") %>% dplyr::pull(gene) 
down_genes <- data %>% filter(direction == "down") %>% dplyr::pull(gene) 

# Get counts of CAGs in each group
cags_in_upreg <- length(intersect(upreg_genes, cags_all))
cags_in_downreg <- length(intersect(down_genes, cags_all))
cags_in_nodegs <- length(intersect(no_degs, cags_all))

# Total number of cocaine-associated genes
total_cags <- length(cags_all)

# Other genes in each group
non_cags_in_upreg <- length(upreg_genes) - cags_in_upreg
non_cags_in_downreg <- length(down_genes) - cags_in_downreg
non_cags_in_nodegs <- length(no_degs) - cags_in_nodegs

## Chi square for upregulated genes
##################################
# Contingency table
contingency_table <- matrix(c(
  cags_in_upreg, non_cags_in_upreg,    # Upregulated
  cags_in_nodegs, non_cags_in_nodegs  # No-DEGs
), nrow = 2, byrow = TRUE)

rownames(contingency_table) <- c("Upregulated", "No-DEGs")
colnames(contingency_table) <- c("CAGs", "Non-CAGs")

print(contingency_table)

# Chi-square test
chi_result <- chisq.test(contingency_table)
print(chi_result)

## Chi square for downregulated genes
####################################
# Contingency table
contingency_table <- matrix(c(
  cags_in_downreg, non_cags_in_downreg,    #Downregulated
  cags_in_nodegs, non_cags_in_nodegs  # No-DEGs
), nrow = 2, byrow = TRUE)

rownames(contingency_table) <- c("Downregulated", "No-DEGs")
colnames(contingency_table) <- c("CAGs", "Non-CAGs")

print(contingency_table)

# Chi-square test
chi_result <- chisq.test(contingency_table)
print(chi_result)



### Significance test for psychiatric_all:
##########################################
no_degs <- genes[genes$gene %notin% data$gene, ] %>% pull(gene) 
upreg_genes <- data %>% filter(direction == "up") %>% dplyr::pull(gene) 
down_genes <- data %>% filter(direction == "down") %>% dplyr::pull(gene) 

# Get counts of psychiatric-associated genes (PAGs) in each group
pags_in_upreg <- length(intersect(upreg_genes, psychiatric_all))
pags_in_downreg <- length(intersect(down_genes, psychiatric_all))
pags_in_nodegs <- length(intersect(no_degs, psychiatric_all))

# Total number of psychiatric-associated genes
total_pags <- length(psychiatric_all)

# Other genes in each group
non_pags_in_upreg <- length(upreg_genes) - pags_in_upreg
non_pags_in_downreg <- length(down_genes) - pags_in_downreg
non_pags_in_nodegs <- length(no_degs) - pags_in_nodegs

## Chi-square for upregulated genes
###################################
# Contingency table
contingency_table <- matrix(c(
  pags_in_upreg, non_pags_in_upreg,    # Upregulated
  pags_in_nodegs, non_pags_in_nodegs  # No-DEGs
), nrow = 2, byrow = TRUE)

rownames(contingency_table) <- c("Upregulated", "No-DEGs")
colnames(contingency_table) <- c("PAGs", "Non-PAGs")

print(contingency_table)

# Chi-square test
chi_result <- chisq.test(contingency_table)
print(chi_result)

## Chi-square for downregulated genes
#####################################
# Contingency table
contingency_table <- matrix(c(
  pags_in_downreg, non_pags_in_downreg,    # Downregulated
  pags_in_nodegs, non_pags_in_nodegs      # No-DEGs
), nrow = 2, byrow = TRUE)

rownames(contingency_table) <- c("Downregulated", "No-DEGs")
colnames(contingency_table) <- c("PAGs", "Non-PAGs")

print(contingency_table)

# Chi-square test
chi_result <- chisq.test(contingency_table)
print(chi_result)

# ========== 8. UpSet plots for figure panels ==========
prop_only_up <- percentages %>%
  mutate(total_minus_count = total_classifier - count) %>%
  filter(direction == "up") %>%
  group_by(classifier) %>%
  mutate(
    count_percentage = count / total_classifier * 100,
    total_minus_count_percentage = total_minus_count / total_classifier * 100
  ) %>%
  dplyr::select(classifier, count_percentage, total_minus_count_percentage) %>%
  pivot_longer(
    cols = c(count_percentage, total_minus_count_percentage),
    names_to = "type",
    values_to = "value"
  ) %>%
  mutate(type = ifelse(type == "count_percentage", "Upregulated", "Remaining"))

# Plotting the data
ggplot(prop_only_up, aes(x = classifier, y = value / 100, fill = type)) +  # Scale the y values to max out at 0.25
  geom_bar(stat = "identity", position = "stack", width = 0.6) +  # use position stack for regular bar stacking
  geom_text(
    data = filter(prop_only_up, type == "Upregulated"), 
    aes(label = sprintf("%.1f%%", value), y = value / 100 + 0.01), size = 5, # Adjust text placement slightly above the bars
    position = position_stack(vjust = 1),
    color = "black"
  ) +
  scale_fill_manual(values = c("Upregulated" = "#EC5D5B", "Remaining" = "grey")) +
  labs(x = "Classifier", y = "Proportion", fill = "Category") +
  ylim(0, 0.25) +  # Set y limits to scale the bars up to 25%
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


prop_only_down <- percentages %>%
  mutate(total_minus_count = total_classifier - count) %>%
  filter(direction == "down") %>%
  group_by(classifier) %>%
  mutate(
    count_percentage = count / total_classifier * 100,
    total_minus_count_percentage = total_minus_count / total_classifier * 100
  ) %>%
  dplyr::select(classifier, count_percentage, total_minus_count_percentage) %>%
  pivot_longer(
    cols = c(count_percentage, total_minus_count_percentage),
    names_to = "type",
    values_to = "value"
  ) %>%
  mutate(type = ifelse(type == "count_percentage", "Downregulated", "Remaining"))

# Plotting the data
ggplot(prop_only_down, aes(x = classifier, y = value / 100, fill = type)) +  # Scale the y values to max out at 0.25
  geom_bar(stat = "identity", position = "stack", width = 0.6) +  # use position stack for regular bar stacking
  geom_text(
    data = filter(prop_only_down, type == "Downregulated"), 
    aes(label = sprintf("%.1f%%", value), y = value / 100 + 0.01), size = 5, # Adjust text placement slightly above the bars
    position = position_stack(vjust = 1),
    color = "black"
  ) +
  scale_fill_manual(values = c("Downregulated" = "#89B8D2", "Remaining" = "grey")) +
  labs(x = "Classifier", y = "Proportion", fill = "Category") +
  ylim(0, 0.25) +  # Set y limits to scale the bars up to 25%
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



##Upset plot: For main figure. 
##############
classifiers <- c("is_CAG", "is_addiction", "is_psychiatric")

# Initialize lists
upregulated_upset <- list()
downregulated_upset <- list()

# Loop through classifiers to create lists
for (classifier in classifiers) {
  upregulated_upset[[classifier]] <- DEGs_complete_info_Classifiers %>%
    filter(direction == "up", .data[[classifier]] != "no") %>%
    pull(gene)
  
  downregulated_upset[[classifier]] <- DEGs_complete_info_Classifiers %>%
    filter(direction == "down", .data[[classifier]] != "no") %>%
    pull(gene)
}

m_up = make_comb_mat(upregulated_upset)
m_down = make_comb_mat(downregulated_upset)


upregulated_upset_plot <- UpSet(t(m_up), set_order =classifiers, comb_order = rev(order(comb_size(m_up))), 
                                top_annotation = upset_top_annotation(t(m_up), add_numbers = TRUE),
                                right_annotation = upset_right_annotation(t(m_up), add_numbers = TRUE))

downregulated_upset_plot <- UpSet(t(m_down), set_order =classifiers, comb_order = rev(order(comb_size(m_down))), 
                                  top_annotation = upset_top_annotation(t(m_down), add_numbers = TRUE),
                                  right_annotation = upset_right_annotation(t(m_down), add_numbers = TRUE))


#Second version
desired_order <- comb_size(m_up)[c("111", "110", "101", "100", "011", "010", "001")]
desired_order_up <- c(1,2,3,5,4,6,7)

upregulated_upset_plot <- UpSet(t(m_up), set_order =classifiers, comb_order = desired_order_up, 
                                top_annotation = upset_top_annotation(t(m_up), add_numbers = TRUE),
                                right_annotation = upset_right_annotation(t(m_up), add_numbers = TRUE))

extract_comb(m_down, "110")

desired_order <- comb_size(m_up)[c("111", "110", "100", "011", "010", "001")]
desired_order_down <- c(1,2,4,3,5,6)
downregulated_upset_plot <- UpSet(t(m_down), set_order =classifiers, comb_order = desired_order_down, 
                                  top_annotation = upset_top_annotation(t(m_down), add_numbers = TRUE),
                                  right_annotation = upset_right_annotation(t(m_down), add_numbers = TRUE))

downregulated_upset_plot



# dev.size()
# width_original <- 4.114286
# height_original <- 3.980952
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter2_upregulated_upset_plot_classifiers"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# upregulated_upset_plot
# dev.off()



# dev.size()
# width_original <- 4.114286
# height_original <- 3.980952
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter2_downregulated_upset_plot_classifiers"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# downregulated_upset_plot
# dev.off()


# ========== 9. Final overlap table ==========
### INSERT CODE HERE ###

# ========== 10. Session information ==========
# saveRDS(sessionInfo(), "sessionInfo.rds")