# ===========================================
# Script Title: 1.4 MAplots and heatmaps
# Author: Luna Zea Redondo
# Date: 2026-05-13
# Description:
#   This script represents the DEG results through MA plots and a heatmap for 
#   all contrasts, based on the results from the previous script (1.3)
#   and highlight some specific genes. It also evaluates significant enrichment of DEGs
#   over several gene lists (CAGs, psychiatric, etc) and plot the results. 
# ===========================================

# ========== Set Environment ==========
rm(list = ls(all.names = TRUE))
gc()
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

source("/fast/AG_Pombo/luna/2026_rebuttal/config.R", local = FALSE)

# ========== Load Required Libraries ==========
library(dplyr)
library(tidyr)
library(tibble)
library(glue)
library(ggplot2)
library(patchwork)
library(reshape2)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(UpSetR)
library(colorspace)


# ========== Set Directory ==========
dir <- "/fast/AG_Pombo/luna/2026_rebuttal/4_DEGs_visualization/"
setwd(dir)


# ========== 1. Load the required files ==========
DEG_complete_results <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/3_muscat/deg_results/260513_DEG_complete_results_Analysis1.tsv")
seu.VTA_DNs <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/2_DN-evaluation/250507_seu.VTA_DNs.rds")

DEG_complete_results <- DEG_complete_results %>%
  filter(!str_detect(gene, "^mt-"))

DEG_complete_results$more_10cells <- "no"

present_genes <- DEG_complete_results$gene %in% rownames(seu.VTA_DNs)

DEG_complete_results$more_10cells[present_genes] <- ifelse(
  rowSums(
    GetAssayData(seu.VTA_DNs, slot = "counts")[
      DEG_complete_results$gene[present_genes], 
    ] > 0
  ) >= 10,
  "yes",
  "no"
)

# ========== 2. Harmonize classifiers and add to DEG info ==========
DEG_complete_results_simplified <- DEG_complete_results %>% 
  dplyr::select(gene, logFC, logCPM, p_val, p_adj.loc, p_adj.loc, contrast, query, control, significant, percent5, significant_and_5per)

#Add mouse to human conversion
mouse_human_conversion <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/human_mouse_conversion.txt") %>% 
  dplyr::select(human_gene_name, mouse_gene_name) 
mouse_human_conversion$human_gene_name[mouse_human_conversion$mouse_gene_name == "Hnrnpa1"] <- "HNRNPA1"

#Harmonize gene lists and add to DEG_complete_results_simplified
gene_lists_dir <- "/fast/AG_Pombo/luna/2026_rebuttal/gene_lists"

#a. List of CAGs
cags_all <- read_tsv(glue("{gene_lists_dir}/dominik_CAGs_harmonized_n584.tsv")) %>%
  dplyr::select(gene_name) %>% dplyr::pull(gene_name)

#2. Psychiatric disorders:
psychiatric_all <- read_tsv(glue("{gene_lists_dir}/240826_psychotic_disorders_final.tsv")) %>%
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

#3. General addiction:
iza_anco_updated <- read_tsv(glue("{gene_lists_dir}/ANCO_genes_curated_to_mouse_genes.txt")) %>% 
  dplyr::select(gene_name, addiction_phenotype) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(addiction_category = if_else(
    addiction_phenotype == "not", "not", "addiction"))
general_addiction_all <- iza_anco_updated %>% dplyr::filter(addiction_category == "addiction") %>% dplyr::pull(gene_name)

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

DEGs_complete_info_simplified <- DEG_complete_results_simplified %>%
  dplyr::mutate(
    is_CAG = ifelse(gene %in% cags_all, "yes", "no"), 
    is_addiction = ifelse(gene %in% general_addiction_all, "yes", "no"), 
    is_psychiatric = ifelse(gene %in% psychiatric_all, "yes", "no"), 
    is_neuropeptide = ifelse(gene %in% neuropeptides$mouse_gene_name, "yes", "no"),
    is_alzheimer = ifelse(gene %in% alzheimer_all$mouse_gene_name, "yes", "no"),
    is_parkinson = ifelse(gene %in% parkinson_all$mouse_gene_name, "yes", "no"),
    is_als = ifelse(gene %in% als_all$mouse_gene_name, "yes", "no")) 

#write_tsv(DEGs_complete_info_simplified, "260513_DEGs_complete_info_simplified.tsv")



# ========== 2. Generate MA Plots ==========
# Manual selection of genes to highlight
downreg_genes <- c("Fancg", "Map2k1", "Ntsr1", "Drd2", "Grik1", "Foxa1")
upreg_genes <- c("Dtx3", "Bdnf", "Cartpt", "Vip", "Penk", "Ptprt", "Npy")

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
    dplyr::filter(contrast == comparison) %>%
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
    ylim(-7, 9) +
    labs(title = "",
         x = "", y = "") +
    theme(legend.position = "none", legend.title = element_blank(), text = element_text(size = 15)) +
    geom_label_repel(data = subset(plot.data, plot.data$to_label != "no"), aes(label = gene), max.overlaps = Inf, col = "black", box.padding = 1, size = 5) +
    labs(x = NULL, y = NULL) +  # Remove axis labels
    theme(
      axis.title.x = element_blank(),  # Remove x-axis title
      axis.title.y = element_blank(),  # Remove y-axis title
      axis.text.x = element_blank(),   # Remove x-axis tick labels
      axis.text.y = element_blank()    # Remove y-axis tick labels
    )
  
  # Store the plot in the list
  MAplots[[paste("MAplot", i, sep = "_")]] <- MA_plot
}

# Example to print a specific plot
MA_row1 <- MAplots$MAplot_1 + MAplots$MAplot_2 + MAplots$MAplot_3 + MAplots$MAplot_4 + MAplots$MAplot_5 +
  plot_layout(nrow = 1)
MA_row2 <- MAplots$MAplot_1 + MAplots$MAplot_6 + MAplots$MAplot_7 + MAplots$MAplot_8 + MAplots$MAplot_9+
  plot_layout(nrow = 1)

# dev.size()
# tiff("chapter1_MAplots_DEGs_v2.tif", width = 19.156250, height = 8.291667, res = 300, units = "in")
# print(MA_row1/MA_row2)
# dev.off()


# ========== 3. Generate Heatmap ==========

#Extablisg contrasts
saline_contrasts <-c("h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline", "h24_cocaine-saline", "d14_cocaine-saline")
time_contrasts <- c("h1_cocaine-saline", "h4_cocaine-h1_cocaine", "h8_cocaine-h4_cocaine", "h24_cocaine-h8_cocaine", "d14_cocaine-h24_cocaine")

samples_totest <- c("h4_saline_R1", "h8_saline_R1", "h24_saline_R1", "d14_saline_R1", #salines
                    "h1_cocaine_R1", "h1_cocaine_R2", "h4_cocaine_R1", "h4_cocaine_R2", "h8_cocaine_R1", "h8_cocaine_R2", #ETP
                    "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3") #LTP
condition <- c("saline", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")


seu.VTA_DNs #1520 cells; "region" separates VTA from SN cells
seu.temporal <- seu.VTA_DNs

#Calculate average expression and format columns
Idents(seu.temporal) <- "orig.ident"
seu.temporal.avgexp = as.data.frame(AverageExpression(seu.temporal, group.by = 'orig.ident')) 
names(seu.temporal.avgexp) <- substring(names(seu.temporal.avgexp), 5)

# DEG_perGroup_results <- DEG_complete_results %>% filter(cluster_id == "VTA", contrast %in% saline_contrasts)
# DEG_perGroup_results <- DEG_complete_results %>% filter(cluster_id == clusterID)

heatmap_genes <- DEG_complete_results %>%
  dplyr::filter(contrast %in% all_contrasts,
                significant != "No significant") %>%   #Changed this bit to include also < 5%
                 dplyr::select("gene") %>% pull()


# seu.temporal@assays[["ATAC"]] <- NULL
# seu.temporal@assays[["SCT"]] <- NULL
# seu.temporal@assays[["peaks"]] <- NULL

#Calculate average expression and format columns
Idents(seu.temporal) <- "orig.ident"
seu.temporal.avgexp = as.data.frame(AverageExpression(seu.temporal, group.by = 'orig.ident')) 
names(seu.temporal.avgexp) <- substring(names(seu.temporal.avgexp), 5)

DEG_perGroup_results <- DEGs_complete_info_simplified %>% dplyr::filter(contrast %in% saline_contrasts) 

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
  dplyr::filter(gene %in% DEG_withIDs$gene) %>% 
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

selected_genes_to_plot <- c("Htr2a", "Ace", "Dtx3", #h1 
                            "Akt1","Bdnf", "Ptprt", # h4 
                            "Slc22a3", "Vip", "Drd1", #h8 
                            "Penk", "Cartpt", "Ucn", "Npy", #h24 
                            "Gdnf", "Il7", "Col23a1",  #d14 
                            "Tacr3", "Xrcc5", "Fancg", "Abca1",  #sal vs h1 
                            "Map2k1", "Nrp1", "Csnk2b", "Ctdp1", "Grk6", #sal vs h4 
                            "Foxa1", "Npy1r", "Chrna3", #sal vs h8  
                            "Egf", "Src", "Drd2", "Ptgds", #sal vs h24
                            "Coq3", "Grik1", "Polr2a")  #sal vs d14

indices <- which(mat_scaled_genes_sorted %in% selected_genes_to_plot)


#2. Get the TRUE/FALSE for each gene in each category:
CAGs_degs <-  DEGs_complete_info_simplified %>% dplyr::filter(is_CAG != "no") %>% pull(gene)
other_addiction_degs <-  DEGs_complete_info_simplified %>% dplyr::filter(is_addiction != "no") %>% pull(gene)
psycho_degs <- DEGs_complete_info_simplified %>% dplyr::filter(is_psychiatric != "no") %>% pull(gene)
neuropeptide_degs <- DEGs_complete_info_simplified %>% dplyr::filter(is_neuropeptide != "no") %>% pull(gene)
alzheimer_degs <- DEGs_complete_info_simplified %>% dplyr::filter(is_alzheimer != "no") %>% pull(gene)
parkinson_degs <- DEGs_complete_info_simplified %>% dplyr::filter(is_parkinson != "no") %>% pull(gene)
als_degs <- DEGs_complete_info_simplified %>% dplyr::filter(is_als != "no") %>% pull(gene)

is_CAG <- mat_scaled_genes_sorted %in% CAGs_degs
is_other_addiction <- mat_scaled_genes_sorted %in% other_addiction_degs
is_psycho <- mat_scaled_genes_sorted %in% psycho_degs  
is_neuropeptide <- mat_scaled_genes_sorted %in% neuropeptide_degs  
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
  Heatmap(is_psycho + 0, name = "psychotic", col = c("0" = "white", "1" = "black"), 
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

#Heatmaps is better to save in tiff format. 
#Calculate the pixel dimensions based on the desired size and DPI
width_pixels <- 16 * 300  # Width in inches * DPI
height_pixels <- 15 * 300   # Height in inches * DPI
 
# tiff("GEX_cocaine_vs_saline_heatmaps.tiff", width = width_pixels, height = height_pixels, res = 300)
# draw(cocaine_vs_saline_heatmap)
# dev.off()


# ========== 4. Classifier Distribution Barplots ==========
# prepare data
ht_row_order_withFeatures <- cbind(ht_row_order, gene = mat_scaled_genes_sorted, isCAG = is_CAG, isPyscho = is_psycho) %>%
  dplyr::mutate(
    direction = ifelse(startsWith(as.character(groupID), "saline"), "down", "up"),
    timepoint = ifelse(startsWith(as.character(groupID), "saline"), sub("^saline_", "", as.character(groupID)), as.character(groupID)))

unique_genes_df <- unique(DEGs_complete_info_simplified[, c("gene", "is_CAG", "is_addiction", "is_psychiatric", "is_neuropeptide", "is_alzheimer", "is_parkinson", "is_als")])

DEGs_complete_info_Classifiers <- ht_row_order_withFeatures[, c("index", "groupID", "gene", "direction", "timepoint")] %>% 
  left_join(unique_genes_df, by = "gene")

data <- DEGs_complete_info_Classifiers
data$direction <- factor(data$direction, levels = c("up", "down"))

upreg_genes <- data %>% dplyr::filter(direction == "up") %>% pull(gene)
downreg_genes <- data %>% dplyr::filter(direction == "down") %>% pull(gene)
no_degs <- unique_genes_df %>% dplyr::filter(gene %notin% c(upreg_genes, downreg_genes))%>% pull(gene)

#Calculate percentages for each category and direction
#First, calculate the total number of upregulated and downregulated genes
total_DEGs <- data %>%
  group_by(direction) %>%
  dplyr::summarise(total = n(), .groups = "drop")

#Then calculate the count of 'yes' statuses for each classifier within each direction
percentages <- data %>%
  pivot_longer(cols = c("gene", "is_CAG", "is_addiction", "is_psychiatric", "is_neuropeptide", "is_alzheimer", "is_parkinson", "is_als"),
               names_to = "classifier", values_to = "status") %>%
  dplyr::filter(status != "no") %>%
  group_by(direction, classifier) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  left_join(total_DEGs, by = "direction") %>%
  mutate(percentage = count / total * 100) %>% 
  dplyr::filter(!is.na(percentage)) %>% 
  complete(direction, classifier, fill = list(count = 0, percentage = 0)) %>%   # Fill missing combinations
  mutate(total = if_else(is.na(total) & direction == "up" & classifier == "in_als", 
                         sum(total_DEGs$total[total_DEGs$direction == "up"]),  # Assuming total_DEGs has total for "up"
                         total),
         percentage = count / total * 100) %>% 
  dplyr::filter(!is.na(percentage), classifier !="gene")  

classifier_order <- c("gene", "is_CAG", "is_addiction", "is_psychiatric", "is_neuropeptide", "is_alzheimer", "is_parkinson", "is_als")
percentages$classifier <- factor(percentages$classifier, levels = classifier_order)

max_percentage <- max(percentages$percentage)

#Percentage of DEGs Associated with Neurological and Addiction-Related Classifiers
percentages <- percentages %>% dplyr::filter(classifier != "is_neuropeptide")
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
# plot_name <- "chapter1_classifiers_upanddown_over_all_v3.pdf"
# pdf(plot_name, width = width_original, height =height_original)
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


# ========== 5. UpSet plots ==========
##Upset plot
##############
classifiers <-  c("is_CAG", "is_addiction", "is_psychiatric", "is_neuropeptide", "is_alzheimer", "is_parkinson", "is_als")

# Initialize lists
upregulated_upset <- list()
downregulated_upset <- list()

# Loop through classifiers to create lists
for (classifier in classifiers) {
  upregulated_upset[[classifier]] <- DEGs_complete_info_Classifiers %>%
    dplyr::filter(direction == "up", .data[[classifier]] != "no") %>%
    pull(gene)
  
  downregulated_upset[[classifier]] <- DEGs_complete_info_Classifiers %>%
    dplyr::filter(direction == "down", .data[[classifier]] != "no") %>%
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


# ========== 6. Significance tests ==========
# Calculated for all gene lists
gene_sets <- list(
  CAGs = cags_all,
  Psychiatric = psychiatric_all,
  Addiction = general_addiction_all,
  Alzheimer = alzheimer_all$mouse_gene_name,
  Parkinson = parkinson_all$mouse_gene_name,
  ALS = als_all$mouse_gene_name
)

test_enrichment <- function(query_genes, background_genes, category_genes, category_name, direction_name) {
  in_cat <- length(intersect(query_genes, category_genes))
  in_bg_cat <- length(intersect(background_genes, category_genes))
  
  non_cat <- length(query_genes) - in_cat
  non_bg_cat <- length(background_genes) - in_bg_cat
  
  tab <- matrix(
    c(in_cat, non_cat,
      in_bg_cat, non_bg_cat),
    nrow = 2, byrow = TRUE,
    dimnames = list(
      c(direction_name, "No-DEGs"),
      c(category_name, paste0("Non-", category_name))
    )
  )
  
  chi <- chisq.test(tab)
  
  list(
    table = tab,
    test = chi,
    counts = data.frame(
      category = category_name,
      direction = direction_name,
      in_category = in_cat,
      in_background_category = in_bg_cat,
      p_value = chi$p.value,
      stringsAsFactors = FALSE
    )
  )
}

results <- list()
summary_df <- list()

for (cat_name in names(gene_sets)) {
  cat_genes <- gene_sets[[cat_name]]
  
  res_up <- test_enrichment(upreg_genes, no_degs, cat_genes, cat_name, "Upregulated")
  res_down <- test_enrichment(downreg_genes, no_degs, cat_genes, cat_name, "Downregulated")
  
  results[[paste0(cat_name, "_up")]] <- res_up
  results[[paste0(cat_name, "_down")]] <- res_down
  
  summary_df[[paste0(cat_name, "_up")]] <- res_up$counts
  summary_df[[paste0(cat_name, "_down")]] <- res_down$counts
  
  print(res_up$table)
  print(res_up$test)
  
  print(res_down$table)
  print(res_down$test)
}

summary_df <- do.call(rbind, summary_df)
summary_df$label <- dplyr::case_when(
  summary_df$p_value < 0.001 ~ "p < 0.001 (***)",
  summary_df$p_value < 0.01  ~ paste0("p = ", format(round(summary_df$p_value, 3), nsmall = 3), " (**)"),
  summary_df$p_value < 0.05  ~ paste0("p = ", format(round(summary_df$p_value, 2), nsmall = 2), " (*)"),
  TRUE ~ paste0("p = ", format(round(summary_df$p_value, 2), nsmall = 2))
)

summary_df$label <- gsub("\\.", ",", summary_df$label)
summary_df


# ========== 7. UpSet plots for figure panels ==========
# For main fifure 
classifiers <- c("is_CAG", "is_addiction", "is_psychiatric")

# Initialize lists
upregulated_upset <- list()
downregulated_upset <- list()

# Loop through classifiers to create lists
for (classifier in classifiers) {
  upregulated_upset[[classifier]] <- DEGs_complete_info_Classifiers %>%
    dplyr::filter(direction == "up", .data[[classifier]] != "no") %>%
    pull(gene)
  
  downregulated_upset[[classifier]] <- DEGs_complete_info_Classifiers %>%
    dplyr::filter(direction == "down", .data[[classifier]] != "no") %>%
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

#View(as.data.frame(extract_comb(m_down, "001")))
extract_comb(m_down, "010")

desired_order <- comb_size(m_up)[c("111", "110", "100", "011", "010", "001")]
desired_order_down <- c(1,2,4,3,5,6)
downregulated_upset_plot <- UpSet(t(m_down), set_order =classifiers, comb_order = desired_order_down, 
                                  top_annotation = upset_top_annotation(t(m_down), add_numbers = TRUE),
                                  right_annotation = upset_right_annotation(t(m_down), add_numbers = TRUE))

downregulated_upset_plot
print(upregulated_upset_plot) + print(downregulated_upset_plot)

# 
# dev.size()
# width_original <- 11.75
# height_original <- 8.36
# 
# plot_name <- "chapter2_upANDdownregulated_upset_plot_classifiers.pdf"
# pdf(plot_name, width = width_original, height =height_original)
# print(upregulated_upset_plot) + print(downregulated_upset_plot)
# dev.off()

# ========== 8. Save session ==========
#load("260514_DEG_visualization.rds")
