# ===========================================
# Script Title: 5.3 DN identity - heatmaps
# Author: Luna Zea Redondo
# Date: 2026-05-26
# Description:
#   Prepares manuscript ready figures for DN identity genes at different resolutions,
#   with classifiers.
# ===========================================


#251114 - Heatmap with DN identity genes
# Create heatmap s in fig2, with identifiers: CAGs, Addiction, Psychitric, FOXA2 targets, Polycomb targets (midbrain)
# Orthogonal data: Phillips et al 2022 
library(biomaRt)
library(dplyr)
library(HGNChelper)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(patchwork)
library(PupillometryR)
library(colorspace)
library(Seurat)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)
library(cluster)
library(glue)
library(tibble)


condition_colors <- c("black", "#9F2365", "#617641", "#C48208", "#326186", "#AE430A", "#564686")
condition_names <- c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")
names(condition_colors) <- condition_names


setwd("/fast/AG_Pombo/luna/2026_rebuttal/15_DN_identity")
mouse_human_conversion <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/mouse_human_conversion.rds")

# 1. Prepare the identifiers:
# --------------------------
#A) CAGs, Addiction and Psychiatric genes are obtained from "250127_classifiers_statistics_tests.R"
cags_all <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/cags_all.rds")
addiction_all <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/general_addiction_all.rds")
psychiatric_all <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/psychiatric_all.rds")

#B) Load FOXA2 targets 
foxa2_human <- read_csv("/fast/AG_Pombo/luna/2026_rebuttal/FOXA_targets_human.csv") #2936
foxa2_mouse <- foxa2_human %>%
  dplyr::left_join(mouse_human_conversion, 
                   by = c("Gene" = "human_gene_name")) %>%
  dplyr::select(mouse_gene_name, everything()) %>%
  dplyr::filter(!is.na(mouse_gene_name))

# Optional: rename and tidy up
foxa2_mouse <- foxa2_mouse %>%
  dplyr::rename(Mouse_Gene = mouse_gene_name,
                Human_Gene = Gene) #2778

#C) Polycomb targets 
toskas22 <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/toskas_table_s1_chromatinStates.csv") %>% 
  dplyr::select(Gene, dat_4mo_wt_state)

genes_K27 <- toskas22 %>%
  filter(grepl("K27", dat_4mo_wt_state)) %>%
  pull(Gene)


#2. Load DN identity genes
# --------------------------
DEG_complete_results <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/3_muscat/deg_results/260513_DEG_complete_results_Analysis1.tsv") 

DN_markers <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/phillip.tsv") %>% 
  dplyr::select(gene, p.adj, log2FC) %>%
  dplyr::filter(gene %in% DEG_complete_results$gene)

top_25_genes <- DN_markers %>%
  arrange(p.adj, desc(log2FC)) %>%  # Sort by p.adj (increasing) and log2FC (decreasing)
  slice_head(n = 25) 

top_50_genes <- DN_markers %>%
  arrange(p.adj, desc(log2FC)) %>%  # Sort by p.adj (increasing) and log2FC (decreasing)
  slice_head(n = 50) 

top_100_genes <- DN_markers %>%
  arrange(p.adj, desc(log2FC)) %>%  # Sort by p.adj (increasing) and log2FC (decreasing)
  slice_head(n = 100) 

top_200_genes <- DN_markers %>%
  arrange(p.adj, desc(log2FC)) %>%  # Sort by p.adj (increasing) and log2FC (decreasing)
  slice_head(n = 200) 

top_300_genes <- DN_markers %>%
  arrange(p.adj, desc(log2FC)) %>%  # Sort by p.adj (increasing) and log2FC (decreasing)
  slice_head(n = 300)

which_top <- top_25_genes

#save table
DN_markers_v2 <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/10_250117_DNs_identity/phillip.tsv") %>%
  dplyr::select(gene, p.adj, log2FC) %>%
  filter(gene %in% DEG_complete_results$gene) %>%
  arrange(p.adj, desc(log2FC)) %>%
  mutate(rank = row_number(),
         top_category = case_when(
           rank <= 25 ~ "top25",
           rank <= 50 ~ "top50",
           rank <= 100 ~ "top100",
           rank <= 200 ~ "top200",
           rank <= 300 ~ "top300",
           TRUE ~ "other"
         ))

#saveRDS(DN_markers_v2, "DN_markers_v2.rds")

#3. Prepare heatmap
# --------------------------
saline_contrasts <-c("h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline", "h24_cocaine-saline", "d14_cocaine-saline")
time_contrasts <- c("h1_cocaine-saline", "h4_cocaine-h1_cocaine", "h8_cocaine-h4_cocaine", "h24_cocaine-h8_cocaine", "d14_cocaine-h24_cocaine")

samples_totest <- c("h4_saline_R1", "h8_saline_R1", "h24_saline_R1", "d14_saline_R1", #salines
                    "h1_cocaine_R1", "h1_cocaine_R2", "h4_cocaine_R1", "h4_cocaine_R2", "h8_cocaine_R1", "h8_cocaine_R2", #ETP
                    "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3") #LTP
condition <- c("saline", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")


seu.VTA_DNs #1520 cells DNs from VTA
seu.temporal <- seu.VTA_DNs

#Calculate average expression and format columns
Idents(seu.temporal) <- "orig.ident"
seu.temporal.avgexp = as.data.frame(AverageExpression(seu.temporal, group.by = 'orig.ident')) 
names(seu.temporal.avgexp) <- substring(names(seu.temporal.avgexp), 5)

# DEG_perGroup_results <- DEG_complete_results %>% filter(cluster_id == "VTA", contrast %in% saline_contrasts)
DEG_perGroup_results <- DEG_complete_results

#top X DN genes
heatmap_genes <- DEG_complete_results %>%
  dplyr::filter(gene %in% which_top$gene) %>% 
  dplyr::select("gene") %>% distinct() %>% pull()



#Calculate average expression and format columns
Idents(seu.temporal) <- "orig.ident"
seu.temporal.avgexp = as.data.frame(AverageExpression(seu.temporal, group.by = 'orig.ident')) 
names(seu.temporal.avgexp) <- substring(names(seu.temporal.avgexp), 5)

DEG_perGroup_results <- DEG_complete_results %>% dplyr::filter(contrast %in% saline_contrasts) 

#Extract DEGs from the big table and create metadata
DEG_withIDs <- DEG_perGroup_results %>% 
  dplyr::filter(gene %in% which_top$gene) %>% 
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


#Create heatmap
cocaine_vs_saline_heatmap <- 
  Heatmap(mat_scaled, split = cluster_rows_vector, cluster_rows = TRUE, column_split = cluster_cols_vector, 
          cluster_row_slices = FALSE, cluster_column_slices = FALSE,  
          top_annotation = ha, 
          right_annotation = row_ha,
          name = "Z-score", row_title = NULL, show_row_names = TRUE,  use_raster = FALSE, 
          col = heatmap_color, 
          show_row_dend = TRUE, show_column_dend = FALSE, 
          heatmap_legend_param = list(title = "Z-score", direction = "horizontal"), 
          show_column_names = FALSE)

ht_row_order <- do.call(rbind, lapply(names(row_order(cocaine_vs_saline_heatmap)), function(group) {
  data.frame(index = row_order(cocaine_vs_saline_heatmap)[[group]], groupID = group)
}))

ht_row_order <- as.data.frame(ht_row_order)
mat_scaled_genes_sorted <- rownames(mat_scaled)[ht_row_order$index]

selected_genes_to_plot <- which_top$gene
indices <- which(mat_scaled_genes_sorted %in% selected_genes_to_plot)

#Get true/false for this cateogires
CAG <- mat_scaled_genes_sorted %in% cags_all
addiction <- mat_scaled_genes_sorted %in% addiction_all
psychiatric <- mat_scaled_genes_sorted %in% psychiatric_all  
foxa2_targets <- mat_scaled_genes_sorted %in% foxa2_mouse$Mouse_Gene
polycomb_targets <- mat_scaled_genes_sorted %in% genes_K27


#Plot everything
ht_opt$ROW_ANNO_PADDING = unit(0.5, "cm")

row_ha = rowAnnotation(Group = gene_groupings, show_legend = FALSE, 
                       col = list(Group = groupings_cols, show_annotation_name = FALSE, simple_anno_size = unit(8, "mm")), 
                       foo = anno_mark(at = indices, 
                                       labels = mat_scaled_genes_sorted[indices]))
cocaine_vs_saline_heatmap <- 
  Heatmap(mat_scaled, split = cluster_rows_vector, cluster_rows = TRUE, column_split = cluster_cols_vector, 
          cluster_row_slices = FALSE, cluster_column_slices = FALSE,  
          top_annotation = ha, 
          right_annotation = row_ha,
          name = "Z-score", row_title = NULL, show_row_names = FALSE,  use_raster = FALSE, 
          col = heatmap_color, 
          show_row_dend = TRUE, show_column_dend = TRUE, 
          heatmap_legend_param = list(title = "Z-score", direction = "horizontal"), 
          show_column_names = TRUE, 
          show_heatmap_legend = FALSE) +
  Heatmap(CAG + 0, name = "CAGs", col = c("0" = "white", "1" = "black"), 
          show_heatmap_legend = FALSE, width = unit(10, "mm")) + 
  Heatmap(addiction + 0, name = "other addiction", col = c("0" = "white", "1" = "black"), 
          show_heatmap_legend = FALSE, width = unit(10, "mm")) +
  Heatmap(psychiatric + 0, name = "psychotic", col = c("0" = "white", "1" = "black"), 
          show_heatmap_legend = FALSE, width = unit(10, "mm")) +
  Heatmap(foxa2_targets + 0, name = "Foxa2_targets", col = c("0" = "white", "1" = "black"), 
          show_heatmap_legend = FALSE, width = unit(10, "mm")) +
  Heatmap(polycomb_targets + 0, name = "Polycomb_targets", col = c("0" = "white", "1" = "black"), 
          show_heatmap_legend = FALSE, width = unit(10, "mm"))
cocaine_vs_saline_heatmap


# dev.size()
# width_original =11.46875
# height_original= 12.12500
# 
# pdf_dir <- "heatmap_plots/"
# plot_name <- "top25_DNmarkers_heatmap_version1"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf")
# 
# dev.size()
# pdf(file_name, width = width_original, height =height_original)
# cocaine_vs_saline_heatmap
# dev.off()




#Second version: hierarchical clustering
#==========================================
# keep your samples ordered as given:
mat_scaled <- mat_scaled[, sample_order]   # ensure columns follow this order

# compute row distances (1 - Pearson correlation)
R <- cor(t(mat_scaled), method = "pearson", use = "pairwise.complete.obs")
d <- as.dist(1 - R)

# hierarchical clustering on rows
hc <- hclust(d, method = "average")

# visualize the dendrogram
plot(as.dendrogram(hc), main = "Row dendrogram", ylab = "1 - Pearson correlation")

# find approximate best k by silhouette (k = 2..10)
best_k <- NA
best_avg <- -Inf

#Less than 10 clusters
for (k in 2:9) {
  cl <- cutree(hc, k = k)
  sil <- silhouette(cl, d)
  avg_sil <- summary(sil)$avg.width
  cat("k =", k, " average silhouette width =", round(avg_sil, 3), "\n")
  
  if (avg_sil > best_avg) {
    best_avg <- avg_sil
    best_k <- k
  }
}

cat("\nBest k =", best_k, "with average silhouette width =", round(best_avg, 3), "\n")

# optional: visualize the silhouette for the best k
plot(silhouette(cutree(hc, k = best_k), d),
     main = paste("Silhouette plot for best k =", best_k))

# optional: heatmap with fixed column order
library(pheatmap)
pheatmap(mat_scaled,
         cluster_rows = as.hclust(hc),
         cluster_cols = FALSE,        # keep your given order
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = paste("Heatmap (best k =", best_k, ")"))









# rows = features to cluster; columns = samples
# If mat_scaled can contain NAs, use pairwise.complete.obs
R <- cor(t(mat_scaled), method = "pearson", use = "pairwise.complete.obs")
d <- as.dist(1 - R)                # correlation-based dissimilarity in [0, 2]
hc <- hclust(d, method = "average")
plot(as.dendrogram(hc))

n <- nrow(mat_scaled)
k_min <- 2
k_max <- min(10, n - 1)            # silhouette requires 2..n-1
best_k   <- NA_integer_
best_avg <- -Inf

for (k in k_min:k_max) {
  cl <- cutree(hc, k = k)
  sil <- silhouette(cl, d)         # from cluster::silhouette
  avg <- summary(sil)$avg.width
  cat("k =", k, " average silhouette width =", avg, "\n")
  if (avg > best_avg) {
    best_avg <- avg
    best_k <- k
  }
}

cat("Best k =", best_k, "with average silhouette width =", best_avg, "\n")

#visualize silhouette for the chosen k
#plot(silhouette(cutree(hc, k = best_k), d), main = paste("Silhouette, k =", best_k))

#Now use this k to make the cluster

# Step 1: Cut the dendrogram into 5 clusters
k <- best_k
cl <- cutree(hc, k)

# Step 2: Draw the heatmap
Heatmap(
  mat_scaled,
  name = "z-score",
  cluster_rows = hc,          # use your hierarchical clustering result
  row_split = k,             # color-separate rows by cluster membership
  clustering_distance_rows = "pearson",
  clustering_method_rows = "average",
  show_row_dend = TRUE,
  cluster_columns = TRUE,     # also cluster columns (optional)
  show_row_names = FALSE,     # cleaner view for many genes
  column_title = "Samples",
  row_title = "Genes"
)







# dev.size()
# width_original =11.46875
# height_original= 12.12500
# 
# pdf_dir <- "heatmap_plots/"
# plot_name <- "top25_DNmarkers_heatmap_version2"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf")
# 
# 
# pdf(file_name, width = width_original, height =height_original)
# draw(ht_main + ht_tracks,
#      heatmap_legend_side = "bottom",
#      annotation_legend_side = "bottom")
# dev.off()
# 





##################################
# Update: average per time point. 
###################################
which_top <- top_25_genes
which_top <- top_100_genes

# version 1:
###################
DEG_perGroup_results <- DEG_complete_results 

#top X DN genes
heatmap_genes <- DEG_complete_results %>%
  dplyr::filter(gene %in% which_top$gene) %>%
  dplyr::select("gene") %>% distinct() %>% pull()

seu.temporal <- seu.VTA_DNs

#Calculate average expression and format columns
Idents(seu.temporal) <- "simpleIdent"
seu.temporal.avgexp = as.data.frame(AverageExpression(seu.temporal, group.by = 'simpleIdent')) 
names(seu.temporal.avgexp) <- substring(names(seu.temporal.avgexp), 5)

DEG_perGroup_results <- DEG_complete_results %>% 
  dplyr::filter(contrast %in% saline_contrasts) 

#Extract DEGs from the big table and create metadata
DEG_withIDs <- DEG_perGroup_results %>% 
  dplyr::filter(gene %in% which_top$gene) %>% 
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
  #dplyr::select(all_of(samples_totest), gene, groupID) %>%
  dplyr::mutate(groupID = factor(groupID, levels = desired_order)) %>%
  dplyr::arrange(groupID) %>% 
  column_to_rownames("gene")

sample_order <- c("saline", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")
cluster_cols_vector <- factor(sample_order, levels = c("saline", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine"))
cluster_rows_vector <- seu.DEGs.avgexp.heatmap$groupID

#Prepare matrix for heatmap (in the order I want to be plotted)
mat <- as.matrix(seu.DEGs.avgexp.heatmap[, sample_order])
mat_scaled = t(scale(t(mat)))[, sample_order]


#Inciso: is data normally distributed
# Flatten all Z-scores into one numeric vector
vals <- as.vector(mat_scaled)

# Histogram
hist(vals, breaks = 50, col = "skyblue", main = "Distribution of Z-scores",
     xlab = "Z-score", border = "white")

# Density plot (smoother)
plot(density(vals, na.rm = TRUE),
     main = "Density of Z-scores across all genes/samples",
     xlab = "Z-score")

# QQ plot
qqnorm(vals, main = "QQ Plot of Z-scores")
qqline(vals, col = "red")


sample_colors <- condition_colors[cluster_cols_vector]

gene_groupings <- seu.DEGs.avgexp.heatmap$groupID
groupings_cols <- c("saline_h1_cocaine" = lighten("#617641", 0.5), "saline_h4_cocaine" = lighten("#C48208", 0.5), "saline_h8_cocaine" = lighten("#326186", 0.5), "saline_h24_cocaine" = lighten("#AE430A", 0.5), "saline_d14_cocaine" = lighten("#564686", 0.5),
                    "h1_cocaine" = "#617641", "h4_cocaine" = "#C48208", "h8_cocaine" = "#326186", "h24_cocaine" = "#AE430A", "d14_cocaine" = "#564686")


ha = HeatmapAnnotation(Condition = cluster_cols_vector, show_legend = FALSE, 
                       col = list(Condition = condition_colors), show_annotation_name = FALSE)

row_ha = rowAnnotation(Group = gene_groupings, show_legend = FALSE, 
                       col = list(Group = groupings_cols, show_annotation_name = FALSE, simple_anno_size = unit(8, "mm")))
heatmap_color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)  


#Create heatmap
cocaine_vs_saline_heatmap <- 
  Heatmap(mat_scaled, split = cluster_rows_vector, cluster_rows = TRUE, column_split = cluster_cols_vector, 
          cluster_row_slices = FALSE, cluster_column_slices = FALSE,  
          top_annotation = ha, 
          right_annotation = row_ha,
          name = "Z-score", row_title = NULL, show_row_names = TRUE,  use_raster = FALSE, 
          col = heatmap_color, 
          show_row_dend = TRUE, show_column_dend = FALSE, 
          heatmap_legend_param = list(title = "Z-score", direction = "horizontal"), 
          show_column_names = FALSE)

ht_row_order <- do.call(rbind, lapply(names(row_order(cocaine_vs_saline_heatmap)), function(group) {
  data.frame(index = row_order(cocaine_vs_saline_heatmap)[[group]], groupID = group)
}))

ht_row_order <- as.data.frame(ht_row_order)
mat_scaled_genes_sorted <- rownames(mat_scaled)[ht_row_order$index]

selected_genes_to_plot <- which_top$gene
indices <- which(mat_scaled_genes_sorted %in% selected_genes_to_plot)

#Get true/false for this cateogires
CAG <- mat_scaled_genes_sorted %in% cags_all
addiction <- mat_scaled_genes_sorted %in% addiction_all
psychiatric <- mat_scaled_genes_sorted %in% psychiatric_all  
foxa2_targets <- mat_scaled_genes_sorted %in% foxa2_mouse$Mouse_Gene
polycomb_targets <- mat_scaled_genes_sorted %in% genes_K27

## Version 2: compact + fixed annotation/gene-name alignment

R <- cor(t(mat_scaled), method = "pearson", use = "pairwise.complete.obs")
d <- as.dist(1 - R); hc <- hclust(d, method = "average"); row_dend <- as.dendrogram(hc)

best_k <- 10
group_levels <- c("saline","h1_cocaine","h4_cocaine","h8_cocaine","h24_cocaine","d14_cocaine")
cluster_cols_vector <- factor(group_levels, levels = group_levels)
heat_cols <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(101)

km_palette <- colorRampPalette(brewer.pal(12, "Set3"))(max(3, best_k))[seq_len(best_k)]
names(km_palette) <- paste0("k", seq_len(best_k))

group2 <- factor(cutree(hc, k = best_k)[rownames(mat_scaled)],
                 levels = seq_len(best_k), labels = paste0("k", seq_len(best_k)))

gene_groupings <- seu.DEGs.avgexp.heatmap[rownames(mat_scaled), "groupID"]
indices <- which(rownames(mat_scaled) %in% which_top$gene)

make_binary_mat <- function(gene_set, label) {
  m <- matrix(as.integer(rownames(mat_scaled) %in% gene_set), ncol = 1,
              dimnames = list(rownames(mat_scaled), label))
  m
}

binary_sets <- list(
  CAGs = cags_all,
  `other addiction` = addiction_all,
  psychotic = psychiatric_all,
  Foxa2_targets = foxa2_mouse$Mouse_Gene,
  Polycomb_targets = genes_K27
)

binary_mats <- Map(make_binary_mat, binary_sets, names(binary_sets))

ha <- HeatmapAnnotation(Condition = cluster_cols_vector, show_legend = FALSE,
                        col = list(Condition = condition_colors),
                        show_annotation_name = FALSE)

row_ha <- rowAnnotation(
  Group = gene_groupings, Group2 = group2,
  col = list(Group = groupings_cols, Group2 = km_palette),
  show_legend = FALSE, show_annotation_name = FALSE,
  simple_anno_size = unit(8, "mm"),
  foo = anno_mark(at = indices, labels = rownames(mat_scaled)[indices])
)

ht_opt$ROW_ANNO_PADDING <- unit(0.5, "cm")

cocaine_vs_saline_heatmap <-
  Heatmap(mat_scaled, name = "Z-score", col = heat_cols, use_raster = TRUE,
          cluster_rows = row_dend, row_split = best_k,
          cluster_row_slices = FALSE, show_row_dend = TRUE,
          row_title = NULL, show_row_names = TRUE,
          column_split = cluster_cols_vector, cluster_columns = FALSE,
          cluster_column_slices = FALSE, show_column_dend = FALSE,
          show_column_names = FALSE, column_title_gp = gpar(fontsize = 12),
          column_title_rot = 45, top_annotation = ha, right_annotation = row_ha,
          heatmap_legend_param = list(title = "Z-score", direction = "horizontal"),
          show_heatmap_legend = FALSE) +
  Reduce(`+`, Map(function(m, nm)
    Heatmap(m, name = nm, col = c("0" = "white", "1" = "black"),
            show_heatmap_legend = FALSE, width = unit(10, "mm"),
            column_names_rot = 45, column_names_side = "bottom"),
    binary_mats, names(binary_mats)))

cocaine_vs_saline_heatmap


dev.size()
plot_name <- "groupAverage_top100_DNmarkers_heatmap_version2.pdf"
pdf(plot_name, width = dev.size()[1], height =dev.size()[2]*2.5)
cocaine_vs_saline_heatmap
dev.off()

#save clustering as table
# assuming 'group2' is a factor with cluster labels (k1..kK)
gene_cluster_df <- data.frame(
  gene     = rownames(mat_scaled),
  cluster  = as.character(group2),
  stringsAsFactors = FALSE
)

#save to CSV
write.csv(gene_cluster_df, "groupAverage_top100_DNmarkers_cluster_assignment.csv", row.names = FALSE)


#Code to find DEGs
significant_genes <- DEG_complete_results %>%
  dplyr::filter(significant != "No significant") %>% 
  dplyr::filter(gene %in% top_100_genes$gene) %>% 
  dplyr::select(gene) %>% pull() %>% unique()



#Plot Neurod6 expression
#########################

seu.VTA_DNs$simpleIdent <- factor(seu.VTA_DNs$simpleIdent,
                                  c("saline","h1_cocaine","h4_cocaine","h8_cocaine","h24_cocaine","d14_cocaine"))

neurod6_umaps<- FeaturePlot(seu.VTA_DNs, "Neurod6", split.by = "simpleIdent",
            cols = c("lightgrey", "red"), order = TRUE, combine = FALSE) |>
  wrap_plots(nrow = 2, ncol = 3) &
  NoAxes() & NoLegend() & theme_void()



# dev.size()
# 
# width_original =10.28
# height_original= 6.63
# 
# plot_name <- "neurod6_umaps_expression.pdf"
# 
# pdf(plot_name, width = width_original, height =height_original)
# neurod6_umaps
# dev.off()
