# ===========================================
# Script Title: 1.5 DEGs visualization: MA plots and Heatmap
# Author: Luna Zea Redondo
# Date: 2025-06-25
# Description:
#   This script generates visual figures for publication, including
#   MA plots, heatmaps, and barplots based on DEG results from 
#   Muscat analysis in 1.4. Additional overlays include gene classifications 
#   such as CAGs, psychiatric risk genes, and neurodegeneration genes.
# ===========================================

# ========== Set Environment ==========
rm(list = ls(all.names = TRUE))
gc()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

# ========== Load Required Libraries ==========
library(muscat)
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(glue)
library(patchwork)
library(readr)
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(UpSetR)
library(scales)
library(limma)
library(reshape2)
library(colorspace)

setwd("/fast/AG_Pombo/luna/2023_vta_multiome/2024/5_figures/240812_DEGs_MAandHeatmap")

# ========== Load Data ==========
#Load environment from previous analysis (1.4)
load("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002.RNAonly.excluding_putative_SN_cells.rds")

#DEGs per contrast
DEG_complete_results <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002_muscat_VTAvsSN_all_subtypes/Analysis2/deg_results/231002_DEG_complete_results_Analysis2.tsv") %>% 
  dplyr::filter(cluster_id == "VTA")
genes <- DEG_complete_results %>% dplyr::select(gene) %>% distinct()

# For convenience, the following table compiles all VTA-DN DEGs along with additional orthogonal annotations, 
# including disease-related classifiers and GAM-derived features.
# It also includes time series classifications from the k-means analysis in section 2.
# Basically, it's a time-traveling table.

DEGs_complete_info <- readRDS("../results/250205_DEG_complete_results_kmeans.corrected.4243.rds")

# The original disease-related gene sets and mouse–human name conversion table 
# are in "data/gene_lists" so you can rebuild the DEG table from scratch


# ========== 2. Generate MA Plots ==========
#Manual selection of genes to highlight
downreg_genes <- c("Fancg", "Map2k1", "Ntsr1", "Drd2", "Grik1")
upreg_genes <- c("Dtx3", "Bdnf", "Cartpt", "Vip", "Penk", "Ucn", "Ptprt", "Npy")
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
    ylim(-6, 9) +
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
# tiff("/fast/AG_Pombo/luna/2025_pdf_files/chapter1_MAplots_DEGs_v2.tif", width = 19.156250, height = 8.291667, res = 300)
# MA_row1/MA_row2
# dev.off()

# ========== 3. Generate Heatmap ==========
#Extablisg contrasts
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


#2. Get the TRUE/FALSE for each gene in each category:
CAGs_degs <-  DEGs_complete_info %>% filter(is_CAG != "no") %>% pull(gene)
other_addiction_degs <-  DEGs_complete_info %>% filter(is_addiction != "no") %>% pull(gene)
psycho_degs <- DEGs_complete_info %>% filter(is_psychiatric != "no") %>% pull(gene)
alzheimer_degs <- DEGs_complete_info %>% filter(is_alzheimer != "no") %>% pull(gene)
parkinson_degs <- DEGs_complete_info %>% filter(is_parkinson != "no") %>% pull(gene)
als_degs <- DEGs_complete_info %>% filter(is_als != "no") %>% pull(gene)

is_CAG <- mat_scaled_genes_sorted %in% CAGs_degs
is_other_addiction <- mat_scaled_genes_sorted %in% other_addiction_degs
is_psycho <- mat_scaled_genes_sorted %in% psycho_degs  
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
# Define the desired reverse order for time points
timepoint_order <- c("d14_cocaine", "h24_cocaine", "h8_cocaine", "h4_cocaine", "h1_cocaine")

# Calculate the percentage of "is_CAG" genes for each time point and direction
data_summary <- ht_row_order_withFeatures %>%
  group_by(timepoint, direction) %>%
  dplyr::summarise(total_genes = n(),
                   isCAG_genes = sum(isCAG == TRUE)) %>%
  mutate(percentage = round((isCAG_genes / total_genes) * 100, 1))

# Ensure the data includes all relevant time points and is ordered correctly
data_summary$timepoint <- factor(data_summary$timepoint, levels = timepoint_order)

# Create the bidirectional bar plot with percentage labels
plota <- ggplot(data_summary, aes(x = factor(timepoint, levels = timepoint_order), y = ifelse(direction == "up", percentage, -percentage), fill = direction)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = abs(percentage)), 
            position = position_stack(vjust = 0.5), 
            color = "black") +
  scale_y_continuous(labels = abs, limits = c(-5, 5)) +  # Set x-axis limits to -5 and +5
  coord_flip() +
  labs(title = "",
       x = "Timepoint",
       y = "Percentage of Cocaine Assiciated Genes (CAGs)") +
  theme_minimal() +
  scale_fill_manual(values = c("up" = "#EC5D5B", "down" = "#89B8D2"))


# Calculate the percentage of "isPyscho" genes for each timepoint and direction
data_summary_psycho <- ht_row_order_withFeatures %>%
  group_by(timepoint, direction) %>%
  dplyr::summarise(total_genes = n(),
                   isPyscho_genes = sum(isPyscho == TRUE)) %>%
  mutate(percentage = round((isPyscho_genes / total_genes) * 100, 1))

# Ensure the data includes all relevant time points and is ordered correctly
data_summary_psycho$timepoint <- factor(data_summary_psycho$timepoint, levels = timepoint_order)

# Create the bidirectional bar plot with percentage labels for isPyscho genes
plotb <- ggplot(data_summary_psycho, aes(x = factor(timepoint, levels = timepoint_order), y = ifelse(direction == "up", percentage, -percentage), fill = direction)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = abs(percentage)), 
            position = position_stack(vjust = 0.5), 
            color = "black") +
  scale_y_continuous(labels = abs, limits = c(-12, 12)) +  # Set x-axis limits to -5 and +5
  coord_flip() +
  labs(title = "",
       x = "Time point",
       y = "Percentage of psychiatric disorder genes") +
  theme_minimal() +
  scale_fill_manual(values = c("up" = "#EC5D5B", "down" = "#89B8D2"))

#Both together
plota + plotb



# Percentages of DEGs in addiction and neurological disorders

DEGs_complete_info_Classifiers <- ht_row_order_withFeatures[, c("index", "groupID", "gene", "direction", "timepoint")] %>% 
  left_join(DEGs_complete_info[, c("gene", "is_CAG", "is_addiction", "is_psychiatric", "is_alzheimer", "is_parkinson", "is_als")])

data <- DEGs_complete_info_Classifiers
data$direction <- factor(data$direction, levels = c("up", "down"))

#Calculate percentages for each category and direction
#First, calculate the total number of upregulated and downregulated genes
total_DEGs <- data %>%
  group_by(direction) %>%
  dplyr::summarise(total = n(), .groups = "drop")

#Second, calculate the count of 'yes' statuses for each classifier within each direction
percentages <- data %>%
  pivot_longer(cols = c("is_CAG", "is_addiction", "is_psychiatric", "is_alzheimer", "is_parkinson", "is_als"),
               names_to = "classifier", values_to = "status") %>%
  filter(status != "no") %>%
  group_by(direction, classifier) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  left_join(total_DEGs, by = "direction") %>%
  mutate(percentage = count / total * 100) %>%  # Calculate percentage based on direction-specific totals
  complete(direction, classifier, fill = list(count = 0, percentage = 0)) %>%   # Fill missing combinations
  mutate(total = if_else(is.na(total) & direction == "up" & classifier == "in_als", 
                         sum(total_DEGs$total[total_DEGs$direction == "up"]),  # Assuming total_DEGs has total for "up"
                         total),
         percentage = count / total * 100, 
         percentage= ifelse(count == 0, 0, percentage)) 
  
classifier_order <- c("is_CAG", "is_addiction", "is_psychiatric", "is_alzheimer", "is_parkinson", "is_als")
percentages$classifier <- factor(percentages$classifier, levels = classifier_order)

max_percentage <- max(percentages$percentage)


DEGs_in_classifiers <- ggplot(percentages, aes(x = direction, y = percentage, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), position = position_dodge(width = 0.7), vjust = -0.25) +
  facet_wrap(~ classifier, scales = "free_y") +
  labs(title = "Percentage of DEGs Associated with Neurological and Addiction-Related Classifiers",
       x = "Direction",
       y = "Percentage (%)") +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(0, max_percentage), breaks = seq(0, max_percentage, 5)) +
  scale_fill_manual(values = c("up" = "#EC5D5B", "down" = "#89B8D2")) +
  theme(legend.position = "none")



# ========== 5. Upset Plot ==========
classifiers <- c("is_CAG", "is_addiction", "is_psychiatric", "is_alzheimer", "is_parkinson", "is_als")

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


# Get gene names for specific combinations:
# comb_size(m_down)
# extract_comb(m_down, "101110")


#Simplified upset plot for main figure. 

classifiers <- classifiers[1:3]

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

# Select preferred order
desired_order <- comb_size(m_up)[c("101", "100", "011", "010", "001")]
desired_order_up <- c(1,3,2,4,5)

upregulated_upset_plot <- UpSet(t(m_up), set_order =classifiers, comb_order = desired_order_up, 
                                top_annotation = upset_top_annotation(t(m_up), add_numbers = TRUE),
                                right_annotation = upset_right_annotation(t(m_up), add_numbers = TRUE)
)

downregulated_upset_plot <- UpSet(t(m_down), set_order =classifiers, comb_order = desired_order_up, 
                                  top_annotation = upset_top_annotation(t(m_down), add_numbers = TRUE),
                                  right_annotation = upset_right_annotation(t(m_down), add_numbers = TRUE))
