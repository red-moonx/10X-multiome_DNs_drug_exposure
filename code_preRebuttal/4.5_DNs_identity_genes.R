# ===========================================
# Script Title: 4.5_DNs_identity_genes.R
# Author: Luna Zea Redondo
# Date: 2025-01-17
# Description:
#   Evaluates the expression of DN identity genes using orthogonal data
#   from Phillips et al. 2022.
#   Main steps:
#     - load DEG tables and DN marker list
#     - define ranked DN identity gene sets (top 25/50/100/200/300)
#     - generate MA plots with highlighted identity genes
#     - compare logFC distributions with Wilcoxon tests
#     - plot boxplots and median trends across contrasts
#     - merge with memory-class DEG annotations for interpretation
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
library(glue)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggExtra)

# ========== Project setup ==========
source("/fast/AG_Pombo/luna/2023_vta_multiome/2024/5_figures/pixels_to_cm.R")
setwd("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/10_250117_DNs_identity")


# ===========================================
# 1. Load inputs and define gene sets
# ===========================================
#DNs_identity <- c("Lmx1a", "En1", "Lmx1b", "Nr4a2", "Foxa2", "Sox2", "En2", "Msx1", "Sox6", "Th", "Slc18a1", "Sox1", "Drd2", "Cck", "Ddc", "Foxa1", "Otx2", "Kcnj6", "Pitx3", "Slc18a2", "Drd1", "Slc6a3")
DEG_complete_results <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002_muscat_VTAvsSN_all_subtypes/Analysis2/deg_results/231002_DEG_complete_results_Analysis2.tsv") %>% 
  dplyr::filter(cluster_id == "VTA") 

DN_markers <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/10_250117_DNs_identity/phillip.tsv") %>% 
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


DEG_identity_toPlot <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002_muscat_VTAvsSN_all_subtypes/Analysis2/deg_results/231002_DEG_complete_results_Analysis2.tsv") %>% 
  dplyr::filter(cluster_id == "VTA") 


# ===========================================
# 2. MA plots for DN identity genes
# ===========================================
#Establish all possible contrasts:
all_contrasts <- c("h1_cocaine-saline",
                   "h4_cocaine-saline",
                   "h8_cocaine-saline",
                   "h24_cocaine-saline",
                   "d14_cocaine-saline")

# Create a list to store all plots
MAplots <- list()
# Loop through each contrast
for (i in seq_along(all_contrasts)) {
  comparison <- all_contrasts[i]
  plot.data <- DEG_complete_results %>%
    filter(contrast == comparison) %>%
    mutate(to_label = ifelse(gene %in% c("Th", "En1", "Slc6a3", "Slc18a2"), "yes", "no"),
           in_top_100 = ifelse(gene %in% top_100_genes$gene, "yes", "no"))
  
  query <- str_split(comparison, "-")[[1]][1]
  control <- str_split(comparison, "-")[[1]][2]
  
  # Define color values for the plot based on top_100_genes
  volcano_colors <- c("gray90", "#D580BA")
  
  # Generate the MA plot
  MA_plot <- ggplot(plot.data, aes(x = logCPM, y = logFC, color = in_top_100)) +
    # Plot gray dots first
    geom_point(size = 1, alpha = 0.5) +
    # Plot pink dots on top
    geom_point(data = filter(plot.data, in_top_100 == "yes"), size = 1, alpha = 0.8) +  # Highlight pink dots
    scale_color_manual(values = volcano_colors, labels = c("No significant", "Top 100 Genes")) +
    theme_classic(base_size = 12) +
    geom_hline(yintercept = 0, colour = "goldenrod") +
    ylim(-6, 8) +
    theme(legend.position = "none", legend.title = element_blank(), text = element_text(size = 15),
          axis.title = element_blank(), plot.title = element_blank()) +
    geom_text_repel(data = subset(plot.data, plot.data$to_label == "yes"), 
                    aes(label = gene), 
                    size = 3,
                    min.segment.length = 0, 
                    seed = 42, 
                    box.padding = 0.5,
                    max.overlaps = Inf,
                    nudge_x = .15,
                    nudge_y = .5,
                    color = "black"
    )
  
  # Add the marginal density plot using ggmarginal
  combined_plot <- ggExtra::ggMarginal(
    MA_plot,
    margins = "y", 
    type = "density", 
    groupColour = TRUE, groupFill = TRUE)
  
  # Store the combined plot in the list
  MAplots[[paste("MAplot", i, sep = "_")]] <- combined_plot
}

# Combine specific plots into a single row
MA_row1 <- wrap_plots(
  MAplots$MAplot_1, 
  MAplots$MAplot_2, 
  MAplots$MAplot_3, 
  MAplots$MAplot_4, 
  MAplots$MAplot_5, 
  nrow = 1
)

MA_row1

#Save
# dev.size()
# width_original =15.108333
# height_original= 2.566667
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter4_MAplots_DN_identity_genes_top100"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# dev.size()
# pdf(file_name, width = width_original*1.5, height =height_original*1.5)
# MA_row1
# dev.off()


# ===========================================
# 3. Wilcoxon tests for DN identity gene sets
# ===========================================
DEGs_for_test <- DEG_complete_results %>%
  dplyr::select(gene, logFC, contrast)

# Initialize a results data frame
contrast = all_contrasts[1]

# Loop through each contrast
contrast_data <- DEGs_for_test[DEGs_for_test$contrast == contrast, ]

# Extract logFC for the 100 genes
lfc_100 <- contrast_data$logFC[contrast_data$gene %in% top_100_genes$gene]

# Extract logFC for all other genes
lfc_all <- contrast_data$logFC

# Perform the Wilcoxon test
#test_result <- ks.test(lfc_500, lfc_all)
test_result <- wilcox.test(lfc_500, lfc_all, alternative = "one.sided")
test_result


ggplot(results_df, aes(x = contrast, y = p_value)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(
    title = "P-values for Wilcoxon Test Across Contrasts",
    x = "Contrast",
    y = "P-value"
  ) +
  theme_minimal()

# ===========================================
# 4. Boxplots for top-ranked identity genes
# ===========================================

#Establish all possible contrasts:
all_contrasts <- c("h1_cocaine-saline",
                   "h4_cocaine-saline",
                   "h8_cocaine-saline",
                   "h24_cocaine-saline",
                   "d14_cocaine-saline")

# top_100_genes <- DN_markers %>%
#   arrange(p.adj, desc(log2FC)) %>%  # Sort by p.adj (increasing) and log2FC (decreasing)
#   slice_head(n = 200) 

#For figure 4
# Precompute numeric x-axis for boxplot
volcano_colors <- c("gray90", "#76A267")

# Make sure 'contrast' is ordered the way you want (2 rows, left to right)
plot.data$contrast <- factor(plot.data$contrast, levels = c(
  "h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline",
  "h24_cocaine-saline", "d14_cocaine-saline"
))

plot.data <- DEG_identity_toPlot %>%
  mutate(to_label = ifelse(gene %in% c("Th", "En1", "Slc6a3", "Slc18a2"), "yes", "no"),
         in_top_100 = ifelse(gene %in% top_100_genes$gene, "yes", "no"), 
         in_top_200 = ifelse(gene %in% top_200_genes$gene, "yes", "no"), 
         in_top_300 = ifelse(gene %in% top_300_genes$gene, "yes", "no")) %>% 
  dplyr::filter(contrast %in% all_contrasts) %>% 
  dplyr::mutate(contrast = factor(contrast, levels =all_contrasts))

# Plot
test <- ggplot(plot.data, aes(x = in_top_200, y = logFC, fill = in_top_200, colour = in_top_200)) +
  #  geom_flat_violin(position = position_nudge(x = 0.25), adjust = 2, trim = TRUE) +
  geom_point(position = position_jitter(width = 0.15), size = 0.25) +
  geom_boxplot(position = position_nudge(x = 0.25), outlier.shape = NA, alpha = 0.3, 
               width = 0.1, colour = "black", lwd = 0.4) +
  # Add median value as a label
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
               position = position_nudge(x = 0.25, y = 0.1), size = 2.5, vjust = 0, colour = "black") +
  labs(title = "", x = "", y = "") +
  #  scale_y_log10() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = volcano_colors) + 
  scale_color_manual(values = volcano_colors) +
  #  geom_text_repel(data = plot.data[plot.data$to_label == "yes", "gene"]) +
  theme_classic() + 
  theme(legend.position = "none") +
  facet_wrap(~ contrast, nrow = 1) 

test

width_original =4.675
# height_original= 6.475
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter4_DN_identity_genes_top500"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# dev.size()
# pdf(file_name, width = width_original*1.5, height =height_original*1.5)
# test
# dev.off()

# ===========================================
# 5. Median logFC trend plot across contrasts
# ===========================================
#More simplified version of the boxplot
# Plot
test <- ggplot(plot.data, aes(x = in_top_200, y = logFC, fill = in_top_200, colour = in_top_200)) +
  # Boxplot
  geom_boxplot(outlier.shape = NA, colour = "black", ) +
  # Median text
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
               position = position_nudge(x = 0.25, y = 0.1), size = 2.5, vjust = 0, colour = "black") +
  # Reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "#DAA628", linewidth = 1) +
  # Colors
  scale_fill_manual(values = volcano_colors) + 
  scale_color_manual(values = volcano_colors) +
  # Set y-axis display range
  coord_cartesian(ylim = c(-2, 3)) +
  labs(title = "", x = "", y = "") +
  theme_bw() + 
  theme(legend.position = "none") +
  facet_wrap(~ contrast, nrow = 1) 
test


# dev.size()
# width_original =11.383333
# height_original= 3.750
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter4_DN_identity_genes_top200_v2"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# 
# pdf(file_name, width = width_original, height =height_original)
# test
# dev.off()


#####################################################
##Significance test for 100, 200, 300 all samples
###################################################
# Select relevant columns
DEGs_for_test <- DEG_complete_results %>%
  dplyr::select(gene, logFC, contrast)

# Named list of gene sets
top_gene_sets <- list(
  top_25 = top_25_genes$gene,
  top_50 = top_50_genes$gene,
  top_100 = top_100_genes$gene,
  top_200 = top_200_genes$gene
)

# Significance formatter
format_significance <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "n.s",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "n.s"
  )
}

# Run Wilcoxon tests for all contrasts and gene sets
final_results <- purrr::map_df(all_contrasts[1:5], function(contrast) {
  contrast_data <- DEGs_for_test %>%
    dplyr::filter(contrast == !!contrast)
  
  purrr::map_df(names(top_gene_sets), function(set_name) {
    genes <- top_gene_sets[[set_name]]
    lfc_top <- contrast_data$logFC[contrast_data$gene %in% genes]
    lfc_rest <- contrast_data$logFC[!contrast_data$gene %in% genes]
    
    p <- if (length(lfc_top) > 0 && length(lfc_rest) > 0) {
      wilcox.test(lfc_top, lfc_rest, alternative = "less")$p.value
    } else {
      NA
    }
    
    dplyr::tibble(
      contrast = contrast,
      gene_set = set_name,
      p_value = p,
      significance = format_significance(p)
    )
  })
})

# View results
# final_results
# write.csv(final_results, "250818.DNsidentity.onesidedwilcoxon.csv")

# Run Wilcoxon tests for all contrasts and gene sets
final_results <- purrr::map_df(all_contrasts[1:5], function(cur_contrast) {
  contrast_data <- DEGs_for_test %>%
    dplyr::filter(.data$contrast == .env$cur_contrast)
  
  purrr::map_df(names(top_gene_sets), function(set_name) {
    genes    <- top_gene_sets[[set_name]]
    lfc_top  <- contrast_data$logFC[contrast_data$gene %in% genes]
    lfc_rest <- contrast_data$logFC[!(contrast_data$gene %in% genes)]
    
    # Remove NA values
    lfc_top  <- lfc_top[!is.na(lfc_top)]
    lfc_rest <- lfc_rest[!is.na(lfc_rest)]
    
    p <- if (length(lfc_top) > 0 && length(lfc_rest) > 0) {
      # Two-sided test for misregulation
      wilcox.test(lfc_top, lfc_rest, alternative = "two.sided")$p.value
    } else {
      NA_real_
    }
    
    dplyr::tibble(
      contrast     = cur_contrast,
      gene_set     = set_name,
      p_value      = p,
      significance = format_significance(p)
    )
  })
})

# View results
# final_results
# write.csv(final_results, "250818.DNsidentity.twosidedwilcoxon.csv")

# ===========================================
# 6. Merge with memory-class DEG annotations
# ===========================================
plot.data <- DEG_identity_toPlot %>%
  mutate(in_top_25 = ifelse(gene %in% top_25_genes$gene, "yes", "no"), 
         in_top_50 = ifelse(gene %in% top_50_genes$gene, "yes", "no"),
         in_top_100 = ifelse(gene %in% top_100_genes$gene, "yes", "no"), 
         in_top_200 = ifelse(gene %in% top_200_genes$gene, "yes", "no"), 
         in_top_300 = ifelse(gene %in% top_300_genes$gene, "yes", "no"),
         background = ifelse(!(gene %in% top_300_genes$gene), "yes", "no")) %>% 
  dplyr::filter(contrast %in% all_contrasts) %>% 
  dplyr::mutate(contrast = factor(contrast, levels =all_contrasts)) %>% 
  dplyr::select(gene, logFC, contrast, in_top_25, in_top_50, in_top_100, in_top_200, in_top_300,background)


# Step 1: Reshape and compute medians
median_data <- plot.data %>%
  pivot_longer(cols = starts_with("in_top_") | background, 
               names_to = "group", values_to = "is_yes") %>%
  dplyr::filter(is_yes == "yes") %>%
  dplyr::group_by(group, contrast) %>%
  dplyr::summarise(median_logFC = median(logFC, na.rm = TRUE), .groups = "drop")


# Step 2: Line plot
ggplot(median_data, aes(x = contrast, y = median_logFC, color = group, group = group)) +
  geom_line(size = 1) +
  geom_point() +
  theme_minimal() +
  labs(title = "Median logFC of 'Yes' Genes by Group and Contrast",
       x = "Contrast",
       y = "Median logFC",
       color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#Join DEG classification
memory_DEGs <- readRDS("../250205_DEG_complete_results_kmeans.corrected.4243.rds") %>% 
  dplyr::mutate(memory_class_dir = glue("{memory_class}_{direction}")) %>% 
  dplyr::select(gene, memory_class_dir, kmeans_cluster, is_CAG, is_addiction, is_psychiatric)

#Top 25
memory_DEGs_with_rank <- plot.data[, c(1, 4:9)] %>% distinct() %>% 
  left_join(memory_DEGs, by = "gene") %>% 
  dplyr::filter(in_top_25 == "yes")

table(memory_DEGs_with_rank$memory_class_dir) 
table(memory_DEGs_with_rank$kmeans_cluster) 
table(memory_DEGs_with_rank$is_CAG) 
