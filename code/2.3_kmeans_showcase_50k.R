# ===========================================
# Script Title: 2.3. 50kmeans clustering - showcase
# Author: Luna Zea Redondo
# Date: 2026-06-14
# Description:
#   This script analyzes DEGs clustered with k-means (50K input set), 
#   generates per-cluster gene lists, and visualizes the distribution 
#   of addiction- and disease-related gene annotations. 
#   It includes barplots, scatter plots, line plots, and prepares input
#   for EnrichR.
# ===========================================

# ========== Set Environment ==========
rm(list = ls(all.names = TRUE))
gc()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

# ========== Load Required Libraries ==========
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(glue)
library(readr)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)
library(colorspace)
library(circlize)
library(reshape2)
library(scales)
library(UpSetR)
library(purrr)


# ========== Set Working Directory ==========
setwd("/fast/AG_Pombo/luna/2026_rebuttal/5_kmeans_GEX")

# ========== 1. Load DEG Results ==========
# Load kmeans-clustered DEGs and add to master table
DEG_complete_results <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/4_DEGs_visualization/260513_DEGs_complete_info_simplified.tsv")
DEGs_kmeans <- read_csv("260514_dynamic_cluster_df_version2.csv")
colnames(DEGs_kmeans) <- c("kmeans_cluster", "gene")

# Create a directory to store the files
output_dir <- "kmeans_clusters"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Save per-cluster and full unique gene lists
DEGs_kmeans %>%
  split(.$kmeans_cluster) %>%
  imap(~ write_lines(.x$gene, file = file.path(output_dir, paste0("cluster_", .y, ".txt"))))

DEGs_kmeans %>%
  distinct(gene) %>%  # Select unique genes
  pull(gene) %>%  # Extract the column as a vector
  write_lines(file = file.path(output_dir, "all_unique_genes.txt"))  # Write to file



# ========== 2. Classifier Barplots ==========
DEGs_kmeans_classifier <- DEG_complete_results %>% 
  left_join(DEGs_kmeans, by = "gene") %>% 
  dplyr::select(gene, kmeans_cluster, is_CAG, is_addiction, is_psychiatric) %>% distinct() %>% na.omit()

# 2.1. Three classifiers on  3X1 grid. 
# Convert classifier columns to logical (TRUE/FALSE)
DEGs_kmeans_classifier_long <- DEGs_kmeans_classifier %>%
  pivot_longer(cols = starts_with("is_"), names_to = "classifier", values_to = "status") %>%
  mutate(status = ifelse(status == "yes", 1, 0))

# Compute proportion of "yes" per cluster and classifier
proportion_data <- DEGs_kmeans_classifier_long %>%
  dplyr::group_by(kmeans_cluster, classifier) %>%
  dplyr::summarise(proportion_yes = mean(status), .groups = "drop")

# 2.1. Three classifiers on  3X1 grid. 
ggplot(proportion_data, aes(x = factor(kmeans_cluster), y = proportion_yes, fill = classifier)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "K-means Cluster", y = "Proportion of 'yes'", fill = "Classifier") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ classifier, ncol = 1, scales = "free_y")

# 2.2. Three classifiers, proportions, on  1X3 grid. 
proportion_data <- proportion_data %>%
  mutate(classifier = factor(classifier, levels = c("is_CAG", "is_addiction", "is_psychiatric"))) %>% 
  dplyr::filter(!(is.na(kmeans_cluster)))

ggplot(proportion_data, aes(x = factor(kmeans_cluster), y = proportion_yes, fill = classifier)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "K-means Cluster", y = "Proportion of 'yes'", fill = "Classifier") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(
    values = c("is_CAG" = "#D21E76", "is_addiction" = "#ED691D", "is_psychiatric" = "#24A8BC")
  ) +
  facet_wrap(~ classifier, ncol = 1, scales = "free_y")  # Allows individual Y-axis scales

# 2.3. Three classifiers, horizontal (1X3 grid)
kmeans_classifier <- ggplot(proportion_data, aes(y = factor(kmeans_cluster), x = proportion_yes, fill = classifier)) +
  geom_bar(stat = "identity", position = "dodge") + 
  labs(y = "K-means Cluster", x = "Proportion of 'yes'", fill = "Classifier") +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +  # Keep y-axis readable
  scale_fill_manual(
    values = c("is_CAG" = "#D21E76", "is_addiction" = "#ED691D", "is_psychiatric" = "#24A8BC")
  ) +
  facet_grid(. ~ classifier, scales = "free_x")


#Save plot:
# dev.size()
# width_original = 7.666667
# height_original= 8.927083
# 
# 
# plot_name <- "chapter2_50kmeans_classifier.pdf"
# 
# pdf(file_name, width = width_original, height =height_original)
# kmeans_classifier
# dev.off()

# 2.4. One unique bar plot, check if any classifier is positive
# Create a new column that checks if any classifier is "yes"
DEGs_kmeans_classifier <- DEGs_kmeans_classifier %>%
  mutate(any_yes = ifelse(is_CAG == "yes" | is_addiction == "yes" | is_psychiatric == "yes", 1, 0))

# Compute proportion of "yes" per cluster
proportion_data <- DEGs_kmeans_classifier %>%
  dplyr::group_by(kmeans_cluster) %>%
  dplyr::summarise(proportion_any_yes = mean(any_yes), .groups = "drop")

# Plot
ggplot(proportion_data, aes(x = factor(kmeans_cluster), y = proportion_any_yes)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "K-means Cluster", y = "Proportion of Genes with Any 'Yes'") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 2.5. Only CAGs and psychiatric, including bidirectional plot (interesting, highlights 12 and 31)

DEGs_kmeans_classifier <- DEGs_kmeans_classifier %>% 
  dplyr::select(gene, kmeans_cluster, is_CAG, is_psychiatric)

# Convert classifier columns to logical (TRUE/FALSE)
DEGs_kmeans_classifier_long <- DEGs_kmeans_classifier %>%
  pivot_longer(cols = starts_with("is_"), names_to = "classifier", values_to = "status") %>%
  mutate(status = ifelse(status == "yes", 1, 0))

# Compute proportion of "yes" per cluster and classifier
proportion_data <- DEGs_kmeans_classifier_long %>%
  dplyr::group_by(kmeans_cluster, classifier) %>%
  dplyr::summarise(proportion_yes = mean(status), .groups = "drop")

# Plot version 1
ggplot(proportion_data, aes(x = factor(kmeans_cluster), y = proportion_yes, fill = classifier)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "K-means Cluster", y = "Proportion of 'yes'", fill = "Classifier") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot version 2 (bidireccional)
plot_data <- proportion_data %>%
  mutate(proportion_plot = ifelse(classifier == "is_psychiatric", -proportion_yes, proportion_yes))

ggplot(plot_data, aes(x = factor(kmeans_cluster), y = proportion_plot, fill = classifier)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_y_continuous(labels = abs) +  # Show absolute values on y-axis
  labs(
    x = "K-means Cluster",
    y = "Proportion",
    title = "Proportion of CAGs and Psychiatric Cases by Cluster"
  ) +
  scale_fill_manual(values = c("is_CAG" = "#D21E76", "is_psychiatric" = "#24A8BC")) +
  theme_minimal()







kmeans_classifier_psychiatric <- proportion_data %>%
  filter(classifier %in% c("is_CAG", "is_psychiatric")) %>%
  mutate(
    cluster_num = as.integer(str_extract(kmeans_cluster, "\\d+")),
    kmeans_cluster = factor(kmeans_cluster, levels = paste0("cluster", 1:50))
  ) %>%
  ggplot(aes(x = kmeans_cluster, y = proportion_yes, fill = classifier)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "K-means Cluster", y = "Proportion of 'yes'", fill = "Classifier") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  scale_fill_manual(
    values = c("is_CAG" = "#D21E76", "is_psychiatric" = "#24A8BC")
  ) +
  facet_wrap(~ classifier, ncol = 1, scales = "free_y")


# dev.size()
# width_original = 11.22
# height_original= 7.21
# 
# plot_name <- "chapter2_kmeans_showcase50.pdf"
# 
# pdf(plot_name, width = width_original, height =height_original)
# kmeans_classifier_psychiatric
# dev.off()




# ========== 4. EnrichR Results Processing ==========
#Add mouse human conversion
mouse_human_conversion <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/5_240219_dynamic_cluster_classifiers/human_mouse_conversion.txt") %>% 
  dplyr::select(human_gene_name, mouse_gene_name) 
mouse_human_conversion$human_gene_name[mouse_human_conversion$mouse_gene_name == "Hnrnpa1"] <- "HNRNPA1"

enrichr_results <- "/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/12_showcase_50K/enrichRresults"

# Get the list of files
files <- list.files(enrichr_results, pattern = "\\.txt$", full.names = TRUE)

# Use Map() to ensure file names are correctly assigned
data_list <- Map(function(file) {
  read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}, files)

# Assign names explicitly to ensure correct correspondence
names(data_list) <- basename(files)

# Print names to check
print(names(data_list))  # Should match the exact file names

# Standardize case & remove extra spaces in the conversion table
mouse_human_conversion <- mouse_human_conversion %>%
  mutate(human_gene_name = str_trim(toupper(human_gene_name)),  # Trim spaces & make uppercase
         mouse_gene_name = str_trim(mouse_gene_name))  # Trim spaces

# Function to convert human genes to mouse genes
convert_to_mouse_genes <- function(gene_string, conversion_table) {
  human_genes <- str_trim(toupper(unlist(strsplit(gene_string, ";"))))  # Standardize input genes
  matched_genes <- conversion_table %>%
    dplyr::filter(human_gene_name %in% human_genes) %>%
    pull(mouse_gene_name)  # Extract corresponding mouse genes
  
  return(ifelse(length(matched_genes) > 0, paste(matched_genes, collapse = ";"), ""))  # Join genes back
}

# Apply the function to the "Genes" column in all dataframes inside data_list
data_list <- lapply(data_list, function(df) {
  if ("Genes" %in% colnames(df)) {  # Ensure "Genes" column exists
    df <- df %>%
      mutate(mouse_gene = sapply(Genes, convert_to_mouse_genes, conversion_table = mouse_human_conversion))
  }
  return(df)
})

# Create a mapping from DEGs_kmeans_classifier
category_mapping <- DEGs_kmeans_classifier %>%
  dplyr::select(gene, is_CAG, is_addiction, is_psychiatric) %>%
  dplyr::rename(is_psycho = is_psychiatric)  # Rename column to match request

# Function to filter genes by category
filter_genes_by_category <- function(gene_string, category_column) {
  genes <- unlist(strsplit(gene_string, ";"))  # Split multiple genes
  matched_genes <- dplyr::filter(category_mapping, gene %in% genes, !!dplyr::sym(category_column) == "yes") %>%  # Match genes in the category
    dplyr::pull(gene)  # Extract the matched gene names
  
  return(ifelse(length(matched_genes) > 0, paste(matched_genes, collapse = ";"), ""))  # Join or return empty
}

# Apply the transformation to each dataframe in data_list
data_list <- lapply(data_list, function(df) {
  if ("mouse_gene" %in% colnames(df)) {  # Ensure "mouse_gene" column exists
    df <- df %>%
      dplyr::mutate(
        is_CAG = sapply(mouse_gene, filter_genes_by_category, category_column = "is_CAG"),
        is_addiction = sapply(mouse_gene, filter_genes_by_category, category_column = "is_addiction"),
        is_psycho = sapply(mouse_gene, filter_genes_by_category, category_column = "is_psycho")
      )
  }
  return(df)
})

# Check the first few rows of the first dataframe
head(data_list[[8]] %>% dplyr::select(Term, Adjusted.P.value, mouse_gene, is_CAG, is_addiction, is_psycho))

#saveRDS(data_list, "enrichr_results_kmeans_selected.rds")



# ========== 5. Line plots for specific examples:  ==========
gene_expression_df <- as.data.frame(read_tsv("obs_average_with_reps.tsv", col_names = TRUE))

# Extract timepoint from timepoint_replicate
gene_expression_df <- gene_expression_df %>%
  mutate(timepoint = gsub("_R\\d+", "", timepoint_replicate))  # Remove replicate suffix

# Define the correct order of timepoints
timepoint_order <- c("control", "h1", "h4", "h8", "h24", "d14")
gene_expression_df$timepoint <- factor(gene_expression_df$timepoint, levels = timepoint_order, ordered = TRUE)

# Convert to long format
df_long <- gene_expression_df %>%
  pivot_longer(cols = -c(timepoint, timepoint_replicate), 
               names_to = "gene", values_to = "expression")

selected_genes <- c("Stom", "Lamtor1", "Src", "Tomm22")
df_filtered <- df_long %>%
  dplyr::filter(gene %in% selected_genes)

df_summary <- df_filtered %>%
  dplyr::group_by(timepoint, gene) %>%
  dplyr::summarise(mean_expression = mean(expression, na.rm = TRUE),
                   sd_expression = sd(expression, na.rm = TRUE),
                   .groups = "drop")


cluster42_gene_examples <- ggplot(df_summary, aes(x = timepoint, y = mean_expression, group = gene, color = gene)) +
  geom_line(size = 0.5, linetype = "solid") +  # Thicker, solid lines for both genes
  geom_point(size = 2) +  # Mean values as points
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), 
                width = 0.1, alpha = 0.5) +  # Error bars with transparency
  scale_color_manual(values = c("#E63946", "#46579B", "#F4A261", "#2A9D8F")) +  # Custom colors
  scale_x_discrete(labels = c("sal", "1h", "4h", "8h", "24h", "14d")) +  # Custom x-axis labels
  theme_classic() +
  labs(title = "",
       x = "",
       y = "Mean Gene Expression ± SD") +
  theme(legend.title = element_blank())
cluster42_gene_examples


#Save
# dev.size()
# width_original = 3.808140
# height_original= 2.174419
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter2_kmeans_cluster42_examples"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# cluster42_gene_examples
# dev.off()



