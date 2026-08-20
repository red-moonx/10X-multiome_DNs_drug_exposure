# ===========================================
# Script Title: GO Terms (as in fig2)
# Author: Luna Zea Redondo
# Date: 2025-12-15
# Description:
#   This script processes and visualizes GO enrichment results from EnrichR.
#   It selects representative GO terms with minimal gene set overlap for 
#   each memory class and direction (up/down), and plots selected terms 
#   per class.
# ===========================================

# ========== Set Environment ==========
rm(list = ls(all.names = TRUE))
gc()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

# ========== Load Required Libraries ==========
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(readr)
library(glue)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(forcats)
library(data.table)
library(patchwork)

# ========== Set Working Directory ==========
setwd("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/11_Fig2_GOterms")

# and define directories
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/2024/5_figures/240910_GOterms_plots/enrichR_tables/"
dirplots <- "/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/11_Fig2_GOterms/"

# ========== 1. Load and Parse EnrichR Results ==========
# Read all GO term tables
# Extract metadata from filenames (memory class, direction, database)

# Get list of files in the directory
files <- list.files(path = dir, pattern = "*.txt", full.names = TRUE)

# Function to read files and add metadata columns
read_enrichr_file <- function(file) {
  data <- fread(file)
  data$memory_class <- str_extract(basename(file), "transient|recovered|memory|delayed")
  data$direction <- ifelse(grepl("UP", basename(file)), "UP", "DOWN")
  data$database <- str_extract(basename(file), "BioPlanet_2019|SynGO_2024")
  return(data)
}

# Read and process EnrichR result files
enrichr_data <- bind_rows(lapply(files, read_enrichr_file)) %>% 
  distinct(Term, memory_class, direction, database, .keep_all = TRUE)

# Function to calculate gene set overlap
gene_overlap <- function(genes1, genes2) {
  length(intersect(genes1, genes2))
}

# ========== 2. Select Top GO Terms by Memory Class ==========
# Define gene set overlap filter
# Select 3 non-overlapping terms per memory class

# Function to select terms with minimal overlap, ensuring the first term is the largest
define_least_overlap_terms <- function(data, n = 3) {
  selected_terms <- list()
  selected_genes <- list()
  
  # Sort terms by gene set size (largest first)
  sorted_data <- data %>% arrange(desc(lengths(Genes)))
  terms <- sorted_data$Term
  genes_list <- sorted_data$Genes
  
  # Select the first (largest) gene set
  selected_terms[[1]] <- terms[1]
  selected_genes[[1]] <- genes_list[[1]]
  
  for (i in 2:length(terms)) {
    term_i <- terms[i]
    genes_i <- genes_list[[i]]
    
    overlap <- sapply(selected_genes, function(genes) gene_overlap(genes, genes_i))
    
    if (all(overlap == 0)) {  # Only select if no overlap
      selected_terms[[length(selected_terms) + 1]] <- term_i
      selected_genes[[length(selected_genes) + 1]] <- genes_i
    }
    
    if (length(selected_terms) == n) {
      break
    }
  }
  
  return(selected_terms)
}

#p.threshold <- 0.05 # modify for more stringent results (e.g. 0.001) or use FDR 

# ========== 3. Plot: Upregulated GO Terms ==========
# Process only UP terms
top_terms_UP <- enrichr_data %>%
  mutate(memory_class = factor(memory_class, levels = c("transient", "recovered", "memory", "delayed"))) %>%
  arrange(memory_class) %>% 
  filter(direction == "UP", `P-value` < p.threshold) %>%
  group_by(memory_class, database) %>%
  mutate(Genes = str_split(Genes, ";")) %>%
  dplyr::summarise(selected_terms = list(define_least_overlap_terms(cur_data(), n = 3)), .groups = "drop") %>%
  unnest(selected_terms)

# Ensure selected terms are included for final dataset
top_terms_final <- enrichr_data %>%
  filter(Term %in% top_terms_UP$selected_terms) %>%
  group_by(memory_class, database, direction) %>%
  mutate(Overlap = parse_ratio(Overlap)) %>%
  arrange(factor(Term, levels = top_terms_UP$selected_terms))

# Reorder factors
top_terms_final$Term <- factor(top_terms_final$Term, levels = unique(top_terms_final$Term))
top_terms_final$memory_class <- factor(top_terms_final$memory_class, levels = c("transient", "recovered", "memory", "delayed"))

# Generate plots for UP terms only
for (db in unique(top_terms_final$database)) {
  enrichR_plot_up <- ggplot(top_terms_final %>% filter(database == db, direction == "UP"),
                            aes(x = memory_class, y = fct_rev(Term))) +
    geom_point(aes(size = Overlap, color = `P-value`)) +
    theme_classic() +
    theme_bw(base_size = 6) +
    ylab(NULL) +
    xlab(NULL) +
    facet_wrap(~ database, scales = "free", ncol = 1) +
    scale_colour_gradient2(limits = c(0, 0.1), low = "#cb181d", mid = "#fb6a4a", high = "#fcbba1",
                           na.value = "grey50", midpoint = 0.05) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title = element_blank())
  print(enrichR_plot_up)
}

# ========== 4. Plot: Downregulated GO Terms ==========
# Same as above, for direction == "DOWN"
# Process DOWN terms only
top_terms_DOWN <- enrichr_data %>%
  mutate(memory_class = factor(memory_class, levels = c("transient", "recovered", "memory", "delayed"))) %>%
  arrange(memory_class) %>% 
  filter(direction == "DOWN", `P-value` < p.threshold) %>%
  group_by(memory_class, database) %>%
  mutate(Genes = str_split(Genes, ";")) %>%
  dplyr::summarise(selected_terms = list(define_least_overlap_terms(cur_data(), n = 3)), .groups = "drop") %>%
  unnest(selected_terms)

# Ensure selected terms are included for final dataset
top_terms_final <- enrichr_data %>%
  filter(Term %in% top_terms_DOWN$selected_terms) %>%
  group_by(memory_class, database, direction) %>%
  mutate(Overlap = parse_ratio(Overlap)) %>%
  arrange(factor(Term, levels = top_terms_DOWN$selected_terms))

# Reorder factors
top_terms_final$Term <- factor(top_terms_final$Term, levels = unique(top_terms_final$Term))
top_terms_final$memory_class <- factor(top_terms_final$memory_class, levels = c("transient", "recovered", "memory", "delayed"))

# Generate plots for DOWN terms only
for (db in unique(top_terms_final$database)) {
  enrichR_plot_down <- ggplot(top_terms_final %>% filter(database == db, direction == "DOWN"),
                              aes(x = memory_class, y = fct_rev(Term))) +
    geom_point(aes(size = Overlap, color = `P-value`)) +
    theme_classic() +
    theme_bw(base_size = 6) +
    ylab(NULL) +
    xlab(NULL) +
    facet_wrap(~ database, scales = "free", ncol = 1) +
    scale_colour_gradient2(limits = c(0, 0.1), low = "#cb181d", mid = "#fb6a4a", high = "#fcbba1",
                           na.value = "grey50", midpoint = 0.05) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title = element_blank())
  print(enrichR_plot_down)
}

# ========== 5. Combine and Save Final Plot ==========
enrichR_plot_up + enrichR_plot_down

# dev.size()
# width_original = 15.395349
# height_original= 9.093023
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter2_DEGs_memroyClass_GOterm"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# enrichR_plot_up + enrichR_plot_down
# dev.off()

