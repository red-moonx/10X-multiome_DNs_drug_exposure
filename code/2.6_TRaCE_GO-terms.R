# ===========================================
# Script Title: 2.6. TRaCE GO analysis
# Author: Luna Zea Redondo
# Date: 2026-05-15
# Description:
#   This script processes and visualizes GO enrichment results from EnrichR.
#   It selects representative GO terms with minimal gene set overlap for 
#   each memory class and direction (up/down), and plots selected terms 
#   per class.
#   - Reactome: define_least_overlap_terms (non-overlapping gene sets)
#   - ChEA: top 3 by p-value (slice_min)
#   - Background: all other DEGs excluding the query class x direction
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
library(glue)
library(ggplot2)
library(readr)
library(stringr)
library(data.table)
library(forcats)
library(patchwork)

source("/fast/AG_Pombo/luna/2026_rebuttal/config.R", local = FALSE)

# ========== Set Working Directory ==========
dir <- "/fast/AG_Pombo/luna/2026_rebuttal/7_TRaCE_GO-terms"
setwd(dir)

parse_ratio <- function(x) {
  sapply(strsplit(x, "/"), function(parts) as.numeric(parts[1]) / as.numeric(parts[2]))
}

# ========== 0. Prepare gene lists and corresponding backgrounds ==========
# Background: all other DEGs excluding the query memory class x direction
# This answers: "what is special about this class vs all other cocaine-responsive genes?"

DEG_complete_results_kmeans <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/6_TRACE/260516_DEG_complete_results_kmeans.tsv")

# Create output directory if it doesn't exist
dir.create("enrichR_files", showWarnings = FALSE)

memory_classes <- c("transient", "recovered", "memory", "delayed")

for (mc in memory_classes) {
  for (dr in c("up", "down")) {
    
    # Target: genes in this specific memory class x direction
    target <- DEG_complete_results_kmeans %>%
      dplyr::filter(memory_class == mc, direction == dr) %>%
      dplyr::select(gene)
    
    # Background: all other DEGs excluding this specific memory class x direction
    background <- DEG_complete_results_kmeans %>%
      dplyr::filter(!(memory_class == mc & direction == dr)) %>%
      dplyr::select(gene)
    
    # Filenames
    target_file     <- file.path("enrichR_files", paste0(mc, toupper(dr), ".tsv"))
    background_file <- file.path("enrichR_files", paste0("bdg_", mc, toupper(dr), ".tsv"))
    
    # Write files
    write_tsv(target,     target_file,     col_names = FALSE)
    write_tsv(background, background_file, col_names = FALSE)
    
    message(glue("Written: {mc}_{dr} — {nrow(target)} query genes, {nrow(background)} background genes"))
  }
}

# From 07/07/26:
# ========== 1. Load and Parse EnrichR Results ==========

files <- list.files(path = glue("{dir}/enrichR_tables"), pattern = "*.txt", full.names = TRUE)

read_enrichr_file <- function(file) {
  data              <- fread(file)
  data$memory_class <- str_extract(basename(file), "transient|recovered|memory|delayed")
  data$direction    <- ifelse(grepl("UP", basename(file)), "UP", "DOWN")
  data$database     <- str_extract(basename(file), "Reactome|ChEA")
  return(data)
}

enrichr_data <- bind_rows(lapply(files, read_enrichr_file)) %>% 
  distinct(Term, memory_class, direction, database, .keep_all = TRUE)


# ========== 2. Helper functions ==========

gene_overlap <- function(genes1, genes2) {
  length(intersect(genes1, genes2))
}

# Select n terms with minimal gene set overlap (largest first)
define_least_overlap_terms <- function(data, n = 3) {
  if (nrow(data) == 0) return(character(0))
  
  data <- data %>%
    arrange(desc(lengths(Genes)))
  
  terms      <- data$Term
  genes_list <- data$Genes
  
  selected_terms <- character(0)
  selected_genes <- list()
  
  for (i in seq_along(terms)) {
    term_i  <- terms[i]
    genes_i <- genes_list[[i]]
    
    if (length(selected_genes) == 0) {
      selected_terms <- c(selected_terms, term_i)
      selected_genes[[1]] <- genes_i
    } else {
      overlap <- sapply(selected_genes, function(genes) {
        length(intersect(genes, genes_i))
      })
      
      if (all(overlap == 0)) {
        selected_terms <- c(selected_terms, term_i)
        selected_genes[[length(selected_genes) + 1]] <- genes_i
      }
    }
    
    if (length(selected_terms) == n) break
  }
  
  selected_terms
}

# Shared plot theme
enrichr_theme <- function() {
  theme_bw(base_size = 6) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title  = element_blank()
    )
}

# Shared color scale
enrichr_color_scale <- function() {
  scale_colour_gradient2(
    limits   = c(0, 0.1),
    low      = "#cb181d",
    mid      = "#fb6a4a",
    high     = "#fcbba1",
    na.value = "grey50",
    midpoint = 0.05
  )
}

# Shared color scale for DOWNregulated terms
enrichr_color_scale_down <- function() {
  scale_colour_gradient2(
    limits   = c(0, 0.1),
    low      = "#08519c",
    mid      = "#3182bd",
    high     = "#bdd7e7",
    na.value = "grey50",
    midpoint = 0.05
  )
}


# ========== 3. Reactome: non-overlapping terms ==========

# UP
reactome_top_UP <- enrichr_data %>%
  dplyr::filter(
    database == "Reactome",
    direction == "UP",
    `P-value` < 0.05
  ) %>%
  mutate(
    memory_class = factor(memory_class, levels = c("transient", "recovered", "memory", "delayed")),
    Genes        = str_split(Genes, ";")
  ) %>%
  arrange(memory_class) %>%
  group_by(memory_class) %>%
  dplyr::summarise(
    selected_terms = list(define_least_overlap_terms(cur_data(), n = 3)),
    .groups = "drop"
  ) %>%
  unnest(selected_terms) %>%
  rename(Term = selected_terms) %>%
  mutate(direction = "UP")


# DOWN
reactome_top_DOWN <- enrichr_data %>%
  dplyr::filter(
    database == "Reactome",
    direction == "DOWN",
    `P-value` < 0.05
  ) %>%
  mutate(
    memory_class = factor(memory_class, levels = c("transient", "recovered", "memory", "delayed")),
    Genes        = str_split(Genes, ";")
  ) %>%
  arrange(memory_class) %>%
  group_by(memory_class) %>%
  dplyr::summarise(
    selected_terms = list(define_least_overlap_terms(cur_data(), n = 3)),
    .groups = "drop"
  ) %>%
  unnest(selected_terms) %>%
  rename(Term = selected_terms) %>%
  mutate(direction = "DOWN") 


# Combine and prepare final Reactome dataset
# FIX: join on Term + direction + memory_class to avoid cross-contamination
reactome_top_all <- bind_rows(reactome_top_UP, reactome_top_DOWN)

reactome_final <- enrichr_data %>%
  dplyr::filter(database == "Reactome") %>%
  dplyr::inner_join(
    reactome_top_all %>%
      dplyr::select(Term, direction) %>%
      dplyr::distinct(),
    by = c("Term", "direction")
  ) %>%
  mutate(
    Overlap      = parse_ratio(Overlap),
    Term         = factor(Term, levels = unique(reactome_top_all$Term)),
    memory_class = factor(memory_class, levels = c("transient", "recovered", "memory", "delayed"))
  )

# Plot Reactome UP
reactome_plot_up <- ggplot(
  reactome_final %>% dplyr::filter(direction == "UP"),
  aes(x = memory_class, y = fct_rev(Term))
) +
  geom_point(aes(size = Overlap, color = `P-value`)) +
  ggtitle("Reactome — Upregulated") +
  ylab(NULL) + xlab(NULL) +
  enrichr_color_scale() +
  enrichr_theme()

# Plot Reactome DOWN
reactome_plot_down <- ggplot(
  reactome_final %>% dplyr::filter(direction == "DOWN"),
  aes(x = memory_class, y = fct_rev(Term))
) +
  geom_point(aes(size = Overlap, color = `P-value`)) +
  ggtitle("Reactome — Downregulated") +
  ylab(NULL) + xlab(NULL) +
  enrichr_color_scale_down() +
  enrichr_theme()

reactome_plot_up + reactome_plot_down


#save:
# dev.size()
# width_original  = 14.55
# height_original = 5.65
# plot_name  <- "chapter2_DEGs_memoryClass_GOterm.pdf"
# pdf(plot_name, width = width_original, height = height_original)
# reactome_plot_up + reactome_plot_down
# dev.off()



# ========== 4. ChEA: top 3 by p-value ==========

chea_data <- enrichr_data %>%
  dplyr::filter(database == "ChEA") %>%
  mutate(
    memory_class = factor(
      memory_class,
      levels = c("transient", "recovered", "memory", "delayed")
    )
  )

# Top 3 per memory class x direction by p-value
chea_top_terms <- chea_data %>%
  group_by(memory_class, direction) %>%
  slice_min(`P-value`, n = 3, with_ties = FALSE) %>%
  ungroup()

# Keep the selected terms, but plot their enrichment across all memory classes
# within the same direction.
#
# Important change:
# join by Term + direction, NOT by Term + direction + memory_class
chea_final <- chea_data %>%
  dplyr::inner_join(
    chea_top_terms %>%
      dplyr::select(Term, direction) %>%
      dplyr::distinct(),
    by = c("Term", "direction")
  ) %>%
  mutate(
    Overlap      = parse_ratio(Overlap),
    Term         = factor(Term, levels = unique(chea_top_terms$Term)),
    memory_class = factor(
      memory_class,
      levels = c("transient", "recovered", "memory", "delayed")
    )
  )

# Plot ChEA UP - reddish scale
chea_plot_up <- ggplot(
  chea_final %>% dplyr::filter(direction == "UP"),
  aes(x = memory_class, y = fct_rev(Term))
) +
  geom_point(aes(size = Overlap, color = `P-value`)) +
  ggtitle("ChEA — Upregulated") +
  ylab(NULL) + xlab(NULL) +
  enrichr_color_scale() +
  enrichr_theme()

# Plot ChEA DOWN - bluish scale
chea_plot_down <- ggplot(
  chea_final %>% dplyr::filter(direction == "DOWN"),
  aes(x = memory_class, y = fct_rev(Term))
) +
  geom_point(aes(size = Overlap, color = `P-value`)) +
  ggtitle("ChEA — Downregulated") +
  ylab(NULL) + xlab(NULL) +
  enrichr_color_scale_down() +
  enrichr_theme()

chea_plot_up + chea_plot_down


# ========== 5. Combine and Save ==========


# save:
# dev.size()
# width_original  = 10
# height_original = 5.65
# plot_name  <- "chapter2_DEGs_memoryClass_ChEA_results.pdf"
# pdf(plot_name, width = width_original, height = height_original)
# chea_plot_up + chea_plot_down
# dev.off()

