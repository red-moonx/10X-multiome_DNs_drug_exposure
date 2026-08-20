# ===========================================
# Script Title: 4.1_mDN-PRC2_vs_DEGs.R
# Author: Luna Zea Redondo
# Date: 2026-05-20
# Description:
#   Compares PRC genes against cocaine DEGs and evaluates:
#     - overlap between PRC status and DEG status
#     - DEG direction vs PRC enrichment
#     - PRC gene proportions across memory classes and contrasts
#     - permutation-based significance for overlap counts
#     - enrichment plots for PRC gene sets
#     - example annotations for external gene categories
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
library(forcats)
library(Seurat)
library(biomaRt)
library(data.table)
library(enrichR)
library(patchwork)
library(gridExtra)

# ========== Project setup ==========
seu.VTA_DNs <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/2_DN-evaluation/250507_seu.VTA_DNs.rds")

dir <- "/fast/AG_Pombo/luna/2026_rebuttal/14_PRC2_vs_DEGs"
setwd(dir)

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  version = 112
)

# ===========================================
# 1) Format datasets to unify
# ===========================================

# Genes included in Toskas paper:
toskas22 <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/toskas_table_s1_chromatinStates.csv") %>% 
  dplyr::select(Gene, dat_4mo_wt_state)

toskas22_genes <- unique(toskas22$Gene)  
toskas22_mapping <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters    = "external_gene_name",
  values     = toskas22_genes,
  mart       = ensembl
) %>% 
  dplyr::rename(toskas_gene = "external_gene_name")

# Genes from multiome: muscat
DEG_complete_results <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/3_muscat/deg_results/260513_DEG_complete_results_Analysis1.tsv") %>% 
  dplyr::select(gene) %>% distinct() 

DEG_complete_genes <- unique(DEG_complete_results$gene)  
DEG_complete_mapping <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters    = "external_gene_name",
  values     = DEG_complete_genes,
  mart       = ensembl
) %>% 
  dplyr::rename(gene = "external_gene_name")

unified_DF <- left_join(DEG_complete_mapping, toskas22_mapping, by = "ensembl_gene_id") %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(gene) %>% distinct()

# FIX: Filter universe to genes detected in >= 10 cells (as we are using this filter throughout the entire manuscript)
# A gene that failed this threshold could never become a DEG in our analysis
# so must be excluded from the "no_deg" background
detection_counts <- rowSums(GetAssayData(seu.VTA_DNs, assay = "RNA", slot = "counts") > 0)
genes_detected_10cells <- names(detection_counts[detection_counts >= 10])

unified_DF <- unified_DF %>%
  dplyr::filter(gene %in% genes_detected_10cells)

message(glue("Genes in universe after detection filter: {nrow(unified_DF)}"))

# DEG datasets:
DEGs_memory_classes <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/6_TRACE/260516_DEG_complete_results_kmeans.tsv")
DEGs_contrast <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/3_muscat/deg_results/260513_DEG_complete_results_Analysis1.tsv")


# ===========================================
# 2) Q1: PRC genes vs cocaine-DEGs
# ===========================================

# 2A. Add PRC and DEG status, run chi-square test
PRC_degs_unified_df <- left_join(unified_DF, DEGs_memory_classes, by = "gene") %>% 
  dplyr::select(gene, memory_class) %>% 
  left_join(toskas22, by = c("gene" = "Gene")) %>% 
  dplyr::mutate(is_deg    = ifelse(is.na(memory_class), "no_deg", "deg")) %>% 
  dplyr::mutate(prc_state = ifelse(grepl("K27", dat_4mo_wt_state), "prc", "no_prc")) 

contingency_table <- table(PRC_degs_unified_df$is_deg, PRC_degs_unified_df$prc_state)

# Perform Chi-square test of independence
chi_square_test <- chisq.test(contingency_table)
print(chi_square_test) 


# 2B. DEG per direction
# Upregulated at any time, downregulated at any time, never deg
all_upregulated <- DEGs_contrast %>%
  dplyr::filter(significant == "Upregulated (pval < 0.05 and logFC > 0.5)") %>%
  dplyr::select(gene) %>% distinct() %>% pull()

all_downregulated <- DEGs_contrast %>%
  dplyr::filter(significant == "Downregulated (pval < 0.05 and logFC < -0.5)") %>%
  dplyr::select(gene) %>% distinct() %>% pull()

never_deg <- PRC_degs_unified_df %>%
  dplyr::filter(is_deg == "no_deg") %>%
  dplyr::select(gene) %>% distinct() %>% pull()

PRC_degs_unified_df.v2 <- PRC_degs_unified_df %>%
  dplyr::mutate(
    any_time_up   = ifelse(gene %in% all_upregulated,   "yes", "no"),
    any_time_down = ifelse(gene %in% all_downregulated, "yes", "no"),
    never_deg     = ifelse(gene %in% never_deg,         "yes", "no")
  ) %>%
  dplyr::select(prc_state, any_time_up, any_time_down, never_deg)

# Test any time up:
contingency_table_up <- table(PRC_degs_unified_df.v2$prc_state, PRC_degs_unified_df.v2$any_time_up)
chi_square_test_up   <- chisq.test(contingency_table_up)
print(chi_square_test_up)

# Test any time down:
contingency_table_down <- table(PRC_degs_unified_df.v2$prc_state, PRC_degs_unified_df.v2$any_time_down)
chi_square_test_down   <- chisq.test(contingency_table_down)
print(chi_square_test_down)

# Test never deg:
contingency_table_never <- table(PRC_degs_unified_df.v2$prc_state, PRC_degs_unified_df.v2$never_deg)
chi_square_test_never   <- chisq.test(contingency_table_never)
print(chi_square_test_never)


# 2C. DEG per memory class
# FIX: renamed to PRC_degs_unified_df_memclass to avoid overwriting the 2A object
PRC_degs_unified_df_memclass <- left_join(unified_DF, DEGs_memory_classes, by = "gene") %>% 
  dplyr::select(gene, memory_class, direction) %>% 
  left_join(toskas22, by = c("gene" = "Gene")) %>% 
  dplyr::mutate(category  = ifelse(is.na(memory_class), "no_deg", glue("{memory_class}_{direction}"))) %>% 
  dplyr::mutate(prc_state = ifelse(grepl("K27", dat_4mo_wt_state), "prc", "no_prc")) %>% 
  dplyr::select(category, prc_state)

# FIX: single omnibus chi-square across all categories (avoids correlated pairwise tests)
contingency_table_omnibus <- table(
  PRC_degs_unified_df_memclass$prc_state,
  PRC_degs_unified_df_memclass$category
)
chi_omnibus <- chisq.test(contingency_table_omnibus)
print(chi_omnibus)

# Post-hoc pairwise tests (only meaningful if omnibus is significant)
# FIX: Fisher fallback for rare categories with small expected counts
# FIX: BH correction applied across all 8 pairwise tests
categories <- c(
  "transient_up", "transient_down",
  "recovered_up", "recovered_down",
  "memory_up",    "memory_down",
  "delayed_up",   "delayed_down"
)

if (chi_omnibus$p.value < 0.05) {
  p_values <- list()
  for (cat in categories) {
    filtered_data     <- PRC_degs_unified_df_memclass %>%
      dplyr::filter(category %in% c("no_deg", cat))
    contingency_table <- table(filtered_data$prc_state, filtered_data$category == cat)
    
    test_result <- tryCatch(
      chisq.test(contingency_table),
      warning = function(w) fisher.test(contingency_table)
    )
    p_values[[cat]] <- test_result$p.value
  }
  
  p_values_df <- data.frame(
    category = names(p_values),
    p_value  = unlist(p_values)
  ) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"))
  
  print(p_values_df)
} else {
  message("Omnibus test not significant — skipping post-hoc pairwise tests.")
}


# 2D. Barplots for PRC genes per contrast + significance tests
all_contrasts <- c(
  "h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline",
  "h24_cocaine-saline", "d14_cocaine-saline"
)

DEGs_contrast.withPRC <- DEGs_contrast %>% 
  dplyr::filter(gene %in% unified_DF$gene) %>% 
  dplyr::select(gene, contrast, significant) %>% 
  left_join(toskas22, by = c("gene" = "Gene")) %>% 
  dplyr::mutate(prc_state = ifelse(grepl("K27", dat_4mo_wt_state), "prc", "no_prc")) 

DEGs_contrast_withPRC.plot <- DEGs_contrast.withPRC %>%
  dplyr::filter(contrast %in% all_contrasts[1:5]) %>%
  dplyr::select(gene, prc_state, contrast, significant) %>%
  dplyr::mutate(category = case_when(
    significant == "No significant"                                ~ "no_degs",
    significant == "Upregulated (pval < 0.05 and logFC > 0.5)"   ~ "upreg",
    significant == "Downregulated (pval < 0.05 and logFC < -0.5)"~ "downreg"
  )) %>%
  dplyr::select(-significant) %>%
  dplyr::group_by(contrast, category, prc_state) 

# Significance tests per contrast x category (upreg/downreg vs no_degs background)
# FIX: Fisher fallback + BH correction across all 10 tests (5 contrasts x 2 directions)
contrasts_test       <- levels(factor(DEGs_contrast_withPRC.plot$contrast))
categories_test      <- c("upreg", "downreg")
p_values_contrast    <- list()

for (cont in contrasts_test) {
  for (cat in categories_test) {
    filtered_data     <- DEGs_contrast_withPRC.plot %>%
      dplyr::filter(contrast == cont, category %in% c("no_degs", cat))
    contingency_table <- table(filtered_data$prc_state, filtered_data$category == cat)
    
    test_result <- tryCatch(
      chisq.test(contingency_table),
      warning = function(w) fisher.test(contingency_table)
    )
    key                      <- paste(cont, cat, sep = "_")
    p_values_contrast[[key]] <- list(
      contrast = cont,
      category = cat,
      p_value  = test_result$p.value
    )
  }
}

p_values_contrast_df <- bind_rows(p_values_contrast) %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"))

print(p_values_contrast_df)

# Calculate summary statistics
summary_data <- DEGs_contrast_withPRC.plot %>%
  group_by(contrast, category) %>%
  dplyr::summarise(
    total_genes = n(),
    prc_genes   = sum(prc_state == "prc", na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(
    percentage = (prc_genes / total_genes) * 100,
    contrast   = factor(contrast, levels = all_contrasts),
    category   = factor(category, levels = c("no_degs", "downreg", "upreg"))
  ) %>%
  left_join(p_values_contrast_df, by = c("contrast", "category")) %>%
  dplyr::mutate(sig_label = case_when(
    is.na(p_adj)  ~ "",
    p_adj < 0.001 ~ "***",
    p_adj < 0.01  ~ "**",
    p_adj < 0.05  ~ "*",
    TRUE          ~ "ns"
  ))

PRC2_contrast <- ggplot(summary_data, aes(x = contrast, y = percentage, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(
    aes(label = sprintf("%.1f%% (%d)", percentage, prc_genes)),
    position = position_dodge(width = 0.9),
    vjust = -0.5, size = 3
  ) +
  geom_text(
    aes(label = sig_label),
    position = position_dodge(width = 0.9),
    vjust = -2, size = 4, fontface = "bold"
  ) +
  scale_fill_manual(values = c("no_degs" = "gray90", "downreg" = "#89B8D2", "upreg" = "#EC5D5B")) +
  labs(
    title = "Percentage of PRC genes per contrast and category",
    x     = "Contrast",
    y     = "Percentage",
    fill  = "Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plotting per memory class
categories_order <- c(
  'no_deg', 'transient_up', 'transient_down', 'recovered_up', 'recovered_down', 
  'memory_up', 'memory_down', 'delayed_up', 'delayed_down'
)

summary_data_memclass <- PRC_degs_unified_df_memclass %>%
  group_by(category) %>%
  dplyr::summarise(
    total_count = n(),
    prc_count   = sum(prc_state == "prc")
  ) %>%
  mutate(
    percentage     = (prc_count / total_count) * 100,
    category       = factor(category, levels = categories_order)
  )

ggplot(summary_data_memclass, aes(x = category, y = percentage)) +
  geom_bar(stat = "identity", fill = "steelblue", position = "dodge") +
  geom_text(aes(label = paste0(round(percentage, 1), "% (", prc_count, ")")),
            vjust = -0.5, color = "black") +
  labs(
    title = "Percentage of PRC genes per memory class",
    x     = "Memory Class",
    y     = "Percentage (%)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine up and down categories into mirrored barplot
summary_data_memclass <- summary_data_memclass %>%
  mutate(
    category_group       = gsub("_up|_down", "", category),
    percentage_adjusted  = ifelse(grepl("down$", category), -percentage,
                                  ifelse(category == "no_deg", NA, percentage)),
    regulation           = ifelse(grepl("up$", category),   "Up",
                                  ifelse(grepl("down$", category), "Down", NA)),
    vjust                = ifelse(percentage_adjusted > 0, -0.5, 1.5)
  )

category_order <- c("transient", "recovered", "memory", "delayed")
summary_data_memclass$category_group <- factor(summary_data_memclass$category_group, levels = category_order)
summary_data_memclass$regulation     <- factor(summary_data_memclass$regulation,     levels = c("Up", "Down"))

no_deg_value <- summary_data_memclass$percentage[summary_data_memclass$category == "no_deg"]

p1 <- ggplot(summary_data_memclass %>% dplyr::filter(!is.na(percentage_adjusted)), 
             aes(x = category_group, y = percentage_adjusted, fill = regulation)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5) +
  geom_text(aes(label = paste0(round(abs(percentage), 1), "% (", prc_count, ")"), vjust = vjust), 
            color = "black") +
  geom_hline(yintercept =  no_deg_value, linetype = "dashed", color = "goldenrod", size = 1) +
  geom_hline(yintercept = -no_deg_value, linetype = "dashed", color = "goldenrod", size = 1) +
  labs(
    title = "Percentage of PRC Genes per Memory Class",
    x     = "Memory Class",
    y     = "Percentage (%)",
    fill  = "DEG"
  ) +
  scale_fill_manual(values = c("Up" = "#EC5D5B", "Down" = "#89B8D2")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1

# plot_name <- "chapter4_PRC_memoryDEGs_barplot.pdf"
# dev.size()
# pdf(plot_name, width = 5.45, height =5.75)
# p1
# dev.off()

# 2E. Permutation / empirical p-values for baseline vs memory classes
# Includes BH/FDR-adjusted p-values across all tested categories
# Flexible y-axis scale for each individual panel

PRC_degs_unified_df_memclass_withgene <- left_join(unified_DF, DEGs_memory_classes, by = "gene") %>% 
  dplyr::select(gene, memory_class, direction) %>% 
  left_join(toskas22, by = c("gene" = "Gene")) %>% 
  dplyr::mutate(category = ifelse(is.na(memory_class), "no_deg", glue("{memory_class}_{direction}"))) %>% 
  dplyr::mutate(prc_state = ifelse(grepl("K27", dat_4mo_wt_state), "prc", "no_prc")) %>% 
  dplyr::select(gene, category, prc_state)

summary_data_iteration <- summary_data_memclass %>%
  dplyr::mutate(category = factor(category, levels = categories_order)) %>% 
  dplyr::arrange(category) %>% 
  dplyr::filter(category != "no_deg") 

num_its <- 10000

perm_results <- list()

for (row in 1:nrow(summary_data_iteration)) {
  
  cat         <- summary_data_iteration[row, "category"]    %>% pull() %>% as.character()
  total_count <- summary_data_iteration[row, "total_count"] %>% pull()
  obs_count   <- summary_data_iteration[row, "prc_count"]   %>% pull()
  
  message(glue("Running: {cat}: {total_count} total DEGs and {obs_count} prc targets"))
  
  baseline_data <- PRC_degs_unified_df_memclass_withgene %>% 
    dplyr::filter(category == "no_deg")
  
  obs <- numeric(num_its)
  
  for (i in 1:num_its) {
    obs[i] <- baseline_data %>%
      dplyr::sample_n(size = total_count) %>%
      dplyr::filter(prc_state == "prc") %>%
      nrow()
  }
  
  obs_df <- data.frame(value = obs)
  
  res <- obs_df %>% 
    dplyr::mutate(
      r_down = value <= obs_count,
      r_up   = value >= obs_count
    ) %>% 
    dplyr::summarize(obs = min(sum(r_up), sum(r_down)))
  
  p_emp <- (res$obs + 1) / (num_its + 1)
  
  perm_results[[cat]] <- list(
    category    = cat,
    total_count = total_count,
    prc_count   = obs_count,
    p_emp       = p_emp,
    obs_df      = obs_df
  )
}

# Collect empirical p-values and calculate BH/FDR-adjusted p-values
pval_results <- data.frame(
  category    = sapply(perm_results, function(x) x$category),
  total_count = sapply(perm_results, function(x) x$total_count),
  prc_count   = sapply(perm_results, function(x) x$prc_count),
  p_emp       = sapply(perm_results, function(x) x$p_emp)
)

pval_results$p_adj <- p.adjust(pval_results$p_emp, method = "BH")

# Make permutation plots using both empirical p-value and FDR-adjusted p-value
# Each plot keeps its own flexible y-axis scale
plot_list <- list()

for (cat in pval_results$category) {
  
  obs_df    <- perm_results[[cat]]$obs_df
  obs_count <- perm_results[[cat]]$prc_count
  
  p_emp <- pval_results %>% 
    dplyr::filter(category == cat) %>% 
    dplyr::pull(p_emp)
  
  p_adj <- pval_results %>% 
    dplyr::filter(category == cat) %>% 
    dplyr::pull(p_adj)
  
  p1 <- ggplot(obs_df, aes(x = value)) +
    geom_density(fill = "gray90", alpha = 0.5) +
    geom_vline(xintercept = obs_count, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = glue("{cat}; p = {round(p_emp, 4)}; FDR = {round(p_adj, 4)}"),
      x = "",
      y = ""
    ) +
    theme_classic() +
    theme(plot.title = element_text(size = 10))
  
  plot_list <- append(plot_list, list(p1))
}

# Combine plots into panel
# Y-axis is flexible because no common ylim is imposed
panel <- plot_list[[1]] + plot_list[[3]] + plot_list[[5]] + plot_list[[7]] + 
  plot_list[[2]] + plot_list[[4]] + plot_list[[6]] + plot_list[[8]] +
  plot_layout(ncol = 4, nrow = 2)

panel

# Optional: view the p-value table
pval_results


# plot_name <- "chapter4_PRC_memoryDEGs_permutations.pdf"
# 
# dev.size()
# width_orig = 3.76
# height_orig = 2.39
# 
# pdf(plot_name, width = (width_orig*3), height =(height_orig*3))
# panel
# dev.off()

# ===========================================
# 3) EnrichR analyses and plots
# ===========================================

background_prc_genes <- unique(
  PRC_degs_unified_df_memclass_withgene$gene[
    PRC_degs_unified_df_memclass_withgene$prc_state == "prc" &
      !is.na(PRC_degs_unified_df_memclass_withgene$category)
  ]
)

writeLines(background_prc_genes, con = "background_all_tested_prc_genes.txt")

# 3A. Export PRC gene sets for each category
# FIX: use PRC_degs_unified_df_memclass_withgene (has gene column)
prc_genes  <- PRC_degs_unified_df_memclass_withgene[PRC_degs_unified_df_memclass_withgene$prc_state == "prc", ]
categories_split <- split(prc_genes, prc_genes$category)

for (cat in names(categories_split)) {
  genes    <- categories_split[[cat]]$gene
  filename <- paste0(cat, "_prc_genes.txt")
  writeLines(genes, con = filename)
}

degs      <- PRC_degs_unified_df_memclass_withgene[PRC_degs_unified_df_memclass_withgene$category != "no_deg", ]
deg_genes <- degs$gene
writeLines(deg_genes, con = "all_degs_genes.txt")

# 3B. Read and parse enrichR tables
dir   <- "/fast/AG_Pombo/luna/2026_rebuttal/14_PRC2_vs_DEGs/enrichR_results"
files <- list.files(path = dir, pattern = "*.txt", full.names = TRUE)

read_enrichr_file <- function(file) {
  data              <- fread(file)
  data$memory_class <- str_extract(basename(file), "transient|recovered|memory|delayed")
  data$direction    <- ifelse(grepl("UP", basename(file)), "UP", "DOWN")
  data$database     <- str_extract(basename(file), "GO_Biological|SynGO|Jensen|Reactome")
  return(data)
}

enrichr_data <- bind_rows(lapply(files, read_enrichr_file)) %>% 
  distinct(Term, memory_class, direction, database, .keep_all = TRUE)

parse_ratio <- function(x) {
  sapply(strsplit(x, "/"), function(parts) as.numeric(parts[1]) / as.numeric(parts[2]))
}

#### UP terms
top_terms_UP <- enrichr_data %>%
  dplyr::filter(direction == "UP") %>% 
  dplyr::group_by(memory_class, database) %>%
  dplyr::slice_min(`P-value`, n = 3, with_ties = FALSE) %>%
  dplyr::select(Term, `P-value`, memory_class) %>% 
  dplyr::ungroup() %>%
  mutate(memory_class = factor(memory_class, levels = c("transient", "recovered", "memory", "delayed"))) %>%
  arrange(database, memory_class, `P-value`) %>%  
  dplyr::distinct(Term, database)

top_terms_UP_data <- enrichr_data %>%
  dplyr::filter(Term %in% top_terms_UP$Term, direction == "UP") %>%
  dplyr::group_by(memory_class, database) %>%
  dplyr::mutate(Overlap = parse_ratio(Overlap)) %>%
  arrange(factor(Term, levels = top_terms_UP$Term))

top_terms_UP_data$Term         <- factor(top_terms_UP_data$Term,         levels = unique(top_terms_UP_data$Term))
top_terms_UP_data$memory_class <- factor(top_terms_UP_data$memory_class, levels = c("transient", "recovered", "memory", "delayed"))

for (db in unique(top_terms_UP_data$database)) {
  plot_file  <- file.path(glue("250108_enrichR_{db}_UP.pdf"))
  pdf(plot_file)
  enrichR_plot <- ggplot(top_terms_UP_data %>% dplyr::filter(database == db),
                         aes(x = memory_class, y = fct_rev(Term))) + 
    geom_point(aes(size = Overlap, color = `P-value`)) +
    theme_bw(base_size = 6) +
    ylab(NULL) + xlab(NULL) +
    facet_wrap(~ database, scales = "free", ncol = 1) +
    scale_colour_gradient2(limits = c(0, 0.1), low = "#cb181d", mid = "#fb6a4a", high = "#fcbba1",
                           na.value = "grey50", midpoint = 0.05) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title  = element_blank())
  print(enrichR_plot)
  dev.off()
}

#### DOWN terms
top_terms_DOWN <- enrichr_data %>%
  dplyr::filter(direction == "DOWN") %>% 
  dplyr::group_by(memory_class, database) %>%
  dplyr::slice_min(`P-value`, n = 3, with_ties = FALSE) %>%
  dplyr::select(Term, `P-value`, memory_class) %>% 
  dplyr::ungroup() %>%
  mutate(memory_class = factor(memory_class, levels = c("transient", "recovered", "memory", "delayed"))) %>%
  arrange(database, memory_class, `P-value`) %>%  
  dplyr::distinct(Term, database)

top_terms_DOWN_data <- enrichr_data %>%
  dplyr::filter(Term %in% top_terms_DOWN$Term, direction == "DOWN") %>%
  dplyr::group_by(memory_class, database) %>%
  dplyr::mutate(Overlap = parse_ratio(Overlap)) %>%
  arrange(factor(Term, levels = top_terms_DOWN$Term))

top_terms_DOWN_data$Term         <- factor(top_terms_DOWN_data$Term,         levels = unique(top_terms_DOWN_data$Term))
top_terms_DOWN_data$memory_class <- factor(top_terms_DOWN_data$memory_class, levels = c("transient", "recovered", "memory", "delayed"))

for (db in unique(top_terms_DOWN_data$database)) {
  plot_file  <- file.path(glue("240608_enrichR_{db}_DOWN.pdf"))
  pdf(plot_file)
  enrichR_plot <- ggplot(top_terms_DOWN_data %>% dplyr::filter(database == db, !is.na(memory_class)),
                         aes(x = memory_class, y = fct_rev(Term))) + 
    geom_point(aes(size = Overlap, color = `P-value`)) +
    theme_bw(base_size = 6) +
    ylab(NULL) + xlab(NULL) +
    facet_wrap(~ database, scales = "free", ncol = 1) +
    scale_colour_gradient2(limits = c(0, 0.1), low = "#cb181d", mid = "#fb6a4a", high = "#fcbba1",
                           na.value = "grey50", midpoint = 0.05) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title  = element_blank())
  print(enrichR_plot)
  dev.off()
}

databases_to_plot <- c("GO_Biological", "SynGO")

plots <- list()
for (db in databases_to_plot) {
  enrichR_plot <- ggplot(top_terms_UP_data %>% dplyr::filter(database == db),
                         aes(x = memory_class, y = fct_rev(Term))) + 
    geom_point(aes(size = Overlap, color = `P-value`)) +
    theme_bw(base_size = 6) +
    ylab(NULL) + xlab(NULL) +
    scale_colour_gradient2(limits = c(0, 0.1), low = "#cb181d", mid = "#fb6a4a", high = "#fcbba1",
                           na.value = "grey50", midpoint = 0.05) +
    theme(
      axis.text.x    = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y    = element_blank(),
      axis.title     = element_blank(),
      legend.position = "bottom",
      legend.box      = "horizontal"
    )
  plots[[db]] <- enrichR_plot
}

ggplot_legend <- function(a_gplot) {
  g      <- ggplotGrob(a_gplot)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

dummy_plot <- ggplot(top_terms_UP_data, aes(x = memory_class, y = fct_rev(Term))) +
  geom_point(aes(size = Overlap, color = `P-value`)) +
  scale_colour_gradient2(limits = c(0, 0.1), low = "#cb181d", mid = "#fb6a4a", high = "#fcbba1",
                         na.value = "grey50", midpoint = 0.05) +
  theme(legend.position = "none")

shared_legend <- ggplot_legend(dummy_plot)

combined_plot <- arrangeGrob(
  grobs         = c(lapply(plots, function(p) p)),
  ncol          = 2,
  layout_matrix = rbind(c(1, 2), c(3, 3)),
  heights       = c(1, 0.1)
)

# plot_name <- "chapter4_PRC_memoryDEGs_prc2_targets_GO.pdf"
# dev.size()
# pdf(plot_name, width = 8.15, height =6.63)
# grid::grid.draw(combined_plot)
# dev.off()

# ===========================================
# 4) Example annotation with external gene categories
# ===========================================

DEGs_for_examples      <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/4_DEGs_visualization/260513_DEGs_complete_info_simplified.tsv")
mouse_human_conversion <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/human_mouse_conversion.txt")

DEGs_for_examples <- DEGs_for_examples %>%
  left_join(mouse_human_conversion, by = c("gene" = "mouse_gene_name"))

DEGs_filtered <- DEGs_for_examples %>%
  dplyr::select(human_gene_name, is_CAG, is_addiction, is_psychiatric, is_neuropeptide, is_parkinson, is_als)

top_terms_with_categories <- top_terms_UP_data %>%
  rowwise() %>%
  mutate(
    is_CAG = paste(
      str_split(Genes, ";")[[1]] %>%
        .[. %in% DEGs_filtered$human_gene_name[DEGs_filtered$is_CAG == "yes"]],
      collapse = ";"
    ),
    is_anyAddiction = paste(
      str_split(Genes, ";")[[1]] %>%
        .[. %in% DEGs_filtered$human_gene_name[DEGs_filtered$is_addiction == "yes"]],
      collapse = ";"
    ),
    is_psychiatric = paste(
      str_split(Genes, ";")[[1]] %>%
        .[. %in% DEGs_filtered$human_gene_name[DEGs_filtered$is_psychiatric == "yes"]],
      collapse = ";"
    ),
    is_parkinson = paste(
      str_split(Genes, ";")[[1]] %>%
        .[. %in% DEGs_filtered$human_gene_name[DEGs_filtered$is_parkinson == "yes"]],
      collapse = ";"
    )
  ) %>%
  ungroup()

print(top_terms_with_categories)

#load("260607_mDN_PRC2-vs-DEGs.rds")
