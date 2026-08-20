# ===========================================
# Script Title: 4.1_mDN-PRC2_vs_DEGs.R
# Author: Luna Zea Redondo
# Date: 2024-11-01
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

library(biomaRt)
library(data.table)

library(enrichR)
library(patchwork)
library(gridExtra)

# ========== Project setup ==========
source("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/7_240510_GEX_memory_only/240510_functions_upset.sh")

project_dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/8_241101_PRCgenes_vs_DEGs"
setwd(project_dir)

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")


# ===========================================
# 1) Format datasets to unify
# ===========================================

#Genes included in Toskas paper:
toskas22 <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/3_240605_PRCdata_Toskas/toskas_table_s1_chromatinStates.csv") %>% 
  dplyr::select(Gene, dat_4mo_wt_state)

toskas22_genes <- unique(toskas22$Gene)  
toskas22_mapping <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = toskas22_genes,
  mart = ensembl
) %>% 
  dplyr::rename(toskas_gene = "external_gene_name")

#Genes from multiome: muscat
DEG_complete_results <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002_muscat_VTAvsSN_all_subtypes/Analysis2/deg_results/231002_DEG_complete_results_Analysis2.tsv") %>% 
  dplyr::filter(cluster_id == "VTA") %>% 
  dplyr::select(gene) %>% distinct() 

DEG_complete_genes <- unique(DEG_complete_results$gene)  
DEG_complete_mapping <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = DEG_complete_genes,
  mart = ensembl
) %>% 
  dplyr::rename(gene = "external_gene_name")

unified_DF<- left_join(DEG_complete_mapping, toskas22_mapping, by = "ensembl_gene_id") %>% 
  filter(complete.cases(.)) %>% 
  dplyr::select(gene) %>% distinct()

#DEG datasets:
DEGs_memory_classes <- read_tsv("../data/DEGs_memory_classes.tsv")
DEGs_contrast <- read_tsv("../data/DEGs_per_contrast.tsv")

# ===========================================
# 2) Q1: PRC genes vs cocaine-DEGs
# ===========================================

# 2A. Add PRC and DEG status, run chi-square test
PRC_degs_unified_df <- left_join(unified_DF, DEGs_memory_classes, by = "gene") %>% 
  dplyr::select(gene, memory_class) %>% 
  left_join(toskas22, by =c("gene"="Gene")) %>% 
  dplyr::mutate(is_deg = ifelse(is.na(memory_class), "no_deg", "deg")) %>% 
  dplyr::mutate(prc_state = ifelse(grepl("K27", dat_4mo_wt_state), "prc", "no_prc")) 
contingency_table <- table(PRC_degs_unified_df$is_deg, PRC_degs_unified_df$prc_state)

#Perform Chi-square test of independence
chi_square_test <- chisq.test(contingency_table)
print(chi_square_test) 

# 2B. DEG per direction
# Upregulated at any time, downregulated at any time, never deg
all_upregulated <- DEGs_contrast %>% dplyr::filter(significant == "Upregulated (pval < 0.05 and logFC > 0.5)") %>%
  dplyr::select(gene) %>% distinct() %>% pull()
all_downregulated <- DEGs_contrast %>% dplyr::filter(significant == "Downregulated (pval < 0.05 and logFC < -0.5)") %>%
  dplyr::select(gene) %>% distinct() %>% pull()
never_deg <- PRC_degs_unified_df %>% dplyr::filter(is_deg == "no_deg") %>% dplyr::select(gene) %>% distinct() %>% pull() 

PRC_degs_unified_df.v2 <- PRC_degs_unified_df %>%
  dplyr::mutate(
    any_time_up = ifelse(gene %in% all_upregulated, "yes", "no"), 
    any_time_down = ifelse(gene %in% all_downregulated, "yes", "no"), 
    never_deg = ifelse(gene %in% never_deg, "yes", "no")) %>%
  dplyr::select(prc_state, any_time_up, any_time_down, never_deg)

#Test the any time up:
contingency_table_up <- table(PRC_degs_unified_df.v2$prc_state, PRC_degs_unified_df.v2$any_time_up)
chi_square_test_up <- chisq.test(contingency_table_up)
print(chi_square_test_up)

#Test the any time down:
contingency_table_down <- table(PRC_degs_unified_df.v2$prc_state, PRC_degs_unified_df.v2$any_time_down)
chi_square_test_down <- chisq.test(contingency_table_down)
print(chi_square_test_down)

#Test the never deg:
contingency_table_never <- table(PRC_degs_unified_df.v2$prc_state, PRC_degs_unified_df.v2$never_deg)
chi_square_test_never <- chisq.test(contingency_table_never)
print(chi_square_test_never)

# 2C. DEG per memory class
PRC_degs_unified_df <- left_join(unified_DF, DEGs_memory_classes, by = "gene") %>% 
  dplyr::select(gene, memory_class, direction) %>% 
  left_join(toskas22, by =c("gene"="Gene")) %>% 
  dplyr::mutate(category = ifelse(is.na(memory_class), "no_deg", glue("{memory_class}_{direction}"))) %>% 
  dplyr::mutate(prc_state = ifelse(grepl("K27", dat_4mo_wt_state), "prc", "no_prc")) %>% 
  dplyr::select(category, prc_state)

# Get the list of unique categories
categories <- unique(PRC_degs_unified_df$category)

p_values <- list()
# Loop over each category
for (cat in categories) {
  #Comparisons only between the memory cat and the "no_degs"
  filtered_data <- PRC_degs_unified_df %>% dplyr::filter(category %in% c("no_deg", cat))
  contingency_table <- table(PRC_degs_unified_df$prc_state, PRC_degs_unified_df$category == cat)
  chi_square_test <- chisq.test(contingency_table)
  p_values[[cat]] <- chi_square_test$p.value
}

p_values_df <- data.frame(
  category = names(p_values),
  p_value = unlist(p_values)
)

# 2D. Barplots for PRC genes per contrast and per memory class
all_contrasts <- c("h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline", "h24_cocaine-saline", "d14_cocaine-saline")

DEGs_contrast.withPRC <- DEGs_contrast %>% 
  dplyr::filter(gene %in% unified_DF$gene) %>% 
  dplyr::select(gene, contrast, significant) %>% 
  left_join(toskas22, by =c("gene"="Gene")) %>% 
  dplyr::mutate(prc_state = ifelse(grepl("K27", dat_4mo_wt_state), "prc", "no_prc")) 

DEGs_contrast_withPRC.plot <- DEGs_contrast.withPRC %>%
  dplyr::filter(contrast %in% all_contrasts[1:5]) %>%
  dplyr::select(gene, prc_state, contrast, significant) %>%
  dplyr::mutate(category = case_when(
    significant == "No significant" ~ "no_degs",
    significant == "Upregulated (pval < 0.05 and logFC > 0.5)" ~ "upreg",
    significant == "Downregulated (pval < 0.05 and logFC < -0.5)" ~ "downreg")) %>%
  dplyr::select(-significant) %>%
  dplyr::group_by(contrast, category, prc_state) 


# Calculate summary statistics
summary_data <- DEGs_contrast_withPRC.plot %>%
  group_by(contrast, category) %>%
  dplyr::summarise(
    total_genes = n(),  # Total genes in this contrast and category
    prc_genes = sum(prc_state == "prc", na.rm = TRUE),  # Count genes with 'prc' state
    .groups = "drop"
  ) %>%
  mutate(
    percentage = (prc_genes / total_genes) * 100,  # Calculate percentage of prc genes
    contrast = factor(contrast, levels = all_contrasts),
    category = factor(category, levels = c("no_degs", "downreg", "upreg"))  # Ensure to use correct levels if they are different
  )

ggplot(summary_data, aes(x = contrast, y = percentage, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.1f%% (%d)", percentage, prc_genes)), 
            position = position_dodge(width = 0.9), vjust = -0.5) +  # Add percentage and count labels
  scale_fill_manual(values = c("no_degs" = "gray90", "downreg" = "#89B8D2", "upreg" = "#EC5D5B")) +
  labs(title = "Percentage of PRC genes per contrast and category",
       x = "Contrast",
       y = "Percentage",
       fill = "Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Per memory class
summary_data <- PRC_degs_unified_df %>%
  group_by(category) %>%
  dplyr::summarise(
    total_count = n(),
    prc_count = sum(prc_state == "prc")
  ) %>%
  mutate(percentage = (prc_count / total_count) * 100) %>% 
  mutate(category = factor(category, levels = categories_order))

# Plotting the data
categories_order <- c('no_deg', 'transient_up', 'transient_down', 'recovered_up', 'recovered_down', 
                      'memory_up', 'memory_down', 'delayed_up', 'delayed_down')

ggplot(summary_data, aes(x = category, y = percentage)) +
  geom_bar(stat = "identity", fill = "steelblue", position = "dodge") +
  geom_text(aes(label = paste0(round(percentage, 1), "% (", prc_count, ")")),
            vjust = -0.5, color = "black") +
  labs(title = "Percentage of PRC genes per memory class",
       x = "Memory Class",
       y = "Percentage (%)",
       fill = "Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine up and down categories into the same x-axis point
summary_data <- summary_data %>%
  mutate(
    category_group = gsub("_up|_down", "", category), # Group by base category name
    percentage_adjusted = ifelse(grepl("down$", category), -percentage, 
                                 ifelse(category == "no_deg", NA, percentage)),
    regulation = ifelse(grepl("up$", category), "Up", 
                        ifelse(grepl("down$", category), "Down", NA)), # Adjust regulation labels
    vjust = ifelse(percentage_adjusted > 0, -0.5, 1.5) # Add vjust column
  )

# Set the custom order for the category_group
category_order <- c("transient", "recovered", "memory", "delayed")
summary_data$category_group <- factor(summary_data$category_group, levels = category_order)
summary_data$regulation <- factor(summary_data$regulation, levels = c("Up", "Down"))

# Extract the y-value for the dashed line (no_deg percentage)
no_deg_value <- summary_data$percentage[summary_data$category == "no_deg"]

# Plotting
p1 <- ggplot(summary_data %>% filter(!is.na(percentage_adjusted)), 
             aes(x = category_group, y = percentage_adjusted, fill = regulation)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5) +
  geom_text(aes(label = paste0(round(abs(percentage), 1), "% (", prc_count, ")"), vjust = vjust), 
            color = "black") +
  # Add the dashed line for no_deg
  geom_hline(yintercept = no_deg_value, linetype = "dashed", color = "goldenrod", size = 1) +
  geom_hline(yintercept = -no_deg_value, linetype = "dashed", color = "goldenrod", size = 1) +
  labs(title = "Percentage of PRC Genes per Memory Class",
       x = "Memory Class",
       y = "Percentage (%)",
       fill = "DEG") + # Update legend title
  scale_fill_manual(values = c("Up" = "#EC5D5B", "Down" = "#89B8D2")) + # Custom colors
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + # Optional: Percent format for y-axis
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter4_PRC_memoryDEGs_barplot"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# dev.size()
# pdf(file_name, width = 5.45, height =5.75)
# p1
# dev.off()

# 2E. Permutation / empirical p-values for baseline vs memory classes
# Initialize list to store p-values and contingency tables
p_values <- list()
contingency_tables <- list()

for (cat in categories[2:9]) {
  # Comparisons only between the memory cat and the "no_degs": the new addition
  filtered_data <- PRC_degs_unified_df %>% dplyr::filter(category %in% c("no_deg", cat))
  contingency_table <- table(filtered_data$prc_state, filtered_data$category == cat)
  
  # Perform chi-square test and store p-value
  chi_square_test <- chisq.test(contingency_table)
  p_values[[cat]] <- chi_square_test$p.value
  
  # Add contingency table to list with the category name
  contingency_tables[[cat]] <- as.data.frame.matrix(contingency_table)
}

# Create p-values DataFrame
p_values_df <- data.frame(
  category = names(p_values),
  p_value = unlist(p_values)
)

# Combine all contingency tables into a single DataFrame
contingency_tables_df <- do.call(rbind, contingency_tables)
rownames(contingency_tables_df) <- NULL  # Reset row names



# Using the baseline (no_degs)
##############################

PRC_degs_unified_df <- left_join(unified_DF, DEGs_memory_classes, by = "gene") %>% 
  dplyr::select(gene, memory_class, direction) %>% 
  left_join(toskas22, by =c("gene"="Gene")) %>% 
  dplyr::mutate(category = ifelse(is.na(memory_class), "no_deg", 
                                  glue("{memory_class}_{direction}"))) %>% 
  dplyr::mutate(prc_state = ifelse(grepl("K27", dat_4mo_wt_state), "prc", "no_prc"))  %>% 
  dplyr::select(gene, category, prc_state)

summary_data_iteration <- summary_data %>%
  dplyr::mutate(category = factor(category, levels = categories_order)) %>% 
  dplyr::arrange(category) %>% 
  dplyr::filter(category != "no_deg") 

plot_list <- list()
for (row in 1:nrow(summary_data_iteration)) {
  #Set parameters:
  cat = summary_data_iteration[row, "category"] %>% pull() %>% as.character()
  total_count = summary_data_iteration[row, "total_count"] %>% pull()
  obs_count = summary_data_iteration[row, "prc_count"] %>% pull()
  
  message(glue("Running: {cat}: {total_count} total DEGs and {obs_count} prc targets"))
  #Loop iterations:
  num_its=10000
  
  #Initate vector:
  obs <- numeric()
  #Initiate loop
  for (i in 1:num_its) {
    baseline_data <- PRC_degs_unified_df %>% dplyr::filter(category == "no_deg")
    num_prc_pos <- baseline_data %>%
      sample_n(size = total_count) %>%
      dplyr::filter(prc_state == "prc") %>%
      nrow()
    obs <- c(obs, num_prc_pos) 
  }
  #Convert obs vector in dataframe
  obs_df <- data.frame(value = obs)
  
  # Calculate the empirical p-value
  res <- obs_df %>% 
    mutate(
      r_down = value <= obs_count,   # Identify values less than or equal to obs_count
      r_up = value >= obs_count      # Identify values greater than or equal to obs_count
    ) %>% 
    dplyr::summarize(
      obs = min(sum(r_up), sum(r_down)) # Take the minimum of the two sums
    )
  
  # Compute the empirical p-value
  p_emp <- (res + 1) / (num_its + 1)
  
  p1 <- ggplot(obs_df, aes(x = value)) +
    geom_density(fill = "gray90", alpha = 0.5) + # Density plot
    geom_vline(xintercept = obs_count, color = "red", linetype = "dashed", size = 1) + # Vertical line at obs_count
    labs(
      title = glue("{cat}; p-val: {round(p_emp, 4)}"),
      x = "",
      y = ""
    ) +
    theme_classic() 
  
  
  
  
  plot_list <- append(plot_list, list(p1))
}

# Define the maximum y value across all plots
max_y <- max(unlist(lapply(plot_list, function(p) ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)))

# Apply consistent tick distribution to all plots
plot_list <- lapply(plot_list, function(p) {
  p + scale_y_continuous(
    breaks = seq(0, max_y, length.out = 6) # 6 evenly spaced ticks
  )
})
# Definir un número fijo de ticks y su distribución uniforme
num_ticks <- 3
tick_positions <- seq(0, 1, length.out = num_ticks) # Normalizar de 0 a 1

# Aplicar los mismos ticks a todos los gráficos
plot_list <- lapply(plot_list, function(p) {
  # Escalar los ticks según el rango del gráfico
  y_range <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
  scaled_ticks <- scales::rescale(tick_positions, to = y_range) # Ajustar los ticks al rango del gráfico
  
  p + scale_y_continuous(
    breaks = scaled_ticks, # Ticks uniformes ajustados al rango
    labels = scales::number_format(accuracy = 0.01) # Máximo 2 decimales
  )
})

plot_list <- lapply(plot_list, function(p) {
  p + theme(
    plot.title = element_text(size = 10) # Adjust the size to your preference
  )
})

panel <- plot_list[[1]] + plot_list[[3]] + plot_list[[5]] + plot_list[[7]] + 
  plot_list[[2]] + plot_list[[4]] + plot_list[[6]] + plot_list[[8]] +
  plot_layout(ncol = 4, nrow = 2)

#load("241211_PRC_empyrical_pval.rds")

# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter4_PRC_memoryDEGs_permutations"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# dev.size()
# width_orig = 3.76
# height_orig = 2.39
# 
# pdf(file_name, width = (width_orig*3), height =(height_orig*3))
# panel
# dev.off()


# ===========================================
# 3) EnrichR analyses and plots
# ===========================================

# 3A. Export PRC gene sets for each category
head(PRC_degs_unified_df)

# Filter rows where prc_state is "prc"
prc_genes <- PRC_degs_unified_df[PRC_degs_unified_df$prc_state == "prc", ]

# Split the data by category
categories <- split(prc_genes, prc_genes$category)

# Write each category to a file
for (cat in names(categories)) {
  genes <- categories[[cat]]$gene
  filename <- paste0(cat, "_prc_genes.txt")
  writeLines(genes, con = filename)
}
#BDG 1:
degs <- PRC_degs_unified_df[PRC_degs_unified_df$category != "no_deg", ]
deg_genes <- degs$gene
writeLines(deg_genes, con = "all_degs_genes.txt")

# 3B. Read and parse enrichR tables
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/8_241101_PRCgenes_vs_DEGs/250108_results_for_enrichR/enrichR_results"
# Get list of files in the directory
files <- list.files(path = dir, pattern = "*.txt", full.names = TRUE) #FROM HERE NEED TO BE DEBUGED. 

# Function to read files and add a column for the memory class, regulation direction, and database
read_enrichr_file <- function(file) {
  data <- fread(file)
  data$memory_class <- str_extract(basename(file), "transient|recovered|memory|delayed")
  data$direction <- ifelse(grepl("UP", basename(file)), "UP", "DOWN")
  data$database <- str_extract(basename(file), "GO_Biological|SynGO|Jensen|Reactome")
  return(data)
  Sys.sleep(1)
}

# 1. Read in the EnrichR result files
enrichr_data <- bind_rows(lapply(files, read_enrichr_file)) %>% 
  # mutate(Term = str_remove_all(Term, "\\d+"),
  #        Term = str_remove_all(Term, "\\s*\\(.*?\\)"))  %>% 
  distinct(Term, memory_class, direction, database, .keep_all = TRUE)  # Keep distinct rows based on specific columns



#### 2. Only UP terms: Prepare the data for plotting
# First, identify the top 3 terms for each groupID
top_terms_UP <- enrichr_data %>%
  dplyr::filter(direction == "UP") %>% 
  dplyr::group_by(memory_class, database) %>%
  dplyr::slice_min(`P-value`, n = 3, with_ties = FALSE) %>%
  dplyr::select(Term, `P-value`, memory_class) %>% 
  dplyr::ungroup()
# Reorder the table
top_terms_UP <- top_terms_UP %>%
  mutate(memory_class = factor(memory_class, levels = c("transient", "recovered", "memory", "delayed"))) %>%
  arrange(database, memory_class, `P-value`) %>%  
  dplyr::distinct(Term, database)  # Keep unique terms across all groups

# Now, ensure these terms are included for all groupID
top_terms_UP_data <- enrichr_data %>%
  dplyr::filter(Term %in% top_terms_UP$Term, 
                direction == "UP") %>%
  dplyr::group_by(memory_class, database) %>%
  dplyr::mutate(Overlap = parse_ratio(Overlap))%>%
  arrange(factor(Term, levels = top_terms_UP$Term))

# Reorder the factors
top_terms_UP_data$Term <- factor(top_terms_UP_data$Term, levels = unique(top_terms_UP_data$Term))
top_terms_UP_data$memory_class <- factor(top_terms_UP_data$memory_class, levels = c("transient", "recovered", "memory", "delayed") )
#top_terms_UP_data <- top_terms_UP_data %>% arrange(database, memory_class)

# Generate the plots
for (db in unique(top_terms_UP_data$database)) {
  plot_file <- file.path(glue("250108_enrichR_{db}_UP.pdf"))
  pdf(plot_file)
  enrichR_plot <- ggplot(top_terms_UP_data %>% dplyr::filter(database == db),
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
  print(enrichR_plot)
  dev.off()
}


#### 2. Only DOWN terms: Prepare the data for plotting
# First, identify the top 3 terms for each groupID
top_terms_DOWN <- enrichr_data %>%
  dplyr::filter(direction == "DOWN") %>% 
  dplyr::group_by(memory_class, database) %>%
  dplyr::slice_min(`P-value`, n = 3, with_ties = FALSE) %>%
  dplyr::select(Term, `P-value`, memory_class) %>% 
  dplyr::ungroup()
# Reorder the table
top_terms_DOWN <- top_terms_DOWN %>%
  mutate(memory_class = factor(memory_class, levels = c("transient", "recovered", "memory", "delayed"))) %>%
  arrange(database, memory_class, `P-value`) %>%  
  dplyr::distinct(Term, database)  # Keep unique terms across all groups

# Now, ensure these terms are included for all groupID
top_terms_DOWN_data <- enrichr_data %>%
  dplyr::filter(Term %in% top_terms_DOWN$Term, 
                direction == "DOWN") %>%
  dplyr::group_by(memory_class, database) %>%
  dplyr::mutate(Overlap = parse_ratio(Overlap))%>%
  arrange(factor(Term, levels = top_terms_DOWN$Term))

# Reorder the factors
top_terms_DOWN_data$Term <- factor(top_terms_DOWN_data$Term, levels = unique(top_terms_DOWN_data$Term))
top_terms_DOWN_data$memory_class <- factor(top_terms_DOWN_data$memory_class, levels = c("transient", "recovered", "memory", "delayed") )
#top_terms_DOWN_data <- top_terms_DOWN_data %>% arrange(database, memory_class)


# Generate the plots
for (db in unique(top_terms_DOWN_data$database)) {
  plot_file <- file.path(dirplots, glue("240608_enrichR_{db}_DOWN.pdf"))
  pdf(plot_file)
  enrichR_plot <- ggplot(top_terms_DOWN_data %>% dplyr::filter(database == db, !is.na(memory_class)),
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
  print(enrichR_plot)
  dev.off()
}



#250109 - Plots only for figure:
# Generate individual plots for each database and store them in a list
plots <- list()
for (db in databases_to_plot) {
  enrichR_plot <- ggplot(top_terms_UP_data %>% filter(database == db),
                         aes(x = memory_class, y = fct_rev(Term))) + 
    geom_point(aes(size = Overlap, color = `P-value`)) +
    theme_classic() +
    theme_bw(base_size = 6) +
    ylab(NULL) +
    xlab(NULL) +
    scale_colour_gradient2(limits = c(0, 0.1), low = "#cb181d", mid = "#fb6a4a", high = "#fcbba1",
                           na.value = "grey50", midpoint = 0.05) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_blank(),  # Remove Y axis text
      axis.title = element_blank(),
      legend.position = "bottom",   # Position the "Overlap" legend below the plot
      legend.box = "horizontal"     # Make the legend horizontal
    )
  plots[[db]] <- enrichR_plot
}

# Extract a shared legend for the P-value scale (horizontal and below)
ggplot_legend <- function(a_gplot) {
  g <- ggplotGrob(a_gplot)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

# Generate a dummy plot to extract the shared P-value legend
dummy_plot <- ggplot(top_terms_UP_data,
                     aes(x = memory_class, y = fct_rev(Term))) +
  geom_point(aes(size = Overlap, color = `P-value`)) +
  scale_colour_gradient2(limits = c(0, 0.1), low = "#cb181d", mid = "#fb6a4a", high = "#fcbba1",
                         na.value = "grey50", midpoint = 0.05) +
  theme(legend.position = "bottom")

shared_legend <- ggplot_legend(dummy_plot)

# Combine all plots side by side and place the shared legend below
combined_plot <- arrangeGrob(
  grobs = c(lapply(plots, function(p) p), list(shared_legend)),
  ncol = 2,  # Arrange plots side by side
  layout_matrix = rbind(c(1, 2), c(3, 3)),  # Shared legend spans both plots
  heights = c(1, 0.1)  # Adjust the height of the shared legend
)
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter4_PRC_memoryDEGs_GOterms"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# dev.size()
# width_orig = 5.81
# height_orig = 5.44
# 
# pdf(file_name, width = (width_orig*1), height =(height_orig*1))
# grid.draw(combined_plot)
# dev.off()


# ===========================================
# 4) Example annotation with external gene categories
# ===========================================

DEGs_for_examples <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/9_241104_CAGs_psychiatric_association/241104_DEG_complete_results_allGeneSets_allGenes.tsv")
mouse_human_conversion <- readRDS("mouse_human_conversion.rds")

DEGs_for_examples <- DEGs_for_examples %>% left_join(mouse_human_conversion, by = c("gene" = "mouse_gene_name"))
# Prepare the DEGs table by selecting relevant columns
DEGs_filtered <- DEGs_for_examples %>%
  dplyr::select(human_gene_name, is_CAG, is_anyAddiction, is_psychiatric, is_alzheimer, is_parkinson)

# Add the new columns to `top_terms_UP_data`
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
        .[. %in% DEGs_filtered$human_gene_name[DEGs_filtered$is_anyAddiction == "yes"]],
      collapse = ";"
    ),
    is_psychiatric = paste(
      str_split(Genes, ";")[[1]] %>% 
        .[. %in% DEGs_filtered$human_gene_name[DEGs_filtered$is_psychiatric == "yes"]],
      collapse = ";"
    ),
    is_alzheimer = paste(
      str_split(Genes, ";")[[1]] %>% 
        .[. %in% DEGs_filtered$human_gene_name[DEGs_filtered$is_alzheimer == "yes"]],
      collapse = ";"
    ),
    is_parkinson = paste(
      str_split(Genes, ";")[[1]] %>% 
        .[. %in% DEGs_filtered$human_gene_name[DEGs_filtered$is_parkinson == "yes"]],
      collapse = ";"
    )
  ) %>%
  ungroup()

# View the updated data
print(top_terms_with_categories)

