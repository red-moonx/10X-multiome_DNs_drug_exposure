# ===========================================
# Script Title: 4.4_PRC_work_exploration_ATAC.R
# Author: Luna Zea Redondo
# Date: 2024-11-13
# Description:
#   Analyses PRC overlap and depletion patterns in midbrain ATAC peaks.
#   Main steps:
#     - summarise PRC peak length and diff/LFC distributions
#     - compare PRC peak presence across contrasts and memory classes
#     - run empirical permutation tests against the no_dar baseline
#     - generate MEME-CHIP input files for PRC-DARs and PRC peak backgrounds
#     - examine peak-type composition across PRC classes
#     - test PRC depletion effects on ATAC signal
#     - integrate PRC-overlapping peaks with DEG annotations
#     - assess TF motif enrichment for PRC vs non-PRC peaks
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
library(scales)

library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(IRanges)

library(patchwork)

# ========== Project setup ==========
source("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/7_240510_GEX_memory_only/240510_functions_upset.sh")
setwd("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/9_PRCwork_with_ATAC")


# ========== Analysis parameters ==========
all_contrasts <- c(
  "h1_cocaine-saline",
  "h4_cocaine-saline",
  "h8_cocaine-saline",
  "h24_cocaine-saline",
  "d14_cocaine-saline",
  "h4_cocaine-h1_cocaine",
  "h8_cocaine-h4_cocaine",
  "h24_cocaine-h8_cocaine",
  "d14_cocaine-h24_cocaine"
)

categories_order <- c(
  "no_dar",
  "transient_up", "transient_down",
  "recovered_up", "recovered_down",
  "memory_up", "memory_down",
  "delayed_up", "delayed_down"
)

memory_classes_colors <- c(
  "no_dar" = "gray90",
  "transient_up" = "#D9481C",
  "recovered_up" = "#F3A712",
  "memory_up" = "#008AB8",
  "delayed_up" = "#954B77",
  "transient_down" = "#EF9A81",
  "recovered_down" = "#F8CD77",
  "memory_down" = "#8FD0E6",
  "delayed_down" = "#C892B3"
)

# ========== 1. Load input tables ==========
#Read files:
all_peaks_withPRC <- read_tsv("241215_all_peaks_withPRC.tsv")
all_peaks_withPRC_withDARs_contrast <- read_tsv("241215_all_peaks_withPRC_withDARs_contrast.tsv")
all_peaks_withPRC_withDARs_memory <- read_tsv("241215_all_peaks_withPRC_withDARs_memory.tsv")


# ========== 2. PRC length distribution and PRC peak overview ==========
all_peaks_withPRC %>%
  dplyr::filter(!is.na(PRC_peakID)) %>%
  mutate(PRC_length = as.numeric(str_extract(PRC_peakID, "(?<=-)[0-9]+")) -
           as.numeric(str_extract(PRC_peakID, "(?<=:)[0-9]+"))) %>%
  ggplot(aes(x = PRC_length)) +
  geom_density(fill = "black", alpha = 0.5) +  # Density curve
  scale_x_log10() +  # Apply log10 transformation to x-axisDEG_complete_results
  labs(title = "Density of PRC Peak Lengths (Log Scale)",
       x = "Length of PRC Peaks (log10)",
       y = "Density") +
  theme_minimal()

#Plot PRC peaks:

all_peaks_withPRC.plot.df <- all_peaks_withPRC %>% 
  dplyr::filter(!is.na(PRC_peakID)) %>% 
  dplyr::distinct(PRC_peakID, .keep_all = TRUE) %>% 
  dplyr::select(WT4mon_vs_Mut4mon_diff, WT4mon_vs_Mut4mon_LFC)


all_peaks_withPRC.plot<- ggplot(all_peaks_withPRC.plot.df, aes(x = WT4mon_vs_Mut4mon_LFC, y = WT4mon_vs_Mut4mon_diff)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(
    x = "Log Fold Change (LFC)",
    y = "Difference",
    title = "Scatter Plot of Diff vs LFC"
  ) +
  theme_classic() +
  xlim(-3, 6) 

#median:
# Load necessary libraries
median_length <- all_peaks_withPRC %>%
  dplyr::filter(!is.na(PRC_peakID)) %>%
  mutate(PRC_length = as.numeric(str_extract(PRC_peakID, "(?<=-)[0-9]+")) -
           as.numeric(str_extract(PRC_peakID, "(?<=:)[0-9]+"))) %>%
  summarise(median_length = median(PRC_length, na.rm = TRUE))


# ========== 3. Distribution of PRC peaks per contrast ==========
all_peaks_withPRC_withDARs_contrast.plot <- all_peaks_withPRC_withDARs_contrast %>%
  dplyr::filter(contrast %in% all_contrasts[1:5]) %>%
  dplyr::mutate(direction = ifelse(logFC > 0, "upreg", "downreg")) %>%
  dplyr::select(peakID, PRC_peakID, contrast, diff, direction) %>%
  dplyr::mutate(category = case_when(
    diff == "No" ~ "no_dars",
    diff == "Yes" & direction == "upreg" ~ "upreg",
    diff == "Yes" & direction == "downreg" ~ "downreg")) %>%
  dplyr::group_by(contrast, category, PRC_peakID_present = !is.na(PRC_peakID)) 

summary_data <- all_peaks_withPRC_withDARs_contrast.plot %>%
  group_by(contrast, category) %>%
  dplyr::summarise(
    total_peakIDs = n(),
    true_peakIDs = sum(PRC_peakID_present == TRUE)
  ) %>%
  mutate(percentage = (true_peakIDs / total_peakIDs) * 100) %>%
  ungroup() %>% 
  mutate(
    contrast = factor(contrast, levels = all_contrasts[1:5]),  # Reorder contrasts
    category = factor(category, levels = c("no_dars", "downreg", "upreg"))  # Reorder categories
  )
ggplot(summary_data, aes(x = contrast, y = percentage, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("no_dars" = "gray90", 
                               "upreg" = "#EC5D5B", 
                               "downreg" = "#89B8D2")) +
  geom_text(aes(label = sprintf("%.1f%% (%d)", percentage, true_peakIDs)), 
            position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Percentage of PRC_peakID present per contrast and category",
       x = "Contrast",
       y = "Percentage",
       fill = "Category") +
  theme_minimal()


# ========== 4. Distribution of DARs and non-DARs per PRC category ==========
dars_vsnondars <- all_peaks_withPRC_withDARs_memory %>% 
  dplyr::select(peakID, PRC_peakID, memory_class, direction, peakType) %>% 
  dplyr::mutate(category = ifelse(is.na(memory_class), "no_dar", direction), 
                PRC_peakID = ifelse(is.na(PRC_peakID), "no_prc", "prc")) %>% 
  dplyr::select(category, PRC_peakID, peakType) %>% 
  dplyr::group_by(category, PRC_peakID) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100) %>%
  dplyr::ungroup() 

dars_nodars_order <- c("no_dar", "down", "up")
prc_order <- c("no_prc", "prc")

dars_vsnondars.plot <- dars_vsnondars %>% 
  dplyr::filter(PRC_peakID == "prc")

PRC_direction <- ggplot(dars_vsnondars.plot, aes(x = category, y = percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = sprintf("%.1f%% (%d)", percentage, count)), 
            vjust = -0.3, position = position_stack(vjust = 1)) + # Label with percentage and count
  scale_fill_manual(values = c("no_dar" = "gray90", 
                               "up" = "#EC5D5B", 
                               "down" = "#89B8D2")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

dars_vsnondars.plot$category = factor(dars_vsnondars.plot$category, levels = c("no_dar", "down", "up"))
PRC_direction <- ggplot(dars_vsnondars.plot, aes(x = category, y = percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = sprintf("%.1f%% (%d)", percentage, count)), 
            vjust = -0.3, position = position_stack(vjust = 1)) + # Label with percentage and count
  #  scale_x_discrete(limits = prc_order) +
  #  scale_fill_manual(values = c("prc" = "darkred", "no_prc" = "gray80")) +  # Assign colors
  scale_fill_manual(values = c("prc" = "black", "no_prc" = "gray80")) +  # Assign colors
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ factor(category, levels = c("no_dar", "down", "up"))) 

dars_up <- dars_vsnondars %>% dplyr::filter(direction == "up") %>% pull(peakID)
dars_down <- dars_vsnondars %>% dplyr::filter(direction == "down") %>% pull(peakID)

PRC_peaks <- all_peaks_withPRC.plot + PRC_direction

# dev.size()
# width_original = 10.19792
# height_original= 4.40625
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter4_PRC_peaks_combined"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf")
# 
# pdf(file_name, width = width_original, height =height_original)
# PRC_peaks
# dev.off()

# ========== 5. Distribution of PRC-marked peaks per memory class ==========

PRCpositive_memory_class.df <- all_peaks_withPRC_withDARs_memory %>% 
  dplyr::select(peakID, PRC_peakID, memory_class, direction) %>% 
  dplyr::mutate(category = ifelse(is.na(memory_class), "no_dar", glue("{memory_class}_{direction}")), 
                PRC_peakID = ifelse(is.na(PRC_peakID), "no_prc", "prc")) %>% 
  dplyr::select(category, PRC_peakID) %>% 
  dplyr::group_by(category, PRC_peakID) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(category = factor(category, levels = categories_order)) %>% 
  dplyr::arrange(category) %>% 
  dplyr::filter(PRC_peakID == "prc")


ggplot(PRCpositive_memory_class.df, aes(x = category, y = percentage, fill = category)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = paste0(round(percentage, 1), "% (", count, ")")), vjust = -0.5, size = 4) +
  scale_fill_manual(values = memory_classes_colors) +
  theme_minimal() +
  labs(title = "PRC Positive Memory Class Distribution", x = "Category", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# Plotting
summary_data <- PRCpositive_memory_class.df %>%
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
no_dar_value <- summary_data$percentage[summary_data$category == "no_dar"]

# Plotting
p1 <- ggplot(summary_data %>% dplyr::filter(!is.na(category_group)), 
             aes(x = category_group, y = percentage_adjusted, fill = regulation)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5) +
  geom_text(aes(label = paste0(round(abs(percentage), 1), "% (", count, ")"), vjust = vjust), 
            color = "black") +
  # Add the dashed line for no_deg
  geom_hline(yintercept = no_dar_value, linetype = "dashed", color = "goldenrod", size = 1) +
  geom_hline(yintercept = -no_dar_value, linetype = "dashed", color = "goldenrod", size = 1) +
  labs(title = "Percentage of PRC-DARs per Memory Class",
       x = "Memory Class",
       y = "Percentage (%)",
       fill = "DAR") + # Update legend title
  scale_fill_manual(values = c("Up" = "#EC5D5B", "Down" = "#89B8D2")) + # Custom colors
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + # Optional: Percent format for y-axis
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter4_PRC_memoryDARs_barplot"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# dev.size()
# pdf(file_name, width = 5.45, height =5.75)
# p1
# dev.off()
# 

# ========== 6. Empirical permutation test against no_dar baseline ==========
PRCpositive_memory_class.p_emp <- all_peaks_withPRC_withDARs_memory  %>% 
  dplyr::mutate(category = ifelse(is.na(memory_class), "no_dar", glue("{memory_class}_{direction}")), 
                PRC_peakID = ifelse(is.na(PRC_peakID), "no_prc", "prc")) %>% 
  dplyr::select(peakID, category, PRC_peakID)

summary_data <- PRCpositive_memory_class.p_emp %>%
  dplyr::group_by(category) %>%
  dplyr::summarise(total_count = n(), prc_count = sum(PRC_peakID == "prc")) %>%
  dplyr::mutate(percentage = (prc_count / total_count) * 100) %>% 
  dplyr::mutate(category = factor(category, levels = categories_order))

summary_data_iteration <- summary_data %>%
  dplyr::mutate(category = factor(category, levels = categories_order)) %>% 
  dplyr::arrange(category) %>% 
  dplyr::filter(!is.na(category)) 

plot_list <- list()
for (row in 1:nrow(summary_data_iteration)) {
  #Set parameters:
  cat = summary_data_iteration[row, "category"] %>% pull() %>% as.character()
  total_count = summary_data_iteration[row, "total_count"] %>% pull()
  obs_count = summary_data_iteration[row, "prc_count"] %>% pull()
  
  message(glue("Running: {cat}: {total_count} total DARs and {obs_count} prc targets"))
  #Loop iterations:
  num_its=10000
  
  #Initate vector:
  obs <- numeric()
  #Initiate loop
  for (i in 1:num_its) {
    baseline_data <- PRCpositive_memory_class.p_emp %>% dplyr::filter(category == "no_dar")
    num_prc_pos <- baseline_data %>%
      sample_n(size = total_count) %>%
      dplyr::filter(PRC_peakID == "prc") %>%
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
    geom_vline(xintercept = obs_count, color = "red", linetype = "dashed", size = 1) + # Vertical line at 99
    labs(
      title = glue("{cat}; p-val: {round(p_emp, 4)}; ATAC "),
      x = "",
      y = ""
    ) +
    theme_classic()
  
  
  plot_list <- append(plot_list, list(p1))
}

# panel <- plot_list[[1]] + plot_list[[3]] + plot_list[[5]] + plot_list[[7]] + 
#   plot_list[[2]] + plot_list[[4]] + plot_list[[6]] + plot_list[[8]] +
#   plot_layout(ncol = 4, nrow = 2)

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

pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
plot_name <- "chapter4_PRC_memoryDARs_permutations"
file_name <- glue("{pdf_dir}/{plot_name}.pdf" )

dev.size()
width_orig = 3.76
height_orig = 2.39

pdf(file_name, width = (width_orig*3), height =(height_orig*3))
panel
dev.off()

#save.image("250115_PRC_empyrical_pval_ATAC.rds")

# ========== 7. Chi-square tests for PRC enrichment across memory classes ==========
PRC_dars_unified_df <- all_peaks_withPRC_withDARs_memory %>% 
  dplyr::mutate(category = ifelse(is.na(memory_class), "no_dar", glue("{memory_class}_{direction}")), 
                prc_state = ifelse(is.na(PRC_peakID), "no_prc", "prc")) %>% 
  dplyr::select(peakID, category, prc_state) %>% 
  dplyr::distinct(peakID, .keep_all = TRUE)

for (cat in categories_order[2:9]) {
  # Comparisons only between the memory cat and the "no_degs": the new addition
  filtered_data <- PRC_dars_unified_df %>% dplyr::filter(category %in% c("no_dar", cat))
  contingency_table <- table(filtered_data$prc_state, filtered_data$category == cat)
  
  # Perform chi-square test and store p-value
  chi_square_test <- chisq.test(contingency_table)
  raw_p_value <- chi_square_test$p.value
  adjusted_p_value <- p.adjust(raw_p_value, method = "bonferroni")
  
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

# ========== 8. MEME-CHIP / FASTA preparation ==========
output_dir <- "250115_MEME_CHIP"

# Set genome
genome <- BSgenome.Mmusculus.UCSC.mm10

# 1. Extract all DARs (category != "no_dar")
all_dars <- PRCpositive_memory_class.p_emp %>%
  filter(category != "no_dar") %>%
  mutate(
    chr = sub(":.*", "", peakID),
    start = as.numeric(sub("-.*", "", sub(".*:", "", peakID))),
    end = as.numeric(sub(".*-", "", peakID))
  ) %>%
  dplyr::select(chr, start, end, peakID, category)

# Write BED file for all DARs
all_dars_bed_path <- file.path(output_dir, "all_dars.bed")
write.table(all_dars, all_dars_bed_path, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Get FASTA sequences for all DARs with IDs
all_dars_gr <- GRanges(seqnames = all_dars$chr, 
                       ranges = IRanges(start = all_dars$start, end = all_dars$end))
all_dars_fasta <- getSeq(genome, all_dars_gr)
names(all_dars_fasta) <- all_dars$peakID  # Use peakID as the sequence ID
all_dars_fasta_path <- file.path(output_dir, "all_dars.fasta")
writeXStringSet(all_dars_fasta, all_dars_fasta_path)

# 2. Extract PRC-DARs directly from the original dataframe
prc_dars <- PRCpositive_memory_class.p_emp %>%
  filter(category != "no_dar" & PRC_peakID != "no_prc") %>%
  mutate(
    chr = sub(":.*", "", peakID),
    start = as.numeric(sub("-.*", "", sub(".*:", "", peakID))),
    end = as.numeric(sub(".*-", "", peakID))
  ) %>%
  dplyr::select(chr, start, end, peakID, category)

# Write BED file for PRC-DARs
prc_dars_bed_path <- file.path(output_dir, "prc_dars.bed")
write.table(prc_dars, prc_dars_bed_path, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Get FASTA sequences for PRC-DARs with IDs
prc_dars_gr <- GRanges(seqnames = prc_dars$chr, 
                       ranges = IRanges(start = prc_dars$start, end = prc_dars$end))
prc_dars_fasta <- getSeq(genome, prc_dars_gr)
names(prc_dars_fasta) <- prc_dars$peakID  # Use peakID as the sequence ID
prc_dars_fasta_path <- file.path(output_dir, "prc_dars.fasta")
writeXStringSet(prc_dars_fasta, prc_dars_fasta_path)

# 3. Process PRC peaks by category and generate backgrounds
categories <- unique(prc_dars$category)

for (cat in categories) {
  # Subset PRC-DARs by category
  prc_dars_subset <- prc_dars %>%
    filter(category == cat)
  
  # Write BED file for this category
  bed_file <- file.path(output_dir, paste0("prc_dars_", cat, ".bed"))
  write.table(prc_dars_subset, bed_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  # Get FASTA sequences for this category with IDs
  prc_dars_cat_gr <- GRanges(seqnames = prc_dars_subset$chr, 
                             ranges = IRanges(start = prc_dars_subset$start, end = prc_dars_subset$end))
  prc_dars_cat_fasta <- getSeq(genome, prc_dars_cat_gr)
  names(prc_dars_cat_fasta) <- prc_dars_subset$peakID  # Use peakID as the sequence ID
  
  # Write FASTA file for this category
  fasta_file <- file.path(output_dir, paste0("prc_dars_", cat, ".fasta"))
  writeXStringSet(prc_dars_cat_fasta, fasta_file)
  
  # Create background files for this category
  bdg_all_dars <- all_dars %>%
    filter(category != cat)
  bdg_prc_dars <- prc_dars %>%
    filter(category != cat)
  
  # Write BED files for backgrounds
  bdg_all_dars_file <- file.path(output_dir, paste0("bdg_all_dars_for_", cat, ".bed"))
  write.table(bdg_all_dars, bdg_all_dars_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  bdg_prc_dars_file <- file.path(output_dir, paste0("bdg_prc_dars_for_", cat, ".bed"))
  write.table(bdg_prc_dars, bdg_prc_dars_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  # Get FASTA sequences for background DARs
  bdg_all_dars_gr <- GRanges(seqnames = bdg_all_dars$chr, 
                             ranges = IRanges(start = bdg_all_dars$start, end = bdg_all_dars$end))
  bdg_all_dars_fasta <- getSeq(genome, bdg_all_dars_gr)
  names(bdg_all_dars_fasta) <- bdg_all_dars$peakID  # Use peakID as the sequence ID
  bdg_all_dars_fasta_file <- file.path(output_dir, paste0("bdg_all_dars_for_", cat, ".fasta"))
  writeXStringSet(bdg_all_dars_fasta, bdg_all_dars_fasta_file)
  
  bdg_prc_dars_gr <- GRanges(seqnames = bdg_prc_dars$chr, 
                             ranges = IRanges(start = bdg_prc_dars$start, end = bdg_prc_dars$end))
  bdg_prc_dars_fasta <- getSeq(genome, bdg_prc_dars_gr)
  names(bdg_prc_dars_fasta) <- bdg_prc_dars$peakID  # Use peakID as the sequence ID
  bdg_prc_dars_fasta_file <- file.path(output_dir, paste0("bdg_prc_dars_for_", cat, ".fasta"))
  writeXStringSet(bdg_prc_dars_fasta, bdg_prc_dars_fasta_file)
}

# 250119 - Alternative background (all prc peaks)
##################################
output_dir <- "250115_MEME_CHIP"

# Set genome
genome <- BSgenome.Mmusculus.UCSC.mm10

# 1. Extract all DARs (category != "no_dar")
all_prc_peaks <- PRCpositive_memory_class.p_emp %>%
  dplyr::filter(PRC_peakID != "no_prc") %>%
  mutate(
    chr = sub(":.*", "", peakID),
    start = as.numeric(sub("-.*", "", sub(".*:", "", peakID))),
    end = as.numeric(sub(".*-", "", peakID))
  ) %>%
  dplyr::select(chr, start, end, peakID, category)

# Write BED file for all prc peaks
all_prc_peaks_bed_path <- file.path(output_dir, "all_prc_peaks.bed")
write.table(all_prc_peaks, all_prc_peaks_bed_path, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Get FASTA sequences for all DARs with IDs
all_prc_peaks_gr <- GRanges(seqnames = all_prc_peaks$chr, 
                            ranges = IRanges(start = all_prc_peaks$start, end = all_prc_peaks$end))
all_prc_peaks_fasta <- getSeq(genome, all_prc_peaks_gr)
names(all_prc_peaks_fasta) <- all_prc_peaks$peakID  # Use peakID as the sequence ID
all_prc_peaks_fasta_path <- file.path(output_dir, "all_prc_peaks.fasta")
writeXStringSet(all_prc_peaks_fasta, all_prc_peaks_fasta_path)

# 2. Extract PRC-DARs directly from the original dataframe
prc_dars <- PRCpositive_memory_class.p_emp %>%
  dplyr::filter(category != "no_dar" & PRC_peakID != "no_prc") %>%
  mutate(
    chr = sub(":.*", "", peakID),
    start = as.numeric(sub("-.*", "", sub(".*:", "", peakID))),
    end = as.numeric(sub(".*-", "", peakID))
  ) %>%
  dplyr::select(chr, start, end, peakID, category)

# Write BED file for PRC-DARs
prc_dars_bed_path <- file.path(output_dir, "prc_dars.bed")
#write.table(prc_dars, prc_dars_bed_path, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Get FASTA sequences for PRC-DARs with IDs
prc_dars_gr <- GRanges(seqnames = prc_dars$chr, 
                       ranges = IRanges(start = prc_dars$start, end = prc_dars$end))
prc_dars_fasta <- getSeq(genome, prc_dars_gr)
names(prc_dars_fasta) <- prc_dars$peakID  # Use peakID as the sequence ID
prc_dars_fasta_path <- file.path(output_dir, "prc_dars.fasta")
#writeXStringSet(prc_dars_fasta, prc_dars_fasta_path)

# 3. Process PRC peaks by category and generate backgrounds
categories <- unique(prc_dars$category)

for (cat in categories) {
  # Subset PRC-DARs by category
  prc_dars_subset <- prc_dars %>%
    dplyr::filter(category == cat)
  
  # Write BED file for this category
  bed_file <- file.path(output_dir, paste0("prc_dars_", cat, ".bed"))
  write.table(prc_dars_subset, bed_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  # Get FASTA sequences for this category with IDs
  prc_dars_cat_gr <- GRanges(seqnames = prc_dars_subset$chr, 
                             ranges = IRanges(start = prc_dars_subset$start, end = prc_dars_subset$end))
  prc_dars_cat_fasta <- getSeq(genome, prc_dars_cat_gr)
  names(prc_dars_cat_fasta) <- prc_dars_subset$peakID  # Use peakID as the sequence ID
  
  # Write FASTA file for this category
  fasta_file <- file.path(output_dir, paste0("prc_dars_", cat, ".fasta"))
  writeXStringSet(prc_dars_cat_fasta, fasta_file)
  
  # Create background files for this category
  bdg_all_prc_peaks <- all_prc_peaks %>%
    dplyr::filter(category != cat)
  
  # Write BED files for backgrounds
  bdg_all_prc_peaks_file <- file.path(output_dir, paste0("bdg_all_prc_peaks_for_", cat, ".bed"))
  write.table(bdg_all_prc_peaks, bdg_all_prc_peaks_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  # Get FASTA sequences for background DARs
  bdg_all_prc_peaks_gr <- GRanges(seqnames = bdg_all_prc_peaks$chr, 
                                  ranges = IRanges(start = bdg_all_prc_peaks$start, end = bdg_all_prc_peaks$end))
  bdg_all_prc_peaks_fasta <- getSeq(genome, bdg_all_prc_peaks_gr)
  names(bdg_all_prc_peaks_fasta) <- bdg_all_prc_peaks$peakID  # Use peakID as the sequence ID
  bdg_all_prc_peaks_fasta_file <- file.path(output_dir, paste0("bdg_all_prc_peaks_for_", cat, ".fasta"))
  writeXStringSet(bdg_all_prc_peaks_fasta, bdg_all_prc_peaks_fasta_file)
}

# IMPORTANT: Final background: prc peaks that are no_dars
all_prc_peaks_no_dars <- all_prc_peaks %>% dplyr::filter(category == "no_dar")
all_prc_peaks_no_dars_bed_path <- file.path(output_dir, "all_prc_peaks_no_dars.bed")
write.table(all_prc_peaks_no_dars, all_prc_peaks_no_dars_bed_path, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Get FASTA sequences for all DARs with IDs
all_prc_peaks_no_dars_gr <- GRanges(seqnames = all_prc_peaks_no_dars$chr, 
                                    ranges = IRanges(start = all_prc_peaks_no_dars$start, end = all_prc_peaks_no_dars$end))
all_prc_peaks_no_dars_fasta <- getSeq(genome, all_prc_peaks_no_dars_gr)
names(all_prc_peaks_no_dars_fasta) <- all_prc_peaks_no_dars$peakID  # Use peakID as the sequence ID
all_prc_peaks_no_dars_fasta_path <- file.path(output_dir, "all_prc_peaks_no_dars.fasta")
writeXStringSet(all_prc_peaks_no_dars_fasta, all_prc_peaks_no_dars_fasta_path)


# ========== 9. Peak type distribution by PRC status and accessibility ==========
dars_vsnondars_peaktype <- all_peaks_withPRC_withDARs_memory %>% 
  dplyr::select(peakID, PRC_peakID, memory_class, direction, peakType) %>% 
  dplyr::mutate(category = ifelse(is.na(memory_class), "no_dar", direction), 
                PRC_peakID = ifelse(is.na(PRC_peakID), "no_prc", "prc"), 
                category = factor(category, levels = c("no_dar", "down", "up")), 
                PRC_peakID = factor(PRC_peakID, levels = c("no_prc", "prc"))) %>% 
  dplyr::select(category, PRC_peakID, peakType) %>% 
  dplyr::group_by(category, PRC_peakID, peakType) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(percentage = count / sum(count) * 100) %>%
  dplyr::ungroup() 

peak_type_colors <- c(Promoter = "#8B0000", Exonic = "#DAA520", Intronic = "#2E8B57", Distal = "#4682B4")

#Calculate percentage out of all DARs
PRC_peaktype_thesis <- ggplot(dars_vsnondars_peaktype, aes(x = PRC_peakID, y = percentage, fill = peakType)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_fill_manual(values = peak_type_colors) +
  facet_wrap(~ category) +
  theme_minimal() +
  labs(
    title = "Peak Type Distribution by PRC Status and Accessibility",
    x = "PRC Status", y = "Percentage", fill = "Peak Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

# dev.size()
# width_original = 10.19792
# height_original= 4.760417
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter4_PRC_peaktype_thesis"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf")
# 
# pdf(file_name, width = width_original, height =height_original)
# PRC_peaktype_thesis
# dev.off()


#For paper, with coord_flip
#Calculate percentage out of all DARs
PRC_peaktype_paper <- ggplot(dars_vsnondars_peaktype, aes(x = PRC_peakID, y = percentage, fill = peakType)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_fill_manual(values = peak_type_colors) +
  facet_wrap(~ category) +
  theme_minimal() +
  labs(
    title = "Peak Type Distribution by PRC Status and Accessibility",
    x = "PRC Status", y = "Percentage", fill = "Peak Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  ) + coord_flip()

# dev.size()
# width_original = 7.562500
# height_original= 4.760417
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter4_PRC_peaktype_paper"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf")
# 
# pdf(file_name, width = width_original, height =height_original)
# PRC_peaktype_paper
# dev.off()

# ========== 10. PRC depletion effects ==========

# Distribution PRC marked peak type per memory class: 
memory_percentage_data_peaktype <- all_peaks_withPRC_withDARs_memory %>% 
  dplyr::select(peakID, PRC_peakID, memory_class, direction, peakType) %>% 
  dplyr::mutate(category = ifelse(is.na(memory_class), "no_dar", glue("{memory_class}_{direction}")), 
                PRC_peakID = ifelse(is.na(PRC_peakID), "no_prc", "prc")) %>% 
  dplyr::select(PRC_peakID, category, peakType) %>% 
  dplyr::group_by(category, PRC_peakID, peakType) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::group_by(category, PRC_peakID) %>% 
  dplyr::mutate(percentage = count / sum(count) * 100) %>%
  dplyr::ungroup()

categories_order <- c('no_dar', 'transient_up', 'transient_down', 'recovered_up', 'recovered_down', 
                      'memory_up', 'memory_down', 'delayed_up', 'delayed_down')

# Plot the data
# plot 1
ggplot(memory_percentage_data_peaktype, aes(x = PRC_peakID, y = percentage, fill = peakType)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            position = position_stack(vjust = 0.5)) +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Percentage of PRC and no_PRC Peaks in Each Category",
       x = "Category",
       y = "Percentage",
       fill = "PRC Peak ID") +
  facet_grid(~ factor(category, levels = categories_order))

#plot 2:
memory_percentage_data_peaktype_v2 <- memory_percentage_data_peaktype %>% filter(category != "no_dar")
ggplot(memory_percentage_data_peaktype_v2, aes(x = PRC_peakID, y = count, fill = peakType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  # Use position_dodge to avoid stacking
  geom_text(aes(label = count), 
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +  # Place the labels above the bars
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Percentage of PRC and no_PRC Peaks in Each Category",
       x = "Category",
       y = "Count",
       fill = "Peak Type") +
  facet_grid(~ factor(category, levels = categories_order))


#Histogram PRCpeaks vs DARs
PRC_depletion_effect.histogram <- all_peaks_withPRC_withDARs_memory %>% 
  dplyr::select(peakID, PRC_peakID, nearestGene, peakType, WT4mon_vs_Mut8mon_diff, WT4mon_vs_Mut8mon_LFC, memory_class, direction) %>% 
  dplyr::filter(!is.na(PRC_peakID)) %>% 
  dplyr::filter(!is.na(memory_class)) %>% 
  dplyr::select(PRC_peakID) %>% 
  dplyr::group_by(PRC_peakID) %>%
  dplyr::summarise(duplicate_count = n())

ggplot(PRC_depletion_effect.histogram, aes(x = duplicate_count)) +
  geom_histogram(binwidth = 1, fill = "steelblue") +
  stat_bin(geom = "text", binwidth = 1, aes(label = ..count..), 
           vjust = -0.5, size = 5) +  
  theme_minimal() +
  labs(title = "Distribution of DARs overlapping the same PRC peak",
       x = "Number of DARs",
       y = "Frequency")

#All memory categories together, diff vs logfc 8 months
PRC_depletion_effect <- all_peaks_withPRC_withDARs_memory %>% 
  dplyr::mutate(category = ifelse(is.na(memory_class), "no_dar", glue("{memory_class}_{direction}"))) %>% 
  dplyr::select(peakID, PRC_peakID, nearestGene, peakType, WT4mon_vs_Mut8mon_diff, WT4mon_vs_Mut8mon_LFC, category) %>% 
  dplyr::filter(!is.na(PRC_peakID)) 


# Create the scatter plot
PRC_depletion_effect$category <- factor(PRC_depletion_effect$category, 
                                        levels = c("no_dar", "transient_up", "recovered_up", 
                                                   "memory_up", "delayed_up", "transient_down", 
                                                   "recovered_down", "memory_down", "delayed_down"))
#Plot1: all together
ggplot(PRC_depletion_effect, aes(x = WT4mon_vs_Mut8mon_LFC , y = WT4mon_vs_Mut8mon_diff, color = category)) +
  geom_point(data = subset(PRC_depletion_effect, category == "no_dar"),
             aes(color = category), size = 1, color = "gray90", alpha=0.5) +
  geom_point(data = subset(PRC_depletion_effect, category != "no_dar"),
             aes(color = category), size = 1.2, alpha=0.75) +
  scale_color_manual(values = memory_classes_colors) +
  theme_minimal() +
  labs(title = "Scatter plot of WT4mon_vs_Mut8mon_diff vs LFC",
       x = "WT4mon_vs_Mut8mon_LFC",
       y = "WT4mon_vs_Mut8mon_diff",
       color = "Category") +
  theme(legend.position = "bottom")

#Plot2: facet
# Create the faceted scatter plot
PRC_depletion_effect <- PRC_depletion_effect %>%
  group_by(category) %>%
  mutate(category_label = paste0(category, " (n=", n(), ")"))  %>% 
  ungroup() %>%
  mutate(category_label = factor(category_label, levels = unique(paste0(levels(category), 
                                                                        " (n=", table(category)[levels(category)], ")"))))

ggplot(PRC_depletion_effect, aes(x = WT4mon_vs_Mut8mon_LFC, y = WT4mon_vs_Mut8mon_diff)) +
  geom_point(data = PRC_depletion_effect %>% filter(category == "no_dar"),
             aes(color = "no_dar"), size = 2) +
  geom_point(data = PRC_depletion_effect %>% filter(category != "no_dar"),
             aes(color = category), size = 2) +
  scale_color_manual(values = memory_classes_colors) +  
  #  scale_y_log10() +
  theme_minimal() +
  labs(title = "PRC peaks overlapping ATAC peaks: coloured by DAR memory class",
       x = "WT / Mut8months (LFC)",
       y = "WT - Mut8months",
       color = "Category") +
  theme(legend.position = "bottom") +
  facet_wrap(~ category_label)


# Peaktype vs PRC depletion
#############################
#All memory categories together, diff vs logfc 8 months
PRC_depletion_effect.peaktype <- all_peaks_withPRC_withDARs_memory %>% 
  dplyr::mutate(category = ifelse(is.na(memory_class), "no_dar", glue("{memory_class}_{direction}"))) %>% 
  dplyr::filter(!is.na(PRC_peakID)) %>% 
  dplyr::select(peakType, WT4mon_vs_Mut8mon_LFC, category) 

PRC_depletion_effect.peaktype$category <- factor(PRC_depletion_effect.peaktype$category, levels = categories_order)
ggplot(PRC_depletion_effect.peaktype, aes(x = category, y = WT4mon_vs_Mut8mon_LFC, fill = peakType, colour = peakType)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust =2, trim = TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25) +
  facet_wrap( ~ peakType) +
  geom_hline(yintercept=2, lty = "dashed", size = 0.5) +
  theme_classic() + theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "LFC effect after PRC depletion on PRC-overlapping cocaine-induced DARs ", x= "", y= "WT/Mut-8months (LFC)") 




# ========== 11. Gene integration ==========
degs_simplified <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/6_data_for_Alessandro/DEGs_memory_classes.tsv") %>% 
  dplyr::mutate(degs_category = glue("{memory_class}_{direction}")) %>% dplyr::select(gene, degs_category)

PRC_depletion_effect.onlyDARs.wDEGs <- PRC_depletion_effect %>%
  dplyr::filter(category != "no_dar") %>% 
  left_join(degs_simplified, by = c("nearestGene" = "gene"))


# ========== 12. TF candidate overview (per contrast and per memory class) ==========
prepare_tf_data <- function(data, tf_name, contrasts) {
  data %>%
    dplyr::group_by(PRC_peakID, contrast, category) %>%
    dplyr::summarise(
      TotalPeaks = n(),
      TFPeaks = sum(!!rlang::sym(tf_name), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(PRC_peakID, contrast) %>%
    dplyr::mutate(PercentageTF = round((TFPeaks / TotalPeaks) * 100, 1)) %>%
    dplyr::mutate(PRC_peakID = factor(PRC_peakID, levels = c("no_prc", "prc")),
                  contrast = factor(contrast, levels = contrasts),
                  category = factor(category, levels = c("no_dars", "downreg", "upreg")))
}


# 12A: per contrast: percentage and counts
###################################################
TFmotifs <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/2024/6_data_for_Alessandro/VTA_DN_allPeaks_motifs_matrix_collapsed.rds")
TFmotifs_selected <- as.data.frame(TFmotifs) %>%
  dplyr::select(NR4A1, AP1_var1, AP1_var2) %>% 
  rownames_to_column("peakID") %>% 
  dplyr::mutate(peakID = sub("-", ":", peakID))

all_peaks_withPRC_withDARs_contrast.v2 <- all_peaks_withPRC_withDARs_contrast %>%
  dplyr::filter(contrast %in% all_contrasts[1:5]) %>%
  dplyr::mutate(direction = ifelse(logFC > 0, "upreg", "downreg")) %>%
  dplyr::select(peakID, PRC_peakID, contrast, diff, direction) %>%
  dplyr::mutate(category = case_when(
    diff == "No" ~ "no_dars",
    diff == "Yes" & direction == "upreg" ~ "upreg",
    diff == "Yes" & direction == "downreg" ~ "downreg")) %>% 
  dplyr::mutate(PRC_peakID = ifelse(is.na(PRC_peakID), "no_prc", "prc")) %>% 
  left_join(TFmotifs_selected, by = "peakID") %>% 
  dplyr::select(PRC_peakID, contrast, category, NR4A1, AP1_var1, AP1_var2)


#Select TF: 
# Variable to specify the TF column name, e.g., 'NR4A1'
tf_column <- "AP1_var2"
all_contrasts <- c("h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline", "h24_cocaine-saline", "d14_cocaine-saline")

# Preparing the data for the plot
all_peaks_withPRC_withDARs_TF.toplot <- all_peaks_withPRC_withDARs_contrast.v2 %>% 
  dplyr::group_by(PRC_peakID, contrast, category) %>% 
  dplyr::summarise(
    TotalPeaks = n(), 
    TFPeaks = sum(!!rlang::sym(tf_column), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    contrast = factor(contrast, levels = all_contrasts)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(PRC_peakID, contrast) %>% 
  dplyr::mutate(
    PercentageTF = round((TFPeaks / TotalPeaks) * 100, 1)
  ) %>%
  tidyr::pivot_longer(cols = c(PercentageTF, TFPeaks), names_to = "Metric", values_to = "Value") %>% 
  dplyr::mutate(Metric = factor(Metric, levels = c("TFPeaks", "PercentageTF")))

# Calculating the total number of unique peaks per PRC_peakID category
binding_totals <- all_peaks_withPRC_withDARs_contrast.v2 %>%
  filter(!!sym(tf_column)) %>%  # Using !!sym to dynamically use the current_tf variable
  group_by(PRC_peakID) %>%
  dplyr::summarize(TF_binding_count = n() / 5) 

# Adding this information to facet labels
facet_labels <- binding_totals %>% 
  dplyr::mutate(Label = ifelse(
    PRC_peakID == "no_prc",
    glue("{PRC_peakID}-{tf_column} ATAC peaks ({TF_binding_count} out of 99080)"),
    glue("{PRC_peakID}-{tf_column} ATAC peaks ({TF_binding_count} out of 5465)" )))


# Mapping the facet labels
label_mapping <- setNames(facet_labels$Label, facet_labels$PRC_peakID)

# Plotting the data
ggplot(all_peaks_withPRC_withDARs_TF.toplot, aes(x = contrast, y = Value, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Value),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3,
            check_overlap = TRUE) +
  theme_minimal() +
  scale_fill_manual(values = c("no_dars" = "gray90", 
                               "downreg" = "#EC5D5B", 
                               "upreg" = "#89B8D2")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = glue("Overview of {tf_column} putative binding "),
       x = "Contrast",
       y = "",
       fill = "DARs category") +
  facet_grid(rows = vars(Metric), cols = vars(PRC_peakID), scales = "free_y", labeller = labeller(PRC_peakID = label_mapping))



# 12B: per memory class: percentage and counts
###############################################
TFmotifs <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/2024/6_data_for_Alessandro/VTA_DN_allPeaks_motifs_matrix_collapsed.rds")
TFmotifs_selected <- as.data.frame(TFmotifs) %>%
  dplyr::select(NR4A1, AP1_var1, AP1_var2) %>% 
  rownames_to_column("peakID") %>% 
  dplyr::mutate(peakID = sub("-", ":", peakID))

all_peaks_withPRC_withDARs_memory.v2 <- all_peaks_withPRC_withDARs_memory %>%
  dplyr::select(peakID, PRC_peakID, memory_class, direction) %>%
  dplyr::mutate(PRC_peakID = ifelse(is.na(PRC_peakID), "no_prc", "prc"), 
                category = ifelse(is.na(memory_class), "no_dar", glue("{memory_class}_{direction}"))) %>% 
  left_join(TFmotifs_selected, by = "peakID") %>% 
  dplyr::select(PRC_peakID, category, NR4A1, AP1_var1, AP1_var2) 


#Select TF: 
# Variable to specify the TF column name, e.g., 'NR4A1'
tf_column <- "AP1_var2"
categories_order

# Preparing the data for the plot
all_peaks_withPRC_withDARs_TF.toplot <- all_peaks_withPRC_withDARs_memory.v2 %>% 
  dplyr::group_by(PRC_peakID, category) %>% 
  dplyr::summarise(
    TotalPeaks = n(), 
    TFPeaks = sum(!!rlang::sym(tf_column), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    category = factor(category, levels = categories_order)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(PRC_peakID, category) %>% 
  dplyr::mutate(
    PercentageTF = round((TFPeaks / TotalPeaks) * 100, 1)
  ) %>%
  tidyr::pivot_longer(cols = c(PercentageTF, TFPeaks), names_to = "Metric", values_to = "Value") %>% 
  dplyr::mutate(Metric = factor(Metric, levels = c("TFPeaks", "PercentageTF")))

# Calculating the total number of unique peaks per PRC_peakID category
binding_totals <- all_peaks_withPRC_withDARs_memory.v2 %>%
  filter(!!sym(tf_column)) %>%  # Using !!sym to dynamically use the current_tf variable
  group_by(PRC_peakID) %>%
  dplyr::summarize(TF_binding_count = n()) 

# Adding this information to facet labels
facet_labels <- binding_totals %>% 
  dplyr::mutate(Label = ifelse(
    PRC_peakID == "no_prc",
    glue("{PRC_peakID}-{tf_column} ATAC peaks ({TF_binding_count} out of 99080)"),
    glue("{PRC_peakID}-{tf_column} ATAC peaks ({TF_binding_count} out of 5465)" )))


# Mapping the facet labels
label_mapping <- setNames(facet_labels$Label, facet_labels$PRC_peakID)

# Plotting the data
ggplot(all_peaks_withPRC_withDARs_TF.toplot, aes(x = category, y = Value, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Value),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3,
            check_overlap = TRUE) +
  theme_minimal() +
  scale_fill_manual(values = memory_classes_colors)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = glue("Overview of {tf_column} putative binding "),
       x = "Contrast",
       y = "",
       fill = "DARs category") +
  facet_grid(rows = vars(Metric), cols = vars(PRC_peakID), scales = "free_y", labeller = labeller(PRC_peakID = label_mapping))



# ========== 13. TF enrichment analysis: PRC vs non-PRC ==========
pval = 0.05
logfc = 0.25
#Whicn bdg? I think for the first question would be the whole list of DARs?
wnn.seu.VTA <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/6_231107_updated_analyses/2_ATAC_alternative_processing/3_231205_ArchR_DARsandMotifs/231219_rpkm_based_DARs/240119_thresholds_setC/240122.wnn.seu.VTA.withMotifs.rds")
bdg_allDARs <- VTA_RPKM_formatted_metadata_withDARs %>% dplyr::filter(complete.cases(.), p_val < pval, abs(logFC) >=logfc) %>% dplyr::select(peakID_v2) %>% distinct() %>% pull() #VTA_RPKM_formatted_metadata_withDARs from: 240315_ATAC_master_table.R !!


# Add motif information to seurat object
pfm <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/1_ATAC_masterTable_continued/231217.pfm.rds")
wnn.seu.VTA <- AddMotifs(
  object = wnn.seu.VTA,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

TFenrichment.complete.df <- data.frame()

for (current_category in categories_order[2:9]) {
  message(glue("Analyzing: {current_category}"))
  
  contrast.group1.dars <- all_peaks_withPRC_withDARs_memory %>%
    dplyr::mutate(category = ifelse(is.na(memory_class), "no_dar", glue("{memory_class}_{direction}")), 
                  peakID_v2 = glue("{seqnames}-{start}-{end}")) %>% 
    dplyr::filter(category == current_category, !is.na(PRC_peakID)) %>% 
    dplyr::select(peakID_v2) %>% pull() %>% as.character()
  
  #B.3. TFenrichment; bdg: allDARs
  contrast.group1vsallDARs.tfenrich <- FindMotifs(
    object = wnn.seu.VTA,
    features = contrast.group1.dars, 
    background = bdg_allDARs) %>% 
    dplyr::mutate(category = current_category, bdg = "allDARs")
  
  TFenrichment.complete.df <-rbind(TFenrichment.complete.df, contrast.group1vsallDARs.tfenrich)
}

#write_tsv(TFenrichment.complete.df, "241028.PRCpeaks_TFenrichment.complete.df.tsv")
