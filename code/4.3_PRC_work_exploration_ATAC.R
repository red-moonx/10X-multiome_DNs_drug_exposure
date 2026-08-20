# ===========================================
# Script Title: 4.3_PRC_work_exploration_ATAC.R
# Author: Luna Zea Redondo
# Date: 2026-05-21
#
# Description:
#   Explores overlap between PRC2-marked regions and ATAC peaks/DARs.
#   Main steps:
#     - load PRC2-annotated ATAC tables
#     - quantify PRC overlap across DAR contrasts and memory classes
#     - test PRC enrichment/depletion against no-DAR background
#     - generate MEME-CHIP BED/FASTA files
#     - explore peak type and PRC depletion effects
# ===========================================

#!/usr/bin/env Rscript

rm(list = ls(all.names = TRUE)); gc()
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(glue)
library(ggplot2)
library(scales)
library(patchwork)
library(GenomicRanges)
library(IRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)

setwd("/fast/AG_Pombo/luna/2026_rebuttal/16_DARs_vs_PRC2")

# ========== Parameters ==========
all_contrasts <- c("h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline", "h24_cocaine-saline", "d14_cocaine-saline",
                   "h4_cocaine-h1_cocaine", "h8_cocaine-h4_cocaine", "h24_cocaine-h8_cocaine", "d14_cocaine-h24_cocaine")

categories_order <- c("no_dar", "transient_up", "transient_down", "recovered_up", "recovered_down",
                      "memory_up", "memory_down", "delayed_up", "delayed_down")

memory_classes_colors <- c("no_dar" = "gray90", "transient_up" = "#D9481C", "recovered_up" = "#F3A712",
                           "memory_up" = "#008AB8", "delayed_up" = "#954B77",
                           "transient_down" = "#EF9A81", "recovered_down" = "#F8CD77",
                           "memory_down" = "#8FD0E6", "delayed_down" = "#C892B3")

peak_type_colors <- c(Promoter = "#8B0000", Exonic = "#DAA520", Intronic = "#2E8B57", Distal = "#4682B4")

# ========== 1. Load input tables ==========
all_peaks_withPRC <- read_tsv("260615_all_peaks_withPRC.tsv")
all_peaks_withPRC_withDARs_contrast <- read_tsv("260615_all_peaks_withPRC_withDARs_contrast.tsv")
all_peaks_withPRC_withDARs_memory <- read_tsv("260615_all_peaks_withPRC_withDARs_memory.tsv")

# ========== 2. PRC peak overview ==========
PRC_length.df <- all_peaks_withPRC %>%
  dplyr::filter(!is.na(PRC_peakID)) %>%
  mutate(PRC_length = as.numeric(str_extract(PRC_peakID, "(?<=-)[0-9]+")) -
           as.numeric(str_extract(PRC_peakID, "(?<=:)[0-9]+")))

PRC_length.plot <- ggplot(PRC_length.df, aes(PRC_length)) +
  geom_density(fill = "black", alpha = 0.5) + scale_x_log10() +
  labs(title = "Density of PRC peak lengths", x = "PRC peak length (log10)", y = "Density") +
  theme_minimal()

median_length <- PRC_length.df %>% summarise(median_length = median(PRC_length, na.rm = TRUE))

PRC_LFC.plot.df <- all_peaks_withPRC %>%
  dplyr::filter(!is.na(PRC_peakID)) %>%
  distinct(PRC_peakID, .keep_all = TRUE) %>%
  dplyr::select(WT4mon_vs_Mut4mon_diff, WT4mon_vs_Mut4mon_LFC)

PRC_LFC.plot <- ggplot(PRC_LFC.plot.df, aes(WT4mon_vs_Mut4mon_LFC, WT4mon_vs_Mut4mon_diff)) +
  geom_point(alpha = 0.5) + geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  xlim(-3, 6) + labs(x = "WT4mon vs Mut4mon LFC", y = "WT4mon vs Mut4mon diff", title = "PRC peaks") +
  theme_classic()

# ========== 3. PRC overlap per contrast ==========
PRC_contrast.df <- all_peaks_withPRC_withDARs_contrast %>%
  dplyr::filter(contrast %in% all_contrasts[1:5]) %>%
  mutate(direction = ifelse(logFC > 0, "upreg", "downreg"),
         category = case_when(diff == "No" ~ "no_dars", diff == "Yes" & direction == "upreg" ~ "upreg",
                              diff == "Yes" & direction == "downreg" ~ "downreg"),
         PRC_present = !is.na(PRC_peakID)) %>%
  dplyr::select(peakID, PRC_peakID, contrast, diff, direction, category, PRC_present)

PRC_contrast.summary <- PRC_contrast.df %>%
  dplyr::group_by(contrast, category) %>%
  dplyr::summarise(total_peakIDs = n(), true_peakIDs = sum(PRC_present), .groups = "drop") %>%
  mutate(percentage = true_peakIDs / total_peakIDs * 100,
         contrast = factor(contrast, levels = all_contrasts[1:5]),
         category = factor(category, levels = c("no_dars", "downreg", "upreg")))

PRC_contrast.plot <- ggplot(PRC_contrast.summary, aes(contrast, percentage, fill = category)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%% (%d)", percentage, true_peakIDs)),
            position = position_dodge(width = 0.9), vjust = -0.5) +
  scale_fill_manual(values = c("no_dars" = "gray90", "upreg" = "#EC5D5B", "downreg" = "#89B8D2")) +
  labs(title = "PRC-positive peaks per contrast", x = "Contrast", y = "% PRC-positive peaks", fill = "Category") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ========== 4. PRC overlap: DARs vs non-DARs ==========
PRC_direction.df <- all_peaks_withPRC_withDARs_memory %>%
  mutate(category = ifelse(is.na(memory_class), "no_dar", direction),
         PRC_state = ifelse(is.na(PRC_peakID), "no_prc", "prc"),
         category = factor(category, levels = c("no_dar", "down", "up"))) %>%
  dplyr::count(category, PRC_state, name = "count") %>%
  dplyr::group_by(category) %>% mutate(percentage = count / sum(count) * 100) %>% ungroup()

PRC_direction.plot <- PRC_direction.df %>%
  dplyr::filter(PRC_state == "prc") %>%
  ggplot(aes(category, percentage, fill = category)) +
  geom_col(color = "black") +
  geom_text(aes(label = sprintf("%.1f%% (%d)", percentage, count)), vjust = -0.4) +
  scale_fill_manual(values = c("no_dar" = "gray90", "up" = "#EC5D5B", "down" = "#89B8D2")) +
  labs(title = "PRC-positive peaks in DARs vs non-DARs", x = "", y = "% PRC-positive") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

PRC_peaks.combined <- PRC_LFC.plot + PRC_direction.plot

PRC_fisher_results <- data.frame(
  comparison = c("down_vs_no_dar", "up_vs_no_dar"),
  prc_cat = c(229, 344),
  no_prc_cat = c(3056, 3250),
  prc_no_dar = 3337,
  no_prc_no_dar = 53833
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    fisher = list(fisher.test(
      matrix(c(prc_cat, no_prc_cat, prc_no_dar, no_prc_no_dar),
             nrow = 2, byrow = TRUE,
             dimnames = list(c("category", "no_dar"), c("prc", "no_prc"))),
      alternative = "greater"
    )),
    odds_ratio = unname(fisher$estimate),
    p_value = fisher$p.value
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    sig_label = dplyr::case_when(
      p_adj_BH < 0.001 ~ "FDR < 0.001 (***)",
      p_adj_BH < 0.005 ~ "FDR < 0.005 (**)",
      p_adj_BH < 0.05  ~ "FDR < 0.05 (*)",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(comparison, odds_ratio, p_value, p_adj_BH, sig_label)

PRC_fisher_results


# ========== 5. PRC-positive peaks per memory class ==========
PRC_memory.df <- all_peaks_withPRC_withDARs_memory %>%
  mutate(category = ifelse(is.na(memory_class), "no_dar", glue("{memory_class}_{direction}")),
         PRC_state = ifelse(is.na(PRC_peakID), "no_prc", "prc")) %>%
  dplyr::count(category, PRC_state, name = "count") %>%
  group_by(category) %>% mutate(percentage = count / sum(count) * 100) %>% ungroup() %>%
  mutate(category = factor(category, levels = categories_order)) %>%
  dplyr::filter(PRC_state == "prc")

PRC_memory.barplot <- ggplot(PRC_memory.df, aes(category, percentage, fill = category)) +
  geom_col(color = "black") +
  geom_text(aes(label = paste0(round(percentage, 1), "% (", count, ")")), vjust = -0.5, size = 4) +
  scale_fill_manual(values = memory_classes_colors) +
  labs(title = "PRC-positive peaks per memory class", x = "", y = "% PRC-positive") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

PRC_memory.mirror.df <- PRC_memory.df %>%
  mutate(category_group = gsub("_up|_down", "", category),
         percentage_adjusted = ifelse(grepl("down$", category), -percentage, percentage),
         regulation = case_when(grepl("up$", category) ~ "Up", grepl("down$", category) ~ "Down"),
         vjust = ifelse(percentage_adjusted > 0, -0.5, 1.5),
         category_group = factor(category_group, levels = c("transient", "recovered", "memory", "delayed")),
         regulation = factor(regulation, levels = c("Up", "Down")))

no_dar_value <- PRC_memory.df$percentage[PRC_memory.df$category == "no_dar"]

PRC_memory.mirror.plot <- PRC_memory.mirror.df %>%
  dplyr::filter(!is.na(regulation)) %>%
  ggplot(aes(category_group, percentage_adjusted, fill = regulation)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = paste0(round(abs(percentage), 1), "% (", count, ")"), vjust = vjust)) +
  geom_hline(yintercept = c(no_dar_value, -no_dar_value), linetype = "dashed", color = "goldenrod", linewidth = 1) +
  scale_fill_manual(values = c("Up" = "#EC5D5B", "Down" = "#89B8D2")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Percentage of PRC-DARs per memory class", x = "Memory class", y = "% PRC-positive", fill = "DAR") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot_name <- "chapter4_PRC_memoryDARss_barplot.pdf"
# dev.size()
# pdf(plot_name, width = 5.45, height =5.75)
# PRC_memory.mirror.plot
# dev.off()

# ========== 6. Empirical permutation test against no_dar baseline ==========
PRC_memory.p_emp.df <- all_peaks_withPRC_withDARs_memory %>%
  mutate(category = ifelse(is.na(memory_class), "no_dar", glue("{memory_class}_{direction}")),
         PRC_state = ifelse(is.na(PRC_peakID), "no_prc", "prc")) %>%
  dplyr::select(peakID, category, PRC_state)

summary_emp <- PRC_memory.p_emp.df %>%
  dplyr::group_by(category) %>%
  dplyr::summarise(total_count = n(), prc_count = sum(PRC_state == "prc"), .groups = "drop") %>%
  mutate(percentage = prc_count / total_count * 100,
         category = factor(category, levels = categories_order)) %>%
  dplyr::arrange(category) %>% dplyr::filter(!is.na(category))

baseline_data <- PRC_memory.p_emp.df %>% dplyr::filter(category == "no_dar")
num_its <- 10000; plot_list <- list(); emp_results <- data.frame()

for (i in seq_len(nrow(summary_emp))) {
  cat <- as.character(summary_emp$category[i])
  total_count <- summary_emp$total_count[i]
  obs_count <- summary_emp$prc_count[i]
  
  message(glue("Running: {cat}: {total_count} total peaks and {obs_count} PRC targets"))
  
  obs <- replicate(num_its, baseline_data %>% sample_n(size = total_count) %>% dplyr::filter(PRC_state == "prc") %>% nrow())
  p_emp <- (min(sum(obs >= obs_count), sum(obs <= obs_count)) + 1) / (num_its + 1)
  emp_results <- rbind(emp_results, data.frame(category = cat, total_count, obs_count, p_emp))
  
  plot_list[[cat]] <- ggplot(data.frame(value = obs), aes(value)) +
    geom_density(fill = "gray90", alpha = 0.5) +
    geom_vline(xintercept = obs_count, color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = glue("{cat}; p = {round(p_emp, 4)}"), x = "", y = "") +
    theme_classic() + theme(plot.title = element_text(size = 10))
}

panel_perm <- plot_list[c("transient_up", "recovered_up", "memory_up", "delayed_up",
                          "transient_down", "recovered_down", "memory_down", "delayed_down")] %>%
  Reduce(`+`, .) + plot_layout(ncol = 4, nrow = 2)

# pdf("chapter4_PRC_memoryDARs_permutations.pdf", width = 11.28, height = 7.17)
# panel_perm
# dev.off()

# ========== 7. Chi-square tests against no_dar ==========
PRC_dars_unified_df <- all_peaks_withPRC_withDARs_memory %>%
  mutate(category = ifelse(is.na(memory_class), "no_dar", glue("{memory_class}_{direction}")),
         prc_state = ifelse(is.na(PRC_peakID), "no_prc", "prc")) %>%
  distinct(peakID, .keep_all = TRUE) %>%
  dplyr::select(peakID, category, prc_state)

p_values <- list(); contingency_tables <- list()

for (cat in categories_order[2:9]) {
  filtered_data <- PRC_dars_unified_df %>%
    dplyr::filter(category %in% c("no_dar", cat))
  contingency_table <- table(filtered_data$prc_state, filtered_data$category == cat)
  chi_square_test <- chisq.test(contingency_table)
  p_values[[cat]] <- chi_square_test$p.value
  contingency_tables[[cat]] <- as.data.frame.matrix(contingency_table)
}

p_values_df <- data.frame(category = names(p_values), p_value = unlist(p_values)) %>%
  mutate(p_adj_bonferroni = p.adjust(p_value, method = "bonferroni"))

contingency_tables_df <- bind_rows(contingency_tables, .id = "category")

# ========== 8. Peak type distribution by PRC status and accessibility ==========
#Add peaktype
VTA_complete_peakset <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/10_DARs_TFenrichment-contrast/260610_VTA-DNs_complete_peakset.rds")

VTA_peaktype <- VTA_complete_peakset %>%
  dplyr::select(peakID, peakType) %>%
  dplyr::distinct() %>%
  dplyr::mutate(peakID = sub("^(chr[^-]+)-", "\\1:", peakID))

PRC_peaktype.df <- all_peaks_withPRC_withDARs_memory %>%
  dplyr::left_join(VTA_peaktype, by = "peakID") %>%
  dplyr::mutate(category = ifelse(is.na(memory_class), "no_dar", direction),
                PRC_state = ifelse(is.na(PRC_peakID), "no_prc", "prc"),
                category = factor(category, levels = c("no_dar", "down", "up")),
                PRC_state = factor(PRC_state, levels = c("no_prc", "prc"))) %>%
  dplyr::count(category, PRC_state, peakType, name = "count") %>%
  dplyr::group_by(category, PRC_state) %>%
  dplyr::mutate(percentage = count / sum(count) * 100) %>%
  dplyr::ungroup()

PRC_peaktype.df <- PRC_peaktype.df %>%
  dplyr::mutate(PRC_state = factor(PRC_state, levels = c("prc", "no_prc")))

PRC_peaktype.df <- PRC_peaktype.df %>%
  dplyr::mutate(
    category_label = dplyr::case_when(
      category == "no_dar" ~ "non-DAR",
      category == "down" ~ "Less\naccessible",
      category == "up" ~ "More\naccessible"
    ),
    category_label = factor(category_label, levels = c("non-DAR", "Less\naccessible", "More\naccessible")),
    PRC_state = factor(PRC_state, levels = c("no_prc", "prc"), labels = c("No", "Yes"))
  )

PRC_peaktype_paper <- ggplot(PRC_peaktype.df, aes(x = percentage, y = PRC_state, fill = peakType)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_fill_manual(values = peak_type_colors) +
  facet_grid(category_label ~ ., switch = "y") +
  labs(x = "Genomic distribution of H3K27me3-marked ATAC peaks\nacross differential accessibility classes",
       y = "H3K27me3\nmarked",
       fill = "Genomic\ndistribution") +
  theme_minimal() +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, size = 11),
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    panel.spacing.y = unit(0.8, "lines"),
    legend.position = "right"
  )

PRC_peaktype_paper

# dev.size()
# width_original = 8.25
# height_original= 3.73
# 
# plot_name <- "chapter4_DAR_PRC2_peaktype.pdf"
# 
# pdf(plot_name, width = width_original, height =height_original)
# PRC_peaktype_paper
# dev.off()

#load("PRC2_vs_ATAC.rds")
