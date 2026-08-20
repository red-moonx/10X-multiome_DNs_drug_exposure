# ===========================================
# Script Title: 5.2 DN identity exploration II
# Author: Luna Zea Redondo
# Date: 2026-05-26
# Description:
#   This workflow evaluates the expression of DN identity genes using orthogonal data
#   from Phillips et al. 2022. Explores wilcoxon test with different backgrounds
# ===========================================

# ===========================================
# 1. Load inputs and define gene sets
# ===========================================
DEG_complete_results <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/3_muscat/deg_results/260513_DEG_complete_results_Analysis1.tsv")
seu.VTA_DNs <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/2_DN-evaluation/250507_seu.VTA_DNs.rds")

DN_markers <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/phillip.tsv") %>%
  dplyr::select(gene, p.adj, log2FC) %>%
  dplyr::filter(gene %in% DEG_complete_results$gene)

top_25_genes  <- DN_markers %>% arrange(p.adj, desc(log2FC)) %>% slice_head(n = 25)
top_50_genes  <- DN_markers %>% arrange(p.adj, desc(log2FC)) %>% slice_head(n = 50)
top_100_genes <- DN_markers %>% arrange(p.adj, desc(log2FC)) %>% slice_head(n = 100)
top_200_genes <- DN_markers %>% arrange(p.adj, desc(log2FC)) %>% slice_head(n = 200)
top_300_genes <- DN_markers %>% arrange(p.adj, desc(log2FC)) %>% slice_head(n = 300)

# ===========================================
# 2. Define background: all expressed genes
#    in DNs, excluding top 300 DN markers
# ===========================================
# expressed_genes <- rownames(seu.VTA_DNs)[
#   rowSums(GetAssayData(seu.VTA_DNs, slot = "counts") > 0) > 0
# ]
detection_counts <- rowSums(GetAssayData(seu.VTA_DNs, assay = "RNA", slot = "counts") > 0)
genes_detected_10cells <- names(detection_counts[detection_counts >= 10])

#expressed_genes <- rownames(seu.VTA_DNs)
expressed_genes <- unique(DEG_complete_results$gene)
background_genes <- expressed_genes[!expressed_genes %in% top_300_genes$gene] 

# ===========================================
# 3. Build plot.data
# ===========================================
DEG_identity_toPlot <- DEG_complete_results %>%
  dplyr::select(gene, logFC, p_val, contrast, query, control, significant)

plot.data <- DEG_identity_toPlot %>%
  mutate(
    in_top_25  = ifelse(gene %in% top_25_genes$gene,  "yes", "no"),
    in_top_50  = ifelse(gene %in% top_50_genes$gene,  "yes", "no"),
    in_top_100 = ifelse(gene %in% top_100_genes$gene, "yes", "no"),
    in_top_200 = ifelse(gene %in% top_200_genes$gene, "yes", "no"),
#    in_top_300 = ifelse(gene %in% top_300_genes$gene, "yes", "no"),
    background = ifelse(gene %in% background_genes,   "yes", "no")
  ) %>%
  dplyr::filter(contrast %in% all_contrasts) %>%
  dplyr::mutate(contrast = factor(contrast, levels = all_contrasts)) %>%
  dplyr::select(gene, logFC, contrast, in_top_25, in_top_50, in_top_100,
                in_top_200, background)

# ===========================================
# 4. Median data and line plot
# ===========================================
group_order  <- c("in_top_25", "in_top_50", "in_top_100", "in_top_200", "background")
group_labels <- c("Top 25", "Top 50", "Top 100", "Top 200", "Background")

median_data <- plot.data %>%
  pivot_longer(cols = c(starts_with("in_top_"), "background"),
               names_to = "group", values_to = "is_yes") %>%
  dplyr::filter(is_yes == "yes") %>%
  dplyr::group_by(group, contrast) %>%
  dplyr::summarise(median_logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  mutate(group = factor(group, levels = group_order, labels = group_labels))

median_line_plot <- ggplot(median_data, aes(x = contrast, y = median_logFC, color = group, group = group)) +
#  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    "Background" = "#7A7C85",
    "Top 25"     = "#9A4A4A",
    "Top 50"     = "#E0BD5A",
    "Top 100"    = "#83AC7E",
    "Top 200"    = "#3A64A3")) +
  coord_cartesian(ylim = c(-0.3, 0.4)) +
  theme_classic() +
  labs(title = "Median logFC of DN identity gene sets across contrasts",
       x = "Contrast", y = "Median logFC", color = "Gene set") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_name <- "chapter5_DNidentity_median_linePlot.pdf"
dev.size()
pdf(plot_name, width = 8.15, height =6.63)
median_line_plot
dev.off()

# ===========================================
# 5. Wilcoxon test 1: against mu = 0
# "Is there repression of DN identity genes?"
# ===========================================

format_significance <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "n.s",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "n.s"
  )
}

wilcox_vs_zero <- plot.data %>%
  pivot_longer(cols = starts_with("in_top_"),
               names_to = "group", values_to = "is_yes") %>%
  dplyr::filter(is_yes == "yes") %>%
  dplyr::group_by(group, contrast) %>%
  dplyr::summarise(
    n_genes         = n(),
    median_logFC    = median(logFC, na.rm = TRUE),
    p_value_vs_zero = wilcox.test(logFC, mu = 0, alternative = "less", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(group = factor(group, levels = group_order, labels = group_labels)) %>%
  arrange(contrast, group) %>% 
  mutate(is_significant = format_significance(p_value_vs_zero))

# ===========================================
# 6. Wilcoxon test 2: against background
# "Is repression specific to DN identity genes?"
# ===========================================
background_logFC <- plot.data %>%
  dplyr::filter(background == "yes") %>%
  dplyr::select(contrast, logFC, gene)

wilcox_vs_background <- plot.data %>%
  pivot_longer(cols = starts_with("in_top_"),
               names_to = "group", values_to = "is_yes") %>%
  dplyr::filter(is_yes == "yes") %>%
  dplyr::group_by(group, contrast) %>%
  dplyr::summarise(
    n_genes                 = n(),
    median_logFC            = median(logFC, na.rm = TRUE),
    median_logFC_background = median(
      background_logFC$logFC[background_logFC$contrast == cur_group()$contrast],
      na.rm = TRUE),
    p_value_vs_background = wilcox.test(
      x           = logFC,
      y           = background_logFC$logFC[background_logFC$contrast == cur_group()$contrast],
      alternative = "less",
      exact       = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(group = factor(group, levels = group_order, labels = group_labels)) %>%
  arrange(contrast, group) %>% 
  mutate(is_significant = format_significance(p_value_vs_background))

# ===========================================
# 7. Combined results table
# ===========================================
wilcox_combined <- wilcox_vs_zero %>%
  dplyr::select(group, contrast, n_genes, median_logFC, p_value_vs_zero) %>%
  left_join(
    wilcox_vs_background %>% dplyr::select(group, contrast, median_logFC_background, p_value_vs_background),
    by = c("group", "contrast")
  )
write_tsv(wilcox_combined, "260607_wilcox_combined.tsv")



# ===============
# Extra. 8h cocaine data
# ===============
#Test: results for 8h cocaine data
median(res8$logFC, na.rm = TRUE)
mean(res8$logFC, na.rm = TRUE)

mean(res8$logFC < 0, na.rm = TRUE)
quantile(res8$logFC, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)

wilcox.test(res8$logFC, mu = 0)
binom.test(sum(res8$logFC < 0, na.rm = TRUE),
           sum(!is.na(res8$logFC)),
           p = 0.5,
           alternative = "greater")


res8 <- DEG_identity_toPlot %>% dplyr::filter(contrast == "h8_cocaine-saline")

res <- res8

res <- res %>%
  dplyr::filter(is.finite(logFC))
res <- res %>%
  dplyr::filter(baseMean > 10)
is.na(logFC)
is.nan(logFC)
is.infinite(logFC)
