# ============================================================
# Script Title: 3.9. ATAC-level QC panels for complete_peaks_table
# Author: Luna Zea Redondo
# Date: 2026-05-13
# Description:
#    This script generates quality control and distribution panels 
#    for ATAC-seq peaks (complete_peaks_table). It evaluates and visualizes 
#    GC content percentages, peak lengths, and contrast-specific RPKM values 
#    across differentially accessible regions (DARs) using faceted boxplots with 
#    median annotations, exporting the resulting figures to PDF format.
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

pdf_dir <- dir

width_original  <- 11.854167
height_original <- 4.666667

direction_colors <- c(
  no_dar  = "gray80",
  downreg = "#4B8EBD",
  upreg   = "#E14C4B"
)

rpkm_group_colors <- c(
  saline       = "gray80",
  h1_cocaine   = "#617641",
  h4_cocaine   = "#C48208",
  h8_cocaine   = "#326186",
  h24_cocaine  = "#AE430A",
  d14_cocaine  = "#564686"
)

# ------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------

complete_peaks_table_qc <- complete_peaks_table %>%
  dplyr::filter(logCPM >= logcpm_threshold) %>%
  tidyr::separate(
    contrast,
    into = c("group1", "group2"),
    sep = "-",
    remove = FALSE
  ) %>%
  dplyr::mutate(
    contrast = factor(contrast, levels = all_contrasts),
    GC_percent = ifelse(GC > 1, GC / 100, GC),
    direction = dplyr::case_when(
      p_val < pval_threshold & logFC >=  logfc_threshold ~ "upreg",
      p_val < pval_threshold & logFC <= -logfc_threshold ~ "downreg",
      TRUE ~ "no_dar"
    ),
    direction = factor(direction, levels = c("no_dar", "downreg", "upreg"))
  )

# ------------------------------------------------------------
# Long data for GC content and peak length
# ------------------------------------------------------------

qc_static_long <- complete_peaks_table_qc %>%
  dplyr::select(
    peakID,
    contrast,
    direction,
    GC_percent,
    peak_length
  ) %>%
  tidyr::pivot_longer(
    cols = c(GC_percent, peak_length),
    names_to = "variable",
    values_to = "value"
  )

# ------------------------------------------------------------
# Long data for contrast-specific RPKM
# ------------------------------------------------------------

rpkm_by_contrast_long <- complete_peaks_table_qc %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("rpkm_"),
    names_to = "rpkm_group",
    values_to = "rpkm_value",
    names_prefix = "rpkm_"
  ) %>%
  dplyr::filter(rpkm_group == group1 | rpkm_group == group2) %>%
  dplyr::filter(rpkm_value > 0) %>%
  dplyr::mutate(
    rpkm_group = factor(
      rpkm_group,
      levels = c("saline", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")
    )
  )

# ------------------------------------------------------------
# Function for GC and peak length panels
# ------------------------------------------------------------

make_static_qc_panel <- function(var, title, ylab, log_y = FALSE) {
  
  plot_data <- qc_static_long %>%
    dplyr::filter(variable == var)
  
  med <- plot_data %>%
    dplyr::group_by(contrast, direction) %>%
    dplyr::summarise(
      median_value = round(stats::median(value, na.rm = TRUE), 3),
      .groups = "drop"
    )
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = direction, y = value, fill = direction)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    ggplot2::stat_summary(
      fun = stats::median,
      geom = "point",
      size = 2,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = med,
      ggplot2::aes(x = direction, y = median_value, label = median_value),
      inherit.aes = FALSE,
      size = 3,
      color = "black",
      vjust = -0.5
    ) +
    ggplot2::facet_wrap(~contrast, nrow = 1) +
    ggplot2::scale_fill_manual(values = direction_colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = title,
      x = "DAR Category",
      y = ylab,
      fill = "Direction"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  if (log_y) {
    p <- p + ggplot2::scale_y_log10()
  }
  
  p
}

# ------------------------------------------------------------
# GC content panel
# ------------------------------------------------------------

gc_panel <- make_static_qc_panel(
  var = "GC_percent",
  title = "GC Content Distribution",
  ylab = "GC Content",
  log_y = FALSE
) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  )

# ------------------------------------------------------------
# Peak length panel
# ------------------------------------------------------------

length_panel <- make_static_qc_panel(
  var = "peak_length",
  title = "Peak Length Distribution",
  ylab = "Peak Length",
  log_y = TRUE
)

# ------------------------------------------------------------
# RPKM panel
# ------------------------------------------------------------

rpkm_medians <- rpkm_by_contrast_long %>%
  dplyr::group_by(contrast, direction, rpkm_group) %>%
  dplyr::summarise(
    median_value = round(stats::median(rpkm_value, na.rm = TRUE), 3),
    .groups = "drop"
  )

rpkm_panel <- ggplot2::ggplot(
  rpkm_by_contrast_long,
  ggplot2::aes(x = direction, y = rpkm_value, fill = rpkm_group)
) +
  ggplot2::geom_boxplot(
    outlier.shape = NA,
    alpha = 0.8,
    position = ggplot2::position_dodge(width = 0.8)
  ) +
  ggplot2::geom_text(
    data = rpkm_medians,
    ggplot2::aes(
      x = direction,
      y = median_value,
      label = median_value,
      group = rpkm_group
    ),
    inherit.aes = FALSE,
    position = ggplot2::position_dodge(width = 0.8),
    size = 3,
    color = "black",
    vjust = -0.5
  ) +
  ggplot2::facet_wrap(~contrast, nrow = 1) +
  ggplot2::scale_fill_manual(values = rpkm_group_colors) +
  ggplot2::scale_y_log10() +
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = "RPKM Distribution by Contrast Group",
    x = "DAR Category",
    y = "RPKM",
    fill = "Group"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

# ------------------------------------------------------------
# Print panels
# ------------------------------------------------------------

print(gc_panel)
print(rpkm_panel)
print(length_panel)

# ------------------------------------------------------------
# Save panels as PDF
# ------------------------------------------------------------

pdf(
  file = glue("{pdf_dir}/chapter3_ATACqc_GC_panel.pdf"),
  width = width_original,
  height = height_original
)
gc_panel
dev.off()

pdf(
  file = glue("{pdf_dir}/chapter3_ATACqc_RPKM_panel.pdf"),
  width = width_original*2,
  height = height_original
)
rpkm_panel
dev.off()

pdf(
  file = glue("{pdf_dir}/chapter3_ATACqc_length_panel.pdf"),
  width = width_original,
  height = height_original
)
length_panel
dev.off()
