#!/usr/bin/env Rscript

# ==============================================================================
# Script: 3.8_prepare_STREME_inputs_and_GCcontent.R
# Author: Luna Zea Redondo
# Updated: 2026-06-17
#
#    This script prepares inputs for motif discovery (STREME) by generating 
#    FASTA sequences for dynamic ATAC-seq peak categories alongside a shared 
#    non-DAR background and export manifests. Additionally, it computes and 
#    visualizes sequence composition metrics (GC content, CpG/CpH dinucleotides, 
#    and CGCG tetranucleotides) across categories using randomized non-DAR control sets.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------

rm(list = ls(all.names = TRUE))
gc()

set.seed(1)

library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)

`%notin%` <- Negate(`%in%`)

# ------------------------------------------------------------------------------
# 1. set working environment parameters
# ------------------------------------------------------------------------------

dir <- "/fast/AG_Pombo/luna/2026_rebuttal/13_STREME_and_GCcontent"
setwd(dir)

# Load input data:
DAR_complete_results_kmeans <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/11_DARs-TRaCE/260612_DAR_complete_results_kmeans.tsv")
VTA_complete_peakset <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/10_DARs_TFenrichment-contrast/260610_VTA-DNs_complete_peakset.rds")

# Background for STREME.

# Recommended default:
streme_background_mode <- "non_DARs"

# Random non-DAR sets used only for GC plots.
# These are plotted as sample_200 / sample_800 / sample_2000 controls.
random_control_sizes <- c(200, 800, 2000)

# Category order for plots.
category_levels <- c(
  "sample_200", "sample_800", "sample_2000",
  "transient_down", "transient_up",
  "recovered_down", "recovered_up",
  "memory_down", "memory_up",
  "delayed_down", "delayed_up"
)

# Colors for plots.
category_colors <- c("transient_up" = "#D9481C", "recovered_up" = "#F3A712", "memory_up"  = "#008AB8", "delayed_up" = "#954B77",
                     "transient_down" = "#EF9A81", "recovered_down" = "#F8CD77", "memory_down" = "#8FD0E6", "delayed_down" = "#C892B3",
                     "sample_200" = "gray70", "sample_800" = "gray70", "sample_2000" = "gray70")


# ------------------------------------------------------------------------------
# Output directories
# ------------------------------------------------------------------------------

fasta_dir <- file.path(dir, "01_FASTA_for_STREME")
category_fasta_dir <- file.path(fasta_dir, "categories")
background_fasta_dir <- file.path(fasta_dir, "background")
per_category_background_dir <- file.path(fasta_dir, "per_category_backgrounds")
region_fasta_dir <- file.path(fasta_dir, "per_region_optional")

gc_dir <- file.path(dir, "02_GC_content")
plot_dir <- file.path(gc_dir, "plots")
table_dir <- file.path(gc_dir, "tables")

for (x in c(
  fasta_dir,
  category_fasta_dir,
  background_fasta_dir,
  per_category_background_dir,
  region_fasta_dir,
  gc_dir,
  plot_dir,
  table_dir
)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

###############################################
# 1. Prepare DAR table for STREME
###############################################

DARs_for_STREME <- DAR_complete_results_kmeans %>%
  dplyr::mutate(category_name = glue("{memory_class}_{direction}")) %>%
  dplyr::select(peakID, category_name) %>%
  dplyr::distinct()

unique_categories <- unique(DARs_for_STREME$category_name)

table(DARs_for_STREME$category_name)


###############################################
# 2. Generate FASTA files for each DAR category
###############################################
for (category in unique_categories) {
  
  subset_data <- DARs_for_STREME %>%
    dplyr::filter(category_name == category) %>%
    dplyr::distinct(peakID, .keep_all = TRUE)
  
  subset_data.gr <- GRanges(subset_data$peakID)
  
  input_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, subset_data.gr)
  names(input_sequences) <- subset_data$peakID
  
  writeXStringSet(
    input_sequences,
    filepath = glue("{category_fasta_dir}/{category}.fa")
  )
}


###############################################
# 3. Generate shared non-DAR background FASTA
###############################################

# VTA_complete_peakset has:
#   peakID    = chr1-3184958-3185226
#   peakID_v2 = chr1:3184958-3185226
#
# DAR_complete_results_kmeans$peakID uses the peakID_v2 format.
# Therefore, use peakID_v2 for matching and GRanges.

background_noDARs <- VTA_complete_peakset %>%
  dplyr::filter(peakID_v2 %notin% DARs_for_STREME$peakID) %>%
  dplyr::select(peakID = peakID_v2) %>%
  dplyr::distinct()

background_noDARs.gr <- GRanges(background_noDARs$peakID)

background_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, background_noDARs.gr)
names(background_sequences) <- background_noDARs$peakID

writeXStringSet(
  background_sequences,
  filepath = glue("{background_fasta_dir}/background_nonDARs.fa")
)


###############################################
# 4. Optional: all-other-DAR backgrounds
###############################################
# These are optional alternative backgrounds.
# For the main STREME analysis, I would use background_nonDARs.fa.

for (category in unique_categories) {
  
  bdg_data <- DARs_for_STREME %>%
    dplyr::filter(category_name != category) %>%
    dplyr::distinct(peakID, .keep_all = TRUE)
  
  bdg_data.gr <- GRanges(bdg_data$peakID)
  
  bdg_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, bdg_data.gr)
  names(bdg_sequences) <- bdg_data$peakID
  
  writeXStringSet(
    bdg_sequences,
    filepath = glue("{background_fasta_dir}/{category}_background_all_other_DARs.fa")
  )
}


###############################################
# 5. Save STREME input table and commands
###############################################

streme_inputs <- data.frame(
  category_name = unique_categories,
  positive_fasta = glue("{category_fasta_dir}/{unique_categories}.fa"),
  background_fasta = rep(
    glue("{background_fasta_dir}/background_nonDARs.fa"),
    length(unique_categories)
  ),
  suggested_command = glue(
    "streme --p {category_fasta_dir}/{unique_categories}.fa --n {background_fasta_dir}/background_nonDARs.fa --dna --oc streme_{unique_categories}"
  )
)

write.csv(
  streme_inputs,
  file = glue("{fasta_dir}/STREME_input_files.csv"),
  row.names = FALSE
)

writeLines(
  streme_inputs$suggested_command,
  con = glue("{fasta_dir}/suggested_STREME_commands.sh")
)


###############################################
# 6. GC content analysis
###############################################
# 6.1 Random non-DAR control sets for plotting

no_DARs.gr <- GRanges(background_noDARs$peakID)

sample_200 <- sample(no_DARs.gr, 200)
sample_800 <- sample(no_DARs.gr, 800)
sample_2000 <- sample(no_DARs.gr, 2000)

random_noDARs_sets <- GRangesList(sample_200, sample_800, sample_2000)
names(random_noDARs_sets) <- c("sample_200", "sample_800", "sample_2000")

saveRDS(
  random_noDARs_sets,
  glue("{table_dir}/random_noDARs_sets.rds")
)


# 6.2 Combine DARs and random non-DAR controls

DARs_forGC.gr <- GRanges(
  DARs_for_STREME$peakID,
  category_name = DARs_for_STREME$category_name
)

noDARs_forGC.gr <- unlist(random_noDARs_sets)
noDARs_forGC.gr$category_name <- names(noDARs_forGC.gr)
names(noDARs_forGC.gr) <- NULL

all_peaks_forGCtesting <- c(DARs_forGC.gr, noDARs_forGC.gr)

table(all_peaks_forGCtesting$category_name)


# 6.3 Extract sequences

all_peaks_forGCtesting_seqs <- getSeq(
  BSgenome.Mmusculus.UCSC.mm10,
  all_peaks_forGCtesting
)

names(all_peaks_forGCtesting_seqs) <- all_peaks_forGCtesting$category_name


# 6.4 Compute GC-related metrics
# Values are proportions from 0 to 1, not percentages.

mononucleotides <- data.frame(
  oligonucleotideFrequency(all_peaks_forGCtesting_seqs, 1),
  length = width(all_peaks_forGCtesting)
) %>%
  dplyr::mutate(GC = (C + G) / length) %>%
  dplyr::select(GC)

dinucleotides <- data.frame(
  oligonucleotideFrequency(all_peaks_forGCtesting_seqs, 2),
  length = width(all_peaks_forGCtesting)
) %>%
  dplyr::mutate(CGdi = CG / length) %>%
  dplyr::select(CGdi)

tetranucleotides <- data.frame(
  oligonucleotideFrequency(all_peaks_forGCtesting_seqs, 4),
  length = width(all_peaks_forGCtesting)
) %>%
  dplyr::mutate(CGCG = CGCG / length) %>%
  dplyr::select(CGCG)

CpH <- data.frame(
  oligonucleotideFrequency(all_peaks_forGCtesting_seqs, 2),
  length = width(all_peaks_forGCtesting)
) %>%
  dplyr::mutate(CpH = (CA + CT + CC) / length) %>%
  dplyr::select(CpH)

all_peaks_forGCtesting.df <- as.data.frame(all_peaks_forGCtesting)

category_order <- c("sample_200", "sample_800", "sample_2000",
                    "transient_down", "transient_up", "recovered_down", "recovered_up", "memory_down", "memory_up", "delayed_down", "delayed_up")
all_peaks_forGCtesting.GCinfo <- cbind(
  all_peaks_forGCtesting.df,
  mononucleotides,
  dinucleotides,
  tetranucleotides,
  CpH
) %>%
  dplyr::select(category_name, GC, CGdi, CGCG, CpH) %>%
  dplyr::mutate(
    category_name = factor(category_name, levels = category_order)
  )

write.csv(
  all_peaks_forGCtesting.GCinfo,
  file = glue("{table_dir}/GC_content_metrics_per_peak.csv"),
  row.names = FALSE
)

GC_summary <- all_peaks_forGCtesting.GCinfo %>%
  dplyr::group_by(category_name) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median_GC = median(GC, na.rm = TRUE),
    median_CGdi = median(CGdi, na.rm = TRUE),
    median_CGCG = median(CGCG, na.rm = TRUE),
    median_CpH = median(CpH, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  GC_summary,
  file = glue("{table_dir}/GC_content_metrics_summary.csv"),
  row.names = FALSE
)


###############################################
# 7. GC content plots
###############################################
# 7.1 GC content plot
plot_GC_metric <- function(metric, title, ylab, outfile) {
  p <- ggplot(all_peaks_forGCtesting.GCinfo, aes(x = category_name, y = .data[[metric]], fill = category_name)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.8, width = 0.5, colour = "black") +
    scale_fill_manual(values = category_colors) +
    theme_classic() + theme(legend.position = "none", plot.title = element_text(size = 16, face = "bold"), axis.title = element_text(size = 14), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks = element_line(linewidth = 1)) +
    labs(title = title, x = "", y = ylab) +
    stat_summary(aes(label = round(after_stat(y), 2)), fun = median, geom = "text", color = "black", size = 4, vjust = -0.6)
  print(p); ggsave(glue("{plot_dir}/{outfile}.pdf"), p, width = 8, height = 5); ggsave(glue("{plot_dir}/{outfile}.png"), p, width = 8, height = 5, dpi = 300); return(p)
}

GCplot   <- plot_GC_metric("GC",   "GC content per peak category",                  "GC content, G+C proportion",      "GC_content_per_category")
CGdiplot <- plot_GC_metric("CGdi", "CpG dinucleotide content per peak category",    "CpG dinucleotide proportion",     "CpG_dinucleotide_per_category")
CGCGplot <- plot_GC_metric("CGCG", "CGCG tetranucleotide content per peak category", "CGCG tetranucleotide proportion", "CGCG_tetranucleotide_per_category")
CpHplot  <- plot_GC_metric("CpH",  "CpH dinucleotide content per peak category",     "CpH dinucleotide proportion",     "CpH_dinucleotide_per_category")

GCplot/CGdiplot

# 7.1 GC content plot
# dev.size()
# width_original = 6.76
# height_original= 8.36
# 
# plot_name <- "chapter3_atac_CG_panels.pdf"
# 
# pdf(plot_name, width = width_original, height =height_original)
# GCplot/CGdiplot
# dev.off()



# ========== Save workspace objects ==========
# save.image(glue("{dir}/STREME_input_files_and_GCcontent.rds"))


