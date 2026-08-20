# ===========================================
# Script Title: 3.12_STREME_and_GCcontent.R
# Author: Luna Zea Redondo
# Date: 2024-12-04
# Description:
#   Performs discriminative motif discovery and GC-content analyses
#   for dynamic ATAC-seq peak classes.
#
#   Main steps:
#     1. Extract FASTA sequences per peak category
#     2. Generate background sets (all DARs / direction-specific)
#     3. Separate promoter vs non-promoter regions
#     4. Prepare random non-DAR control sets
#     5. Compute GC-related metrics
#     6. Generate comparative plots
#
#   Inputs:
#     - DARs_simp
#     - DAR_complete_results_kmeans
#     - VTA_RPKM_formatted_metadata_withDARs
#     - proj6.VTA (optional, for peak annotation)
#
#   Outputs:
#     - FASTA files (per category and region)
#     - GC metrics tables
#     - Plots
# ===========================================

#!/usr/bin/env Rscript

# ========== Environment setup ==========
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() 
#.rs.restartR()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))

# ========== Libraries ==========
library(dplyr)
library(tidyr)
library(glue)
library(ggplot2)
library(scales)

library(GenomicRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)

library(ArchR)   
library(ggpubr) 

# ========== Parameters & paths ==========
set.seed(1)

addArchRThreads(threads = 16)
addArchRGenome("mm10")

condition_colors <- c("black", "#9F2365", "#617641", "#C48208", "#326186", "#AE430A", "#564686")
condition_names <- c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")
names(condition_colors) <- condition_names

dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/7_241204_STREME_and_GCcontent"
setwd(dir)


###############################################
# 1. Discriminative TF enrichment (ATAC classes)
###############################################

# 1.1 Extract sequences per category and backgrounds
# Iterate through each unique category
for (category in unique_categories) {
  #a. Subset the data for the current category
  subset_data <- DARs_simp[DARs_simp$category_name == category, ]
  subset_data.gr <- GRanges(subset_data$peakID)
  
  #b. Get fasta sequences and save:
  input_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, subset_data.gr)
  names(input_sequences) <- subset_data$peakID
  writeXStringSet(input_sequences, filepath = glue("{region_dir}{category}_sequences.tsv"))
  
  #c. Generate two backgrounds:
  #c.1 All other DARs:
  bdg_data <- DARs_simp[DARs_simp$category_name != category, ]
  bdg_data.gr <- GRanges(bdg_data$peakID)
  bdg_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, bdg_data.gr)
  names(bdg_sequences) <- bdg_data$peakID
  writeXStringSet(bdg_sequences, filepath = glue("{category}_BDG_all_other_DARs.tsv"))
  
  #c.1 UP/DOWN DARs:
  specific_bdg <- DARs_simp[DARs_simp$category_name != category, ] %>% 
    tidyr::separate(category_name, into=c("category", "direction"), remove = FALSE) %>% 
    dplyr::filter(ifelse(category %in% c("transient_up", "recovered_up", "memory_up", "delayed_up"), direction == "up", direction == "down"), 
                  category_name != category)
  
  specific_bdg.gr <- GRanges(specific_bdg$peakID)
  specific_bdg_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, specific_bdg.gr)
  names(specific_bdg_sequences) <- specific_bdg$peakID
  writeXStringSet(specific_bdg_sequences, filepath = glue("{category}_BDG_specific_direction.tsv"))
}

# 1.2 Separate promoters vs other regions
region_dir <- glue("{dir}/per_region")

#Separate promoters / other regulatory 
all_peaks_region <- as.data.frame(getPeakSet(proj6.VTA)) %>% 
  dplyr::mutate(peakID = glue("{seqnames}:{start}-{end}")) %>% 
  dplyr::select(peakID, peakType)
DARs_simp_region <- left_join(DARs_simp, all_peaks_region, by = "peakID")

for (category in unique_categories) {
  # Subset the data for the current category
  subset_data <- DARs_simp_region[DARs_simp_region$category_name == category, ]
  
  # Separate Promoters and Other Peaks
  promoters <- subset_data[subset_data$peakType == "Promoter", ]
  others <- subset_data[subset_data$peakType != "Promoter", ]
  
  # Convert to GRanges
  promoters.gr <- GRanges(promoters$peakID)
  others.gr <- GRanges(others$peakID)
  
  # Get FASTA sequences for Promoters
  promoters_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, promoters.gr)
  names(promoters_sequences) <- promoters$peakID
  writeXStringSet(promoters_sequences, filepath = glue("{region_dir}/{category}_promoters_sequences.tsv"))
  
  # Get FASTA sequences for Other Peaks
  others_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, others.gr)
  names(others_sequences) <- others$peakID
  writeXStringSet(others_sequences, filepath = glue("{region_dir}/{category}_others_sequences.tsv"))
  
  # Backgrounds for Promoters and Others
  bdg_data <- DARs_simp_region[DARs_simp_region$category_name != category, ]
  
  # Separate Promoters and Other Peaks in Background
  bdg_promoters <- bdg_data[bdg_data$peakType == "Promoter", ]
  bdg_others <- bdg_data[bdg_data$peakType != "Promoter", ]
  
  # Convert Backgrounds to GRanges
  bdg_promoters.gr <- GRanges(bdg_promoters$peakID)
  bdg_others.gr <- GRanges(bdg_others$peakID)
  
  # Get Background Sequences for Promoters
  bdg_promoters_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, bdg_promoters.gr)
  names(bdg_promoters_sequences) <- bdg_promoters$peakID
  writeXStringSet(bdg_promoters_sequences, filepath = glue("{region_dir}/{category}_BDG_promoters.tsv"))
  
  # Get Background Sequences for Other Peaks
  bdg_others_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, bdg_others.gr)
  names(bdg_others_sequences) <- bdg_others$peakID
  writeXStringSet(bdg_others_sequences, filepath = glue("{region_dir}/{category}_BDG_others.tsv"))
  
  # Specific Background for UP/DOWN DARs
  specific_bdg <- bdg_data %>%
    tidyr::separate(category_name, into = c("category", "direction"), remove = FALSE) %>%
    dplyr::filter(ifelse(category %in% c("transient_up", "recovered_up", "memory_up", "delayed_up"), 
                         direction == "up", direction == "down"),
                  category_name != category)
  
  # Separate Promoters and Others in Specific Background
  specific_bdg_promoters <- specific_bdg[specific_bdg$peakType == "Promoter", ]
  specific_bdg_others <- specific_bdg[specific_bdg$peakType != "Promoter", ]
  
  # Convert to GRanges
  specific_bdg_promoters.gr <- GRanges(specific_bdg_promoters$peakID)
  specific_bdg_others.gr <- GRanges(specific_bdg_others$peakID)
  
  # Get Specific Background Sequences for Promoters
  specific_bdg_promoters_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, specific_bdg_promoters.gr)
  names(specific_bdg_promoters_sequences) <- specific_bdg_promoters$peakID
  writeXStringSet(specific_bdg_promoters_sequences, filepath = glue("{region_dir}/{category}_BDG_specific_promoters.tsv"))
  
  # Get Specific Background Sequences for Other Peaks
  specific_bdg_others_sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, specific_bdg_others.gr)
  names(specific_bdg_others_sequences) <- specific_bdg_others$peakID
  writeXStringSet(specific_bdg_others_sequences, filepath = glue("{region_dir}/{category}_BDG_specific_others.tsv"))
}



###############################################
# 3. GC content analysis
###############################################
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/8_241205_GCcontent"
setwd(dir)

# 2.1 Generate random non-DAR control sets

#A) Calculate median of number of peaks and length in all categories:
median_number_peaks <- DAR_complete_results_kmeans %>%
  dplyr::mutate(category_name = glue("{memory_class}_{direction}")) %>%
  dplyr::count(category_name) %>%
  dplyr::summarise(median_peaks = median(n)) %>%
  dplyr::pull(median_peaks)

median_length_peaks <- DAR_complete_results_kmeans %>%
  dplyr::mutate(
    start = as.numeric(sub(".*:(\\d+)-.*", "\\1", peakID)),
    end = as.numeric(sub(".*-(\\d+)", "\\1", peakID)),
    peak_length = end - start) %>% 
  dplyr::summarise(median_length = median(peak_length, na.rm = TRUE)) %>%
  pull(median_length)


#Sanity check: Calculate the median length per category
# median_length_per_category <- DARs_length %>%
#   dplyr::group_by(category_name) %>%
#   dplyr::summarise(median_length = median(length))


#B) Generate three random sets with the previous values for the median:
no_DARs <- VTA_RPKM_formatted_metadata_withDARs %>% 
  dplyr::filter(diff == "No", 
                peakID %notin% DAR_complete_results_kmeans$peakID) %>% 
  dplyr::select(peakID) %>% distinct() 

no_DARs.gr <- GRanges(no_DARs$peakID)

#sample: 200,800,200 peaks
sample_200 <- sample(no_DARs.gr, 200)
median(width(sample_200)) #406 bp

sample_800 <- sample(no_DARs.gr, 800)
median(width(sample_800)) #401 bp

sample_2000 <- sample(no_DARs.gr, 2000)
median(width(sample_2000)) #410 bp

random_noDARs_sets <- GRangesList(sample_200, sample_800, sample_2000)
names(random_noDARs_sets) <- c("sample_200", "sample_800", "sample_2000")
#saveRDS(random_noDARs_sets, glue("241205_random_noDARs_sets.rds"))

# 2.2 Prepare combined dataset for GC testing
DARs_forGC <- DAR_complete_results_kmeans %>% 
  dplyr::mutate(category_name = glue("{memory_class}_{direction}")) %>% 
  dplyr::select(peakID, category_name)

DARs_forGC.gr <- GRanges(DARs_forGC$peakID, category_name = DARs_forGC$category_name)
noDARs_forGC.gr <- unlist(random_noDARs_sets)
noDARs_forGC.gr$category_name = names(noDARs_forGC.gr)
names(noDARs_forGC.gr) <- NULL

all_peaks_forGCtesting <- c(DARs_forGC.gr, noDARs_forGC.gr)
#table(all_peaks_forGCtesting$category_name) #sanity check

# 2.3 Extract sequences
all_peaks_forGCtesting_seqs <- BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, all_peaks_forGCtesting)
names(all_peaks_forGCtesting_seqs) <- all_peaks_forGCtesting$category_name


# 2.4 Compute GC-related metrics
mononucleotides <- data.frame(oligonucleotideFrequency(all_peaks_forGCtesting_seqs, 1),
                              length = width(all_peaks_forGCtesting)) %>% 
  dplyr::mutate(GC = round(100 * (C + G)/length, 1)) %>% 
  dplyr::select(GC)

dinucleotides <- data.frame(oligonucleotideFrequency(all_peaks_forGCtesting_seqs, 2), 
                            length = width(all_peaks_forGCtesting)) %>% 
  dplyr::mutate(CGdi = round(100 * CG/length, 1)) %>% 
  dplyr::select(CGdi)

tetranucleotides <- data.frame(oligonucleotideFrequency(all_peaks_forGCtesting_seqs, 4),
                               length = width(all_peaks_forGCtesting)) %>% 
  dplyr::mutate(CGCG = round(100 * CGCG/length, 1)) %>% 
  dplyr::select(CGCG)

CpH <- data.frame(oligonucleotideFrequency(all_peaks_forGCtesting_seqs, 2), 
                  length = width(all_peaks_forGCtesting)) %>% 
  dplyr::mutate(CpH = round(100 * (CA+CT+CC)/length, 1)) %>% 
  dplyr::select(CpH)

all_peaks_forGCtesting.df <- as.data.frame(all_peaks_forGCtesting)
all_peaks_forGCtesting.GCinfo <- cbind(all_peaks_forGCtesting.df, mononucleotides, dinucleotides, tetranucleotides, CpH) %>% 
  dplyr::select(category_name, GC, CGdi, CGCG, CpH)

# 2.5 Plot GC metrics
GCplot <- ggplot(all_peaks_forGCtesting.GCinfo, aes(x = category_name, y = GC, fill = category_name, colour = category_name)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust =2, trim = TRUE) +
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(aes(x = as.integer(category_name) + 0.25, y = GC), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_fill_manual(values = memory_classes_colors) +
  scale_color_manual(values = memory_classes_colors) +
  theme_classic() + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1), # Incline x-axis labels,
        axis.ticks = element_line(size = 1)) + 
  labs(title = "GC content per peak category", x= "", y= "%GC content (G+C)") +
  stat_summary(aes(label = round(..y.., 1)), fun = median, geom = "text", color = "black", size = 4, position = position_nudge(x = 0.5))
GCplot

CGdiplot <- ggplot(all_peaks_forGCtesting.GCinfo, aes(x = category_name, y = CGdi, fill = category_name, colour = category_name)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust =2, trim = TRUE) +
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(aes(x = as.integer(category_name) + 0.25, y = CGdi), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_fill_manual(values = memory_classes_colors) +
  scale_color_manual(values = memory_classes_colors) +
  theme_classic() + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1), # Incline x-axis labels,
        axis.ticks = element_line(size = 1)) + 
  labs(title = "CpG (dinucleotide) content per peak category", x= "", y= "%CpG dinucleotide") +
  stat_summary(aes(label = round(..y.., 1)), fun = median, geom = "text", color = "black", size = 4, position = position_nudge(x = 0.4))
CGdiplot

CGCGplot <- ggplot(all_peaks_forGCtesting.GCinfo, aes(x = category_name, y = CGCG, fill = category_name, colour = category_name)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust =2, trim = TRUE) +
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(aes(x = as.integer(category_name) + 0.25, y = CGCG), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_fill_manual(values = memory_classes_colors) +
  scale_color_manual(values = memory_classes_colors) +
  theme_classic() + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1), # Incline x-axis labels,
        axis.ticks = element_line(size = 1)) + 
  labs(title = "CGCG (tetranucleotide) content per peak category", x= "", y= "%CGCG tetranucleotide") +
  stat_summary(aes(label = round(..y.., 1)), fun = median, geom = "text", color = "black", size = 4, position = position_nudge(x = 0.4))
CGCGplot

CpHplot <- ggplot(all_peaks_forGCtesting.GCinfo, aes(x = category_name, y = CpH, fill = category_name, colour = category_name)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust =2, trim = TRUE) +
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(aes(x = as.integer(category_name) + 0.25, y = CpH), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_fill_manual(values = memory_classes_colors) +
  scale_color_manual(values = memory_classes_colors) +
  theme_classic() + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1), # Incline x-axis labels,
        axis.ticks = element_line(size = 1)) + 
  labs(title = "CpH dinucleotide content per peak category", x= "", y= "%[CpH] dinucleotide") +
  stat_summary(aes(label = round(..y.., 1)), fun = median, geom = "text", color = "black", size = 4, position = position_nudge(x = 0.4))
CpHplot

CGdiplot/CpHplot
#load("241205_CGcontent.rds")

#Compact plot (24/12/11:)
all_peaks_forGCtesting.GCinfo.v2 <- all_peaks_forGCtesting.GCinfo %>% dplyr::select(-CGCG)

# Prepare the data
medians_long <- all_peaks_forGCtesting.GCinfo %>%
  dplyr::select(-CGCG) %>%
  dplyr::group_by(category_name) %>%
  dplyr::summarise(across(c(GC, CGdi, CpH), ~ median(.x, na.rm = TRUE), .names = "median_{.col}")) %>%
  pivot_longer(cols = starts_with("median_"), names_to = "metric", values_to = "median_value") %>%
  dplyr::mutate(metric = gsub("median_", "", metric)) %>% 
  dplyr::mutate(metric = factor(metric, levels = c("GC", "CGdi", "CpH"))) 

# Plot
ggplot(medians_long, aes(x = category_name, y = median_value)) +
  geom_point(aes(color = category_name, shape = metric), size = 4) +
  geom_text(aes(label = sprintf("%.1f", median_value), y = median_value + 0.5), vjust = -0.5, show.legend = FALSE) +
  scale_color_manual(values = memory_classes_colors) +
  scale_shape_manual(values = c("GC" = 16, "CGdi" = 17, "CpH" = 18)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",         # Position legends at the bottom
    legend.box = "horizontal"           # Arrange legends horizontally
  ) +
  labs(
    title = "CG content metrics",
    x = NULL,
    y = "Median %",
    color = "Peak dyanmic class",
    shape = "Metric"
  )
#save.image("241211.GCcontent.finalPlot.rds")

#Extra: compact figure
# Define factor levels in desired order
all_peaks_forGCtesting.GCinfo$category_name <- factor(
  all_peaks_forGCtesting.GCinfo$category_name,
  levels = c("sample_200", "sample_800", "sample_2000",
             "transient_down", "transient_up", "recovered_down", "recovered_up",
             "memory_down", "memory_up", "delayed_down", "delayed_up")
)

# Define your custom colors
memory_classes_colors <- c(
  "transient_up"="#D9481C", "recovered_up"="#F3A712", "memory_up" ="#008AB8", "delayed_up" = "#954B77",
  "transient_down"="#EF9A81", "recovered_down"="#F8CD77", "memory_down" ="#8FD0E6", "delayed_down" = "#C892B3",
  "sample_200"= "gray70", "sample_800"= "gray70", "sample_2000"= "gray70"
)

# Melt to long format
long_df <- melt(all_peaks_forGCtesting.GCinfo, id.vars = "category_name")

# Calculate group means for annotation
mean_labels <- long_df %>%
  group_by(category_name, variable) %>%
  summarise(mean_value = mean(value), .groups = 'drop')

# Plot with annotations
GCmetrics <- ggplot(long_df, aes(x = category_name, y = value, fill = category_name)) +
  geom_boxplot() +
  geom_text(data = mean_labels, aes(x = category_name, y = mean_value + 0.05 * max(long_df$value), 
                                    label = round(mean_value, 1)), 
            inherit.aes = FALSE, size = 3, angle = 90, vjust = -0.2) +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = memory_classes_colors) +
  theme_minimal() +
  labs(x = "Category", y = "Value", title = "GC-related Metrics by Category") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )



