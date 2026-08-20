# ===========================================
# Script Title: ATAC exploratory analysis (VTA DNs)
# Author: Luna Zea Redondo
# Date: 2024-03-18
#
# Description:
#   This script performs a exploratory analysis on VTA-DN peakset,
#   integrating chromatin accessibility metrics, differential accessibility results (DARs),
#   motif enrichment, and transcription factor binding information.
#
#
#   The workflow includes:
#     - Integration of ArchR peak metadata with RPKM-normalized accessibility values
#     - Incorporation of differential accessibility statistics across cocaine timepoints
#     - Generation of MA plots for DAR visualization
#     - Motif enrichment analysis using Signac
#     - Identification of candidate transcription factors and binding sites
#     - Construction of TF–DAR matrices and upset plots
#     - Entropy-based accessibility analyses across RNA-defined clusters
# ===========================================


# ========== 0. Environment Setup ==========

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() 
#.rs.restartR()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))

# ========== 1. Load Required Libraries ==========

library(Seurat)
library(Signac)
library(ArchR)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(tidyr)
library(readr)
library(glue)
library(ggplot2)
library(ggside)
library(scales)
library(ComplexHeatmap)
library(UpSetR)
library(limma)
library(dendextend)
library(dendsort)
library(gtools)
library(patchwork)


# ========== 2. ArchR Configuration ==========
addArchRThreads(threads = 16)
addArchRGenome("mm10")

# ========== 3. Define Paths and Parameters ==========
#source("/fast/AG_Pombo/luna/2023_vta_multiome/6_231107_updated_analyses/3_ATAC_master_table/GTFgetAnnotation.R")
#DNs.RNA.seu.VTA <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/6_231107_updated_analyses/3_ATAC_master_table/240123.DNs.RNA.seu.VTA.rds")
genes_tracked <- c("Vip", "Cartpt", "Nr4a1", "Rbfox1", "Pcdh9", "Bdnf")

dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal"
setwd(dir)

all_contrasts <- c("h1_cocaine-saline",
                   "h4_cocaine-saline",
                   "h8_cocaine-saline",
                   "h24_cocaine-saline",
                   "d14_cocaine-saline", 
                   "h4_cocaine-h1_cocaine",
                   "h8_cocaine-h4_cocaine", 
                   "h24_cocaine-h8_cocaine", 
                   "d14_cocaine-h24_cocaine")


# ========== 4. Load Input Data ==========

# ArchR and seurat objects
proj6.VTA <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/6_231107_updated_analyses/2_ATAC_alternative_processing/2_231123_condition_specific_peaks/231218.proj6.VTA.annotated.rds")
wnn.seu.VTA <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/6_231107_updated_analyses/2_ATAC_alternative_processing/3_231205_ArchR_DARsandMotifs/231219_rpkm_based_DARs/240119_thresholds_setC/240122.wnn.seu.VTA.withMotifs.rds")

wnn.seu.VTA.1291 <- subset(wnn.seu.VTA, subset = orig.ident %notin% c("m30_cocaine_R1", "h1_saline_R1"))

# RPKMs work (formatted for new DARs):
VTA_DARs_RPKM_complete <-  read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/6_231107_updated_analyses/2_ATAC_alternative_processing/3_231205_ArchR_DARsandMotifs/231219_rpkm_based_DARs/240119_thresholds_setC/240122.VTA_DARs_RPKM_realClass_TFmatches_complete.tsv")
VTA_RPKM_formatted <- VTA_DARs_RPKM_complete[, c(5,1,3,4,22,23,21,17,19,20,18,16,29,25,27,28,26,24,30:44,6,45:49)]

# Peakset annotation:
peaks_metadata <- data.frame(getPeakSet(proj6.VTA)) %>% dplyr::mutate(peakID = glue("{seqnames}:{start}-{end}")) %>% 
  dplyr::select(width, distToGeneStart, nearestGene, peakType, distToTSS, nearestTSS, GC, peakID)

VTA_RPKM_formatted_metadata <- left_join(VTA_RPKM_formatted, peaks_metadata, by="peakID")

#Set parameters:
pval = 0.05
logfc = 0.25
logcpm = 3

# DARs from Vedran: UPDATED LIST
DARs_corrected <- read_tsv("241113.peaks-edgeR.tsv") %>% 
  dplyr::mutate(contrast = ifelse(contrast == "h4_cocaine-h8_cocaine", "h8_cocaine-h4_cocaine",
                                  ifelse(contrast == "h1_cocaine-h4_cocaine", "h4_cocaine-h1_cocaine", contrast)),
                logFC = ifelse(contrast %in% c("h8_cocaine-h4_cocaine", "h4_cocaine-h1_cocaine"), -logFC, logFC), 
                unshrunk.logFC = ifelse(contrast %in% c("h8_cocaine-h4_cocaine", "h4_cocaine-h1_cocaine"), -unshrunk.logFC, unshrunk.logFC)) %>% 
  dplyr::filter(contrast %in% all_contrasts) %>% 
  dplyr::mutate(peakID = str_replace(gene_id, "(\\w+)-(\\d+)-(\\d+)", "\\1:\\2-\\3"), 
                diff = ifelse(abs(logFC) >= logfc & p_val < pval & logCPM >= logcpm, "Yes", "No")) %>% 
  dplyr::rename(peakID_v2 = "gene_id")

# Final "basic" table"
VTA_RPKM_formatted_metadata_withDARs <- left_join(VTA_RPKM_formatted_metadata, DARs_corrected, by = c("peakID", "contrast"))
#saveRDS(VTA_RPKM_formatted_metadata_withDARs, "/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/241121_VTA_RPKM_formatted_metadata_withDARs.rds")


# ========== 5. MA Plot Visualization ==========
dirplots <- dir #In this case. 
condition_colors <- c("black", "#9F2365", "#617641", "#C48208", "#326186", "#AE430A", "#564686")
condition_names <- c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")
names(condition_colors) <- condition_names

MAplots <- list()
for (contrast in all_contrasts) {
  message(glue("Analyzing: {contrast}"))
  group1=str_split(contrast, "-")[[1]][1]
  group2=str_split(contrast, "-")[[1]][2]
  contrast_label=contrast
  
  volcano_colors <- c(condition_colors[group1], condition_colors[group2], "gray80")
  names(volcano_colors) <- c("Upregulated", "Downregulated", "No significant")
  
  DARS_to_plot <- DARs_corrected %>% #load
    dplyr::filter(contrast == contrast_label) %>% 
    dplyr::mutate(is_DAR = ifelse(p_val < pval & logFC >= logfc & logCPM >= logcpm, "Upregulated", 
                                  ifelse(p_val < pval & logFC <= -logfc & logCPM >= logcpm, "Downregulated", "No significant")))
  
  downreg_num <- table(DARS_to_plot$is_DAR)["Downregulated"] 
  upreg_num <- table(DARS_to_plot$is_DAR)["Upregulated"] 
  
  MA_plot <- ggplot(DARS_to_plot, aes(x = logCPM, y = logFC)) +
    geom_point(alpha=0.6, aes(color = is_DAR)) +
    geom_point(alpha=0.6,aes(color = is_DAR), data=DARS_to_plot[DARS_to_plot$is_DAR != "No significant",]) +
    scale_color_manual(values = volcano_colors) +
    theme_classic(base_size = 12) + theme(legend.position = "bottom") +
    labs(title = glue("{contrast_label}"),
         subtitle = glue("{upreg_num} up- and {downreg_num} down-regulated"), 
         x= "logCPM", 
         y="LogFC") +
    guides(fill=guide_legend(title="is DARs")) +
    geom_vline(xintercept = logcpm, linetype="dotted", color = "darkred", linewidth=1.5) + 
    theme(legend.position = "none")
  MA_plot
  #ggsave(filename = glue("{dirplots}/{contrast_label}_peakANDmotifs_results_MAplot.png"), MA_plot, units = "px", device = "png", width=6000,height=2000,dpi = 300)  
  MAplots[[paste("MAplot", contrast, sep = "_")]] <- MA_plot
}

MA_row1 <- grid.arrange(grobs = MAplots[1:5], nrow = 1)
MA_row2 <- grid.arrange(grobs = MAplots[6:9], nrow = 1)

#load("241113.ATAC_master_table.rds")
blank <- ggplot() + theme_void() 

# Arrange the plots
combined_layout <- (MAplots[[1]] + MAplots[[2]] + MAplots[[3]] + MAplots[[4]] + MAplots[[5]] +
                      MAplots[[6]] + MAplots[[7]] + MAplots[[8]] + MAplots[[9]] + blank) + plot_layout(ncol = 5)


# ========== 6. Motif Enrichment Analysis ==========
#All backgrounds: In this iteration, fitlered peaks  == all peaks
bdg_allPeaks <- VTA_RPKM_formatted_metadata_withDARs %>% dplyr::select(peakID) %>%  mutate(peakID = str_replace(peakID, ":", "-")) %>% distinct() %>% pull()
bdg_filteredPeaks <- VTA_RPKM_formatted_metadata_withDARs %>% dplyr::filter(complete.cases(.)) %>% dplyr::select(peakID_v2) %>% distinct() %>% pull()
bdg_allDARs <- VTA_RPKM_formatted_metadata_withDARs %>% dplyr::filter(complete.cases(.), p_val < pval, abs(logFC) >=logfc, logCPM >= logcpm) %>% dplyr::select(peakID_v2) %>% distinct() %>% pull()

#Save DARs in table (for kmeans):
DARs.table <- VTA_RPKM_formatted_metadata_withDARs %>% dplyr::filter(complete.cases(.), p_val < pval, abs(logFC) >=logfc, logCPM >= logcpm) %>%
  dplyr::distinct(peakID)
#write_tsv(DARs.table, "241113_vta.edgeR.DARs.txt")

# Add motif information to seurat object
pfm <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/1_ATAC_masterTable_continued/231217.pfm.rds")
wnn.seu.VTA.1291 <- AddMotifs(
  object = wnn.seu.VTA.1291,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
TFenrichment.complete.df <- data.frame()

#Extract enriched TFs
for (contrast_label in all_contrasts) {
  message(glue("Analyzing: {contrast_label}"))
  group1=str_split(contrast_label, "-")[[1]][1]
  group2=str_split(contrast_label, "-")[[1]][2]
  print(contrast_label)
  
  #A. Extract DARs and divide by contrast and direction
  contrast.group1.dars <- DARs_corrected %>% dplyr::filter(contrast == contrast_label, logFC >= logfc, p_val < pval, logCPM >= logcpm) %>% dplyr::select(peakID_v2) %>% pull()
  contrast.group2.dars <- DARs_corrected %>%dplyr::filter(contrast == contrast_label, logFC <= -logfc, p_val < pval, logCPM >= logcpm) %>% dplyr::select(peakID_v2) %>% pull()
  
  
  #B.1. TFenrichment; bdg: allPeaks
  contrast.group1vsallPeaks.tfenrich <- FindMotifs(
    object = wnn.seu.VTA.1291,
    features = contrast.group1.dars, 
    background = bdg_allPeaks) %>% 
    dplyr::mutate(contrast = contrast_label, upreg_in = group1, bdg = "allPeaks")
  
  contrast.group2vsallPeaks.tfenrich <- FindMotifs(
    object = wnn.seu.VTA.1291,
    features = contrast.group2.dars, 
    background = bdg_allPeaks) %>% 
    dplyr::mutate(contrast = contrast_label, upreg_in = group2, bdg = "allPeaks")
  
  #B.2. TFenrichment; bdg: filteredPeaks
  contrast.group1vsfilteredPeaks.tfenrich <- FindMotifs(
    object = wnn.seu.VTA.1291,
    features = contrast.group1.dars, 
    background = bdg_filteredPeaks) %>% 
    dplyr::mutate(contrast = contrast_label, upreg_in = group1, bdg = "filteredPeaks")
  
  contrast.group2vsfilteredPeaks.tfenrich <- FindMotifs(
    object = wnn.seu.VTA.1291,
    features = contrast.group2.dars, 
    background = bdg_filteredPeaks) %>% 
    dplyr::mutate(contrast = contrast_label, upreg_in = group2, bdg = "filteredPeaks")
  
  #B.3. TFenrichment; bdg: allDARs
  contrast.group1vsallDARs.tfenrich <- FindMotifs(
    object = wnn.seu.VTA.1291,
    features = contrast.group1.dars, 
    background = bdg_allDARs) %>% 
    dplyr::mutate(contrast = contrast_label, upreg_in = group1, bdg = "allDARs")
  
  contrast.group2vsallDARs.tfenrich <- FindMotifs(
    object = wnn.seu.VTA.1291,
    features = contrast.group2.dars, 
    background = bdg_allDARs) %>% 
    dplyr::mutate(contrast = contrast_label, upreg_in = group2, bdg = "allDARs")
  
  TFenrichment.complete.df <-rbind(TFenrichment.complete.df,
                                   contrast.group1vsallPeaks.tfenrich, contrast.group2vsallPeaks.tfenrich,
                                   contrast.group1vsfilteredPeaks.tfenrich, contrast.group2vsfilteredPeaks.tfenrich, 
                                   contrast.group1vsallDARs.tfenrich, contrast.group2vsallDARs.tfenrich)
}
#write_tsv(TFenrichment.complete.df, glue("{dirplots}/241113.TFenrichment.complete.df.tsv"))



# ========== 7. TF–DAR Matrix Construction ==========
#Candidate TFs to study (only comparisons vs saline):
preselected_TFs.df <- TFenrichment.complete.df %>%
  group_by(contrast, upreg_in, bdg) %>% 
  dplyr::select(pvalue, motif.name, contrast, upreg_in, bdg) %>% 
  dplyr::filter(pvalue < 0.01, bdg == "allDARs") %>% 
  dplyr::filter(contrast %in% all_contrasts[1:5]) 

preselected_TFs_unique_vec <- selected_TFs.df %>%
  dplyr::ungroup() %>% 
  dplyr::select(motif.name) %>% distinct() %>% pull()
#saveRDS(preselected_TFs_unique_vec, "241113_preselected_TFs_unique_vec.rds")

#Extract motif matrix and convert motifID to motif names. 
motif_motifName_conversion <- TFenrichment.complete.df %>% dplyr::select(motif, motif.name)
motif_matrix <- as.matrix(wnn.seu.VTA@assays[["peaks"]]@motifs@data)

motif_names <- motif_motifName_conversion$motif.name
colnames(motif_matrix) <- motif_names[match(colnames(motif_matrix), motif_motifName_conversion$motif)]

#Reduce matrix to candidate TFs:
motif_matrix_reduced <- motif_matrix[, colnames(motif_matrix) %in% selected_TFs_unique_vec] # 104545 X 184
motif_matrix_DARs <- motif_matrix_reduced[bdg_allDARs ,] #8201 X 184


# Important: save the DARs_corrected file, will come back to it often. 
# save.image("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/241113.ATAC_master_table.rds")
# saveRDS(DARs_corrected, "/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/241113.DARs_corrected.rds")



############################################################
## PART 2 — ATAC peak composition & quality control
############################################################

## ---- Setup / assumptions ----
# - proj6.VTA exists
# - required libraries already loaded (dplyr, ggplot2, tidyr, glue, patchwork)
# - VTA_RPKM_formatted_metadata_withDARs already in memory
# - contrast_order already defined


# ========== 7. Peak type distributions ==========
all_peaks <- as.data.frame(getPeakSet(proj6.VTA))

summary_data <- all_peaks %>%
  group_by(peakType) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  mutate(
    percentage = round((count / sum(count)) * 100, 1),
    label = paste0(percentage, "% (n = ", count, ")"),
    peakType = factor(peakType, levels = c("Distal", "Intronic", "Exonic", "Promoter"))
  )

# Define colors
peak_type_colors <- c(Promoter = "#8B0000", Exonic   = "#DAA520", Intronic = "#2E8B57",Distal   = "#4682B4")

# Bar chart
bar_chart <- ggplot(summary_data, aes(x = count, y = peakType, fill = peakType)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = peak_type_colors) +
  theme_classic() +
  labs(title = "Peak Type Distribution", x = "Count", y = "Peak Type") +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

# save
# dev.size()
# width_original = 4.479167
# height_original= 4.281250
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter3_peaktype_barchart"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# bar_chart
# dev.off()


# ========== 8. Simplification of metadata & DAR class ==========

VTA_RPKM_formatted_metadata_withDARs <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/241121_VTA_RPKM_formatted_metadata_withDARs.rds")
VTA_RPKM_formatted_metadata_withDARs_simplified <- VTA_RPKM_formatted_metadata_withDARs %>% 
  dplyr::select(peakID, contrast, length, bdgGroups_rpkm_value,  useGroups_rpkm_value, GC, diff) %>% 
  dplyr::mutate(direction = ifelse(diff == "Yes" & useGroups_rpkm_value > bdgGroups_rpkm_value, "upreg", 
                                   ifelse(diff == "Yes" & useGroups_rpkm_value < bdgGroups_rpkm_value, "downreg", "no_dar")))


# ========== 9. Evaluation of GC content ==========
gc_per_peakType <- VTA_RPKM_formatted_metadata_withDARs %>% 
  dplyr::select(peakID, GC, peakType) %>% 
  dplyr::distinct(peakID, .keep_all = TRUE) %>% 
  dplyr::mutate(peakType = factor(peakType, levels = c("Distal", "Intronic", "Exonic", "Promoter")))

gc_per_peakType.plot <- ggplot(gc_per_peakType, aes(x = GC, y = peakType, fill = peakType, colour = peakType)) +
  #  geom_flat_violin(position = position_nudge(x = 0.25, y = 0), adjust = 2, trim = TRUE) +  # Adjusted for horizontal layout
  geom_point(position = position_jitter(height = 0.15), size = 0.25) +  # Jitter points horizontally
  geom_boxplot(aes(x = GC, y = as.numeric(peakType) + 0.25), 
               outlier.shape = NA, alpha = 0.3, width = 0.1, colour = "BLACK") +  # Adjusted boxplot alignment
  scale_x_log10() +  # Log scale for GC content (now on X-axis)
  labs(title = "", x = "GC Content", y = "Peak Type") +
  scale_fill_manual(values = peak_type_colors) + 
  scale_color_manual(values = peak_type_colors) +
  theme_classic() +
  theme(legend.position = "none") 

# Print the plot
bar_chart + gc_per_peakType.plot

# dev.size()
# width_original = 8.968750
# height_original= 4.666667
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter3_allPeaks_number_and_GCcontent"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# bar_chart + gc_per_peakType.plot
# dev.off()

# ========== 10. ATAC-level QC big panel ==========

# Big panel of quality control
# Ensure 'contrast' and 'direction' have the correct factor levels for ordering
VTA_RPKM_formatted_metadata_withDARs_simplified <- VTA_RPKM_formatted_metadata_withDARs_simplified %>%
  mutate(contrast = factor(contrast, levels = all_contrasts),
         direction = factor(direction, levels = c("no_dar", "downreg", "upreg")))

# Define new color palette
direction_colors <- c("upreg" = "#E14C4B", "downreg" = "#4B8EBD", "no_dar" = "gray80")

# Compute median values per contrast and direction
compute_medians <- function(data, variable_name) {
  data %>%
    filter(variable == variable_name) %>%
    dplyr::group_by(contrast, direction) %>%
    dplyr::summarise(median_value = round(median(value, na.rm = TRUE), 3), .groups = "drop")
}

# Reshape data to long format
qc_data_long <- VTA_RPKM_formatted_metadata_withDARs_simplified %>%
  pivot_longer(cols = c("GC", "length", "useGroups_rpkm_value"), 
               names_to = "variable", 
               values_to = "value")

# --- 1. GC CONTENT PANEL ---
gc_median_values <- compute_medians(qc_data_long, "GC")

gc_panel <- ggplot(filter(qc_data_long, variable == "GC"), aes(x = direction, y = value, fill = direction)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  geom_text(data = gc_median_values, aes(x = direction, y = median_value, label = median_value),
            size = 3, color = "black", vjust = -0.5) +
  facet_wrap(~contrast, nrow = 1) +
  scale_fill_manual(values = direction_colors) +
  theme_classic() +
  labs(title = "GC Content Distribution", x = "DAR Category", y = "GC Content", fill = "Direction") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 2. RPKM PANEL ---
rpkm_median_values <- compute_medians(qc_data_long, "useGroups_rpkm_value")

rpkm_panel <- ggplot(filter(qc_data_long, variable == "useGroups_rpkm_value"), aes(x = direction, y = value, fill = direction)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  geom_text(data = rpkm_median_values, aes(x = direction, y = median_value, label = median_value),
            size = 3, color = "black", vjust = -0.5) +
  facet_wrap(~contrast, nrow = 1) +
  scale_fill_manual(values = direction_colors) +
  theme_classic() +
  labs(title = "RPKM Value Distribution", x = "DAR Category", y = "RPKM Value", fill = "Direction") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 3. PEAK LENGTH PANEL ---
length_median_values <- compute_medians(qc_data_long, "length")

length_panel <- ggplot(filter(qc_data_long, variable == "length"), aes(x = direction, y = value, fill = direction)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  geom_text(data = length_median_values, aes(x = direction, y = median_value, label = median_value),
            size = 3, color = "black", vjust = -0.5) +
  facet_wrap(~contrast, nrow = 1) +
  scale_fill_manual(values = direction_colors) +
  theme_classic() +
  labs(title = "Peak Length Distribution", x = "DAR Category", y = "Peak Length", fill = "Direction") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print all plots
print(gc_panel)
print(rpkm_panel)
print(length_panel)


# To be run for each panel
# dev.size()
# width_original = 11.854167
# height_original= 4.666667
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter3_ATACqc_length_panel"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# length_panel
# dev.off()