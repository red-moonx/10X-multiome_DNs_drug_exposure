# ===========================================
# Script Title: 3.3 VTA-DN ATAC Core Analysis
# Author: Luna Zea Redondo
# Date: 2026-05-27
# Description:
#    This script processes dopaminergic neuron (DN) ATAC-seq data using ArchR 
#    and Signac. It performs quality control and filtering of cells based on 
#    fragment counts and TSS enrichment, extracts sample-specific fragments, 
#    and runs MACS2 peak calling (both pooled and per condition) to build a 
#    final non-overlapping peak universe. It then constructs a Signac chromatin 
#    assay object, annotates peaks relative to genomic features (promoters, 
#    exons, introns, distal regions) using gene annotations, computes GC content 
#    metrics, and visualizes peak type distributions and GC content summaries.
# ===========================================


# ========== Set Environment ==========
rm(list = ls())
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230310/site-library"))
set.seed(1)

# ========== Load Required Libraries ==========
library(ArchR)
library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(glue)
library(ComplexHeatmap)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(patchwork)
library(readr)
library(forcats)
library(GenomeInfoDb)
library(ggpubr)
library(PupillometryR)
library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library(Matrix)

source("/fast/AG_Pombo/luna/2026_rebuttal/config.R", local = FALSE)

# ========== Set Working Directory ==========
dir <- "/fast/AG_Pombo/luna/2026_rebuttal/9_ATAC-QCandIntegration"
setwd(dir)

# ========== ArchR Setup ==========
addArchRThreads(threads = 16)
addArchRGenome("mm10")


# ========== 0) Load Data ==========

# Load ArchR project with all cells
proj_raw <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/8_ATAC-preprocessing/250520_VTA-DNs_proj_raw.rds")
seu.VTA_DNs <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/2_DN-evaluation/250507_seu.VTA_DNs.rds")


# ========== 1) ATAC QC — All Cells ==========
# QC scatterplots by sample and full dataset
# log10(nFrags) vs TSS enrichment
proj <- proj_raw
df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))

atac_qc_all_plot <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
#ggsave(filename = "atac_qc_all.png", atac_qc_all_plot, units = "px", device = "png", width=6000,height=2000,dpi = 300)

df2 <- as.data.frame(df) %>% rownames_to_column("cellNames") %>% 
  separate(cellNames, c("Sample", "CB"), sep = "#") %>% 
  dplyr::select(Sample, log10.nFrags., TSSEnrichment) %>% 
  dplyr::rename(log10nFrags = log10.nFrags.)


sample_ids <- c(
  "m30_cocaine_R1", "h1_saline_R1", "h1_cocaine_R1", "h1_cocaine_R2",
  "h4_saline_R1", "h4_cocaine_R1", "h4_cocaine_R2",
  "h8_saline_R1", "h8_cocaine_R1", "h8_cocaine_R2",
  "h24_saline_R1", "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3",
  "d14_saline_R1", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3"
)

allplots <- list()
for (sample in sample_ids) {
  df3 <- df2 %>% dplyr::filter(Sample == sample)
  name <- sample
  plot <- ggPoint(
    x = df3$log10nFrags,
    y = df3$TSSEnrichment,
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    title = sample,
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df3$log10nFrags, probs = 0.99)),
    ylim = c(0, quantile(df3$TSSEnrichment, probs = 0.99))) +
    geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed") +
    theme(plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0)) # Left margin 
  assign(name, plot)  
  allplots <- c(allplots, list(plot))
}

atac_qc_all <- ggarrange(plotlist = allplots, ncol = 6, nrow = 3)
#ggsave(filename = "atac_qc_all_persample.png", atac_qc_all, units = "px", device = "png", width=6000,height=2000,dpi = 300)



# ========== 2) Subset ArchR project to DN cells only ==========
# Cell names match format already, use directly
dn_cell_names <- colnames(seu.VTA_DNs)

# Check overlap
valid_cells <- dn_cell_names[dn_cell_names %in% getCellNames(proj_raw)]
cat("DN cells in Seurat:        ", length(dn_cell_names), "\n")
cat("DN cells found in ArchR:   ", length(valid_cells), "\n")
cat("Lost (not in ArchR):       ", length(dn_cell_names) - length(valid_cells), "\n")

# Subset DNs from proj_raw into proj_DNs
proj_DNs <- proj_raw[valid_cells, ] 
df <- getCellColData(proj_DNs, select = c("log10(nFrags)", "TSSEnrichment"))

atac_qc_all_plot <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
#ggsave(filename = "atac_qc_onlyDNs.png", atac_qc_all_plot, units = "px", device = "png", width=6000,height=2000,dpi = 300)

#Save
# width_original = 6.760417
# height_original= 5.083333
# 
# plot_name <- "chapter3_DNs_atac_allQC.pdf"
# 
# pdf(plot_name, width = width_original, height =height_original)
# atac_qc_all_plot
# dev.off()

df2 <- as.data.frame(df) %>% rownames_to_column("cellNames") %>% 
  separate(cellNames, c("Sample", "CB"), sep = "#") %>% 
  dplyr::select(Sample, log10.nFrags., TSSEnrichment) %>% 
  dplyr::rename(log10nFrags = log10.nFrags.)

allplots <- list()
for (sample in sample_ids) {
  df3 <- df2 %>% dplyr::filter(Sample == sample)
  name <- sample
  plot <- ggPoint(
    x = df3$log10nFrags,
    y = df3$TSSEnrichment,
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    title = sample,
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df3$log10nFrags, probs = 0.99)),
    ylim = c(0, quantile(df3$TSSEnrichment, probs = 0.99))) +
    geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed") +
    theme(plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0)) # Left margin 
  assign(name, plot)  
  allplots <- c(allplots, list(plot))
}

atac_qc_all <- ggarrange(plotlist = allplots, ncol = 6, nrow = 3)
#ggsave(filename = "atac_qc_all_persample_onlyDNs.png", atac_qc_all, units = "px", device = "png", width=6000,height=2000,dpi = 300)

#Another plot, more "paper-friendly"
# ========== QC Bar Plot — DN cells per sample ==========
df_summary <- df2 %>%
  dplyr::mutate(QC_pass = ifelse(TSSEnrichment >= 4 & log10nFrags >= log10(1000), "Pass", "Fail")) %>%
  group_by(Sample, QC_pass) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(total = sum(n), pct = n / total * 100) %>%
  ungroup() %>%
  mutate(Sample = factor(Sample, levels = sample_ids))

barPlot_atac_passing <- ggplot(df_summary, aes(x = Sample, y = pct, fill = QC_pass)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = distinct(df_summary, Sample, total),
            aes(x = Sample, y = 102, label = total), inherit.aes = FALSE, size = 3) +
  scale_fill_manual(values = c("Pass" = "grey80", "Fail" = "#c0392b"), name = "QC") +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25)) +
  labs(x = "Sample", y = "% of cells", title = "DN cells passing ATAC QC per sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "top")

# Save
# dev.size()
# width_original = 9.56
# height_original= 4.73
# 
# plot_name <- "chapter3_DNs_atac_barPlot_atac_passing.pdf"
# 
# pdf(plot_name, width = width_original, height =height_original)
# barPlot_atac_passing
# dev.off()

table_passing <- df2 %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::summarise(total = n(), 
                   passingQC = sum(log10nFrags > 3 & TSSEnrichment > 4), 
                   percentage = round((passingQC/total)*100))

passQCplot <- df2 %>% dplyr::mutate(PassQC = ifelse(log10nFrags >= 3 & TSSEnrichment >= 4, "yes", "no")) %>% dplyr::select(Sample, PassQC)
passQCplot.cast <- melt(dcast(passQCplot, Sample~PassQC))
passQCplot.cast$variable <- ifelse(passQCplot.cast$variable == "yes", passQCplot.cast$Sample, "no")
passQCplot.cast$Sample <- factor(passQCplot.cast$Sample, levels = samples)

passQCplot_colors <- c(sample_colors, no="gray")

passQCplot.plot <- ggplot(passQCplot.cast) + 
  geom_bar(aes(y = value, x = Sample, fill = variable), stat="identity") +
  geom_text(data=passQCplot.cast, aes(x = Sample, y = as.numeric(value), label = value, size=4), position = position_stack(vjust = 0.5), show.legend = F) + 
  scale_fill_manual(values= passQCplot_colors) + theme_bw()
passQCplot.plot

# Filter out bad cells:
proj_DNs.filtered <- proj_DNs[proj_DNs$TSSEnrichment >= 4 & proj_DNs$nFrags >= 1000]
proj <- proj_DNs.filtered

# Add (RNA) Group information 
proj.metadata <- as.data.frame(proj@cellColData) %>% rownames_to_column("cellNames")
DNs_fromRNA <- seu.VTA_DNs@meta.data %>% dplyr::select(orig.ident, simpleIdent, seurat_clusters, sample_prep, FACS, GEXlibrary, Sequencing) %>% rownames_to_column("cellNames")
proj.metadata <- left_join(proj.metadata,DNs_fromRNA, by = "cellNames") 


# ========== 3) Compare QC metrics after filtering (only DNs) ==========
ATACmetrics <- as.data.frame(proj@cellColData)
ATACmetrics$Sample <- factor(ATACmetrics$Sample, levels = sample_ids)
ATACmetrics.TSS <- ggplot(ATACmetrics, aes(x=Sample, y=TSSEnrichment, fill = Sample, colour = Sample))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  geom_boxplot(aes(x = as.numeric(Sample)+0.25, y = TSSEnrichment),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_y_log10() +
  labs(title = "DNs: TSS enrichment", x= "", y= "") +
  scale_fill_manual(values = sample_colors) + 
  scale_color_manual(values = sample_colors) +
  theme_minimal() +
  theme(legend.position = "none")

ATACmetrics.nFrags <- ggplot(ATACmetrics, aes(x=Sample, y=nFrags, fill = Sample, colour = Sample))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  geom_boxplot(aes(x = as.numeric(Sample)+0.25, y = nFrags),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_y_log10() +
  labs(title = "DNs: Number of unique fragments", x= "", y= "") +
  scale_fill_manual(values = sample_colors) + 
  scale_color_manual(values = sample_colors) +
  theme_minimal() +
  theme(legend.position = "none") 
ATACmetrics.TSS /ATACmetrics.nFrags

DNs_fromRNA <-  seu.VTA_DNs@meta.data %>%
  tibble::rownames_to_column("cellNames") %>% 
  dplyr::select(cellNames, orig.ident, simpleIdent, seurat_clusters, sample_prep, FACS, GEXlibrary, Sequencing) %>%
  dplyr::filter(cellNames %in% proj$cellNames) %>%
  tibble::column_to_rownames("cellNames") 

proj <- addCellColData(ArchRProj = proj,
                       data = as.character(DNs_fromRNA$seurat_clusters),
                       cells = as.character(rownames(DNs_fromRNA)),
                       name = "seurat_clusters")


#proj is the object with the final set of VTA-DNs, with high quality cells (in total, 1353 VTA-DNs)

#Last QC plot
# ========== Classify cells by ATAC status ==========

atac_all     <- getCellNames(proj_DNs)
atac_passing <- getCellNames(proj_DNs.filtered)

# Classify each Seurat cell
seu.VTA_DNs$ATAC_status <- case_when(
  colnames(seu.VTA_DNs) %in% atac_passing ~ "Pass QC",
  colnames(seu.VTA_DNs) %in% atac_all     ~ "Low quality ATAC",
  TRUE                                     ~ "Not detected in ATAC"
)

# Check counts
table(seu.VTA_DNs$ATAC_status)

# ========== DimPlot ==========
seu.VTA_DNs_atac_info <-DimPlot(seu.VTA_DNs, 
                                group.by = "ATAC_status",
                                order    = c("Not detected in ATAC", "Low quality ATAC", "Pass QC"),
                                cols     = c("Pass QC"               = "gray80",
                                             "Low quality ATAC"      = "#2980b9",
                                             "Not detected in ATAC"  = "#c0392b")) +
  labs(title = "ATAC quality per DN cell") +
  theme_void()

# dev.size()
# width_original = 9.03
# height_original= 6.35
# 
# plot_name <- "chapter3_seu.VTA_DNs_atac_info.pdf"
# 
# pdf(plot_name, width = width_original, height =height_original)
# seu.VTA_DNs_atac_info
# dev.off()


#saveRDS(proj_DNs.filtered, "proj_DNs.filtered_backup.rds")



# ========== 4) Process ATAC data needed for peak calling ==========
# cells_peakcalling <- getCellNames(proj_DNs.filtered)[
#   proj_DNs.filtered$Sample %notin% c("h1_saline_R1", "m30_cocaine_R1")
# ]
# proj_DNs.peakcalling <- proj_DNs.filtered[cells_peakcalling, ]
# 
# # sanity checks
# table(proj_DNs.peakcalling$Sample)  # 16 samples
# sum(table(proj_DNs.peakcalling$Sample)) # should be 1353 - 41 - 68 = 1244
# 
# # Convert ArchR project to Signac/Seurat object
# seurat_atac <- getMatrixFromProject(
#   ArchRProj  = proj_DNs.peakcalling,
#   useMatrix  = "TileMatrix"
# )

data_dir <- "/fast/AG_Pombo/luna/2026_rebuttal/9_ATAC-QCandIntegration/data_dir/"
dir.create("DN_fragments")
sample_ids.peakcalling <- sample_ids[sample_ids %notin% c("h1_saline_R1", "m30_cocaine_R1")]

for (sample in sample_ids.peakcalling) {
  cat("Processing:", sample, "\n")
  
  input_file  <- glue("{data_dir}{sample}.atac_fragments.tsv.gz")
  temp_file   <- glue("DN_fragments/{sample}.temp.tsv")
  output_file <- glue("DN_fragments/{sample}.atac_fragments.tsv.gz")
  
  input_data <- data.table::fread(input_file, header = FALSE) %>%
    dplyr::mutate(V4 = glue("{sample}#{V4}")) %>%
    dplyr::filter(V4 %in% getCellNames(proj_DNs.peakcalling))
  
  write.table(input_data, file = temp_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Compress and index with Rsamtools
  bgzip(temp_file, dest = output_file, overwrite = TRUE)
  indexTabix(output_file, format = "bed")
  file.remove(temp_file)
  
  cat("Done:", sample, "->", nrow(input_data), "fragments\n")
}
#load("260531.backup.rds")

# ========== Create fragment objects per sample ==========
fragment_list <- list()
for (sample in sample_ids.peakcalling) {
  frag_path    <- glue("DN_fragments/{sample}.atac_fragments.tsv.gz")
  sample_cells <- getCellNames(proj_DNs.peakcalling)[grepl(sample, getCellNames(proj_DNs.peakcalling))]
  frag         <- CreateFragmentObject(path = frag_path, cells = sample_cells,
                                       verbose = FALSE, tolerance = 0.5, validate.fragments = FALSE)
  fragment_list <- append(fragment_list, frag)
}
names(fragment_list) <- sample_ids.peakcalling


# ========== 5. Peak calling through Signac  ==========

pathToMacs2 <- "/fast/home/l/lzeared/profiles/r_multiome230913/bin/macs2"

# ------------------------------------------------------------
# 5.1 Check fragment list
# ------------------------------------------------------------

fragment_list
sapply(fragment_list, class)
sapply(fragment_list, function(x) length(Cells(x)))

# Extract fragment file paths from Fragment objects
fragment_paths <- sapply(fragment_list, function(x) {
  GetFragmentData(object = x, slot = "path")
})

fragment_paths
file.exists(fragment_paths)
file.exists(paste0(fragment_paths, ".tbi"))


# ------------------------------------------------------------
# 5.2 Call peaks — pooled across all DNs
# ------------------------------------------------------------
pathToMacs2 <- "/fast/home/l/lzeared/profiles/r_multiome230913/bin/macs2"

# Call MACS2 using all fragment files together
peaks_pool <- CallPeaks(
  object = fragment_paths,
  macs2.path = pathToMacs2,
  effective.genome.size = 1.87e9,
  name = "DNs_pool"
)

peaks_pool <- keepStandardChromosomes(
  peaks_pool,
  pruning.mode = "coarse"
)

saveRDS(peaks_pool, "260601_peaks_pool.rds")


# ------------------------------------------------------------
# 5.3 Call peaks — per condition
# ------------------------------------------------------------

#Assign each fragment file to a sampleIdent group

condition_by_sample <- rep(NA_character_, length(fragment_paths))
names(condition_by_sample) <- names(fragment_paths)

# Cocaine timepoints
condition_by_sample[startsWith(names(fragment_paths), "h1_cocaine")]  <- "h1_cocaine"
condition_by_sample[startsWith(names(fragment_paths), "h4_cocaine")]  <- "h4_cocaine"
condition_by_sample[startsWith(names(fragment_paths), "h8_cocaine")]  <- "h8_cocaine"
condition_by_sample[startsWith(names(fragment_paths), "h24_cocaine")] <- "h24_cocaine"
condition_by_sample[startsWith(names(fragment_paths), "d14_cocaine")] <- "d14_cocaine"
condition_by_sample[startsWith(names(fragment_paths), "m30_cocaine")] <- "m30_cocaine"

# Pool all saline samples together
condition_by_sample[grepl("saline", names(fragment_paths), ignore.case = TRUE)] <- "saline"

# Check assignment
print(condition_by_sample)
print(table(condition_by_sample, useNA = "ifany"))

if (any(is.na(condition_by_sample))) {
  print(names(fragment_paths)[is.na(condition_by_sample)])
  stop("Some fragment paths could not be assigned to a sampleIdent group.")
}



#Split fragment paths by sampleIdent group
fragment_paths_by_condition <- split(fragment_paths, condition_by_sample)

# Check how many fragment files per condition
print(sapply(fragment_paths_by_condition, length))
print(fragment_paths_by_condition)


peaks_per_condition <- lapply(names(fragment_paths_by_condition), function(condition_name) {
  message("Calling peaks for condition: ", condition_name)
  
  peaks <- CallPeaks(
    object = unname(fragment_paths_by_condition[[condition_name]]),
    macs2.path = pathToMacs2,
    effective.genome.size = 1.87e9,
    name = paste0("DNs_", condition_name)
  )
  
  peaks <- keepStandardChromosomes(
    peaks,
    pruning.mode = "coarse"
  )
  
  return(peaks)
})

names(peaks_per_condition) <- names(fragment_paths_by_condition)

# Check peak counts per condition
print(sapply(peaks_per_condition, length))

saveRDS(peaks_per_condition, "260601_peaks_per_condition.rds")


# ------------------------------------------------------------
# 5.4 Add condition-specific peaks missing from pooled peak set
# ------------------------------------------------------------
# Label pooled peaks
peaks_pool$peak_found_in <- "pool"

# For each condition, find peaks that do NOT overlap the pooled peak set
peaks_condition_only <- lapply(names(peaks_per_condition), function(condition_name) {
  
  message("Checking condition-specific peaks for: ", condition_name)
  
  cond_peaks <- peaks_per_condition[[condition_name]]
  
  cond_peaks <- keepStandardChromosomes(
    cond_peaks,
    pruning.mode = "coarse"
  )
  
  # Keep only condition peaks with no overlap to pooled peaks
  no_pool_overlap <- !overlapsAny(cond_peaks, peaks_pool)
  cond_only <- cond_peaks[no_pool_overlap]
  
  # Annotate
  cond_only$peak_found_in <- paste0(condition_name, "_only")
  
  message(condition_name, ": ", length(cond_only), " peaks not found in pool")
  
  return(cond_only)
})

names(peaks_condition_only) <- names(peaks_per_condition)

# Check how many condition-only peaks will be added
sapply(peaks_condition_only, length)


# ------------------------------------------------------------
# 5.5 Build final peak universe
# ------------------------------------------------------------

# Combine pool peaks plus condition-only peaks
all_peaks_with_labels <- c(
  peaks_pool,
  unlist(GRangesList(peaks_condition_only), use.names = FALSE)
)

# Reduce to a non-overlapping final peak universe
peaks_final <- reduce(all_peaks_with_labels)

# Annotate final peaks:
# If a final peak overlaps the pooled set, label as "pool".
# Otherwise, label by the condition-specific peak(s) it overlaps.

# Annotate final peaks
peak_found_in <- rep(NA_character_, length(peaks_final))

# Mark peaks overlapping the pooled peak set
pool_hits <- overlapsAny(peaks_final, peaks_pool)
peak_found_in[pool_hits] <- "pool"

# Mark peaks that do not overlap pool but overlap condition-only peaks
condition_only_all <- unlist(GRangesList(peaks_condition_only), use.names = FALSE)

condition_hits <- findOverlaps(peaks_final, condition_only_all)

for (i in seq_along(peaks_final)) {
  
  # Skip pooled peaks
  if (isTRUE(peak_found_in[i] == "pool")) {
    next
  }
  
  hits_i <- subjectHits(condition_hits)[queryHits(condition_hits) == i]
  
  if (length(hits_i) > 0) {
    peak_found_in[i] <- paste(
      unique(condition_only_all$peak_found_in[hits_i]),
      collapse = ";"
    )
  }
}

peaks_final$peak_found_in <- peak_found_in

# Check final annotation
table(peaks_final$peak_found_in, useNA = "ifany")

length(peaks_pool)
sapply(peaks_condition_only, length)
length(peaks_final)

saveRDS(peaks_final, "260601_peaks_pool_plus_condition_only.rds")

# ------------------------------------------------------------
# 5.6 Export peak table with peak_found_in column
# ------------------------------------------------------------

peak_table <- data.frame(
  peak_id = paste0(
    seqnames(peaks_final), "-",
    start(peaks_final), "-",
    end(peaks_final)
  ),
  seqnames = as.character(seqnames(peaks_final)),
  start = start(peaks_final),
  end = end(peaks_final),
  width = width(peaks_final),
  peak_found_in = peaks_final$peak_found_in
)

head(peak_table)
table(peak_table$peak_found_in, useNA = "ifany")

write.csv(peak_table,"260601_peaks_pool_plus_condition_only_table.csv",row.names = FALSE)

#save.image("260601_ATAC_QC_Integration_PeakCalling.rds")


# ----------------------------------------
# 6 Creat signac object with peak matrix
# ----------------------------------------

# Use final non-overlapping peak universe from section 5
peaks_use <- peaks_final

# Keep standard chromosomes and remove problematic ranges if needed
peaks_use <- keepStandardChromosomes(peaks_use, pruning.mode = "coarse")
peaks_use <- sort(peaks_use)

# Create stable peak IDs
peak_ids <- paste0(seqnames(peaks_use), "-", start(peaks_use), "-", end(peaks_use))
names(peaks_use) <- peak_ids

# Check fragment objects and paths
stopifnot(exists("fragment_list"))
stopifnot(length(fragment_list) > 0)

fragment_paths <- sapply(fragment_list, function(x) GetFragmentData(object = x, slot = "path"))
stopifnot(all(file.exists(fragment_paths)))
stopifnot(all(file.exists(paste0(fragment_paths, ".tbi"))))

# Name samples if needed
if (is.null(names(fragment_list)) || any(names(fragment_list) == "")) names(fragment_list) <- paste0("sample_", seq_along(fragment_list))
sample_ids <- names(fragment_list)

# Create one Signac/Seurat object per sample
objs <- lapply(sample_ids, function(sample_id) {
  message("Creating peak matrix for: ", sample_id)
  
  frag_obj <- fragment_list[[sample_id]]
  frag_path <- GetFragmentData(object = frag_obj, slot = "path")
  cell_barcodes <- Cells(frag_obj)
  
  counts <- FeatureMatrix(fragments = frag_obj, features = peaks_use, cells = cell_barcodes)
  rownames(counts) <- peak_ids
  
  chrom_assay <- CreateChromatinAssay(counts = counts, ranges = peaks_use, fragments = frag_path, genome = "mm10", min.cells = 1, min.features = 1)
  obj <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC", project = sample_id)
  
  obj$sample <- sample_id
  if (exists("condition_by_sample")) obj$condition <- rep(unname(condition_by_sample[sample_id]), ncol(obj))
  
  return(obj)
})

names(objs) <- sample_ids

# Merge all samples into one ATAC object
atac <- merge(x = objs[[1]], y = objs[-1], add.cell.ids = names(objs))

# Set default assay
DefaultAssay(atac) <- "ATAC"

# Basic checks
atac
dim(atac)
table(atac$sample)
if ("condition" %in% colnames(atac@meta.data)) table(atac$condition, useNA = "ifany")

# Save object
# saveRDS(atac, "260601_signac_atac_peak_matrix_object.rds")


# ----------------------------------------
# 7. Annotation of peaks 
# ----------------------------------------
atac <- readRDS("260601_signac_atac_peak_matrix_object.rds")

# Get peaks from the Seurat/Signac ATAC assay
peaks <- granges(atac[["ATAC"]])
names(peaks) <- rownames(atac[["ATAC"]])

# Keep original peak names so metadata can be aligned back to Seurat
all_peak_names <- names(peaks)

# Get ArchR annotations from your ArchR project
geneAnnotation <- getGeneAnnotation(proj_DNs.filtered)
genomeAnnotation <- getGenomeAnnotation(proj_DNs.filtered)

# Your custom promoter definition
promoterRegion <- c(1000, 500)

# Validate inputs, ArchR-style
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- ArchR:::.validGRanges(peaks)

peakSummits <- GenomicRanges::resize(peaks, 1, "center")

geneAnnotation$genes <- ArchR:::.validGRanges(geneAnnotation$genes)
geneAnnotation$exons <- ArchR:::.validGRanges(geneAnnotation$exons)
geneAnnotation$TSS <- ArchR:::.validGRanges(geneAnnotation$TSS)

BSgenome <- eval(parse(text = genomeAnnotation$genome))
BSgenome <- ArchR:::validBSgenome(BSgenome)

message("Annotating Peaks : Nearest Gene")

distPeaks <- distanceToNearest(
  peakSummits,
  GenomicRanges::resize(geneAnnotation$genes, 1, "start"),
  ignore.strand = TRUE
)

distToGeneStart <- rep(NA_integer_, length(peaks))
nearestGene <- rep(NA_character_, length(peaks))

distToGeneStart[queryHits(distPeaks)] <- mcols(distPeaks)$distance
nearestGene[queryHits(distPeaks)] <-
  mcols(geneAnnotation$genes)$symbol[subjectHits(distPeaks)]


message("Annotating Peaks : Gene")

promoters <- ArchR:::extendGR(
  GenomicRanges::resize(geneAnnotation$genes, 1, "start"),
  upstream = promoterRegion[1],
  downstream = promoterRegion[2]
)

op <- overlapsAny(peakSummits, promoters, ignore.strand = TRUE)
og <- overlapsAny(peakSummits, geneAnnotation$genes, ignore.strand = TRUE)
oe <- overlapsAny(peakSummits, geneAnnotation$exons, ignore.strand = TRUE)

peakType <- rep("Distal", length(peaks))
peakType[which(og & oe)] <- "Exonic"
peakType[which(og & !oe)] <- "Intronic"
peakType[which(op)] <- "Promoter"


message("Annotating Peaks : TSS")

distTSS <- distanceToNearest(
  peakSummits,
  GenomicRanges::resize(geneAnnotation$TSS, 1, "start"),
  ignore.strand = TRUE
)

distToTSS <- rep(NA_integer_, length(peaks))
nearestTSS <- rep(NA_character_, length(peaks))

distToTSS[queryHits(distTSS)] <- mcols(distTSS)$distance

if ("symbol" %in% colnames(mcols(geneAnnotation$TSS))) {
  nearestTSS[queryHits(distTSS)] <-
    mcols(geneAnnotation$TSS)$symbol[subjectHits(distTSS)]
} else if ("tx_name" %in% colnames(mcols(geneAnnotation$TSS))) {
  nearestTSS[queryHits(distTSS)] <-
    mcols(geneAnnotation$TSS)$tx_name[subjectHits(distTSS)]
}


message("Annotating Peaks : GC")

nucFreq <- BSgenome::alphabetFrequency(
  BSgenome::getSeq(BSgenome, peaks)
)

GC <- round(rowSums(nucFreq[, c("G", "C"), drop = FALSE]) / rowSums(nucFreq), 4)
N <- round(nucFreq[, "N"] / rowSums(nucFreq), 4)

peak_metadata <- data.frame(
  peakType = peakType,
  distToGeneStart = distToGeneStart,
  nearestGene = nearestGene,
  distToTSS = distToTSS,
  nearestTSS = nearestTSS,
  GC = GC,
  N = N,
  row.names = names(peaks)
)

# Full metadata table aligned to the original Seurat feature order
# This avoids problems if non-standard chromosomes were removed above.
n_peaks <- length(all_peak_names)

full_metadata <- data.frame(
  peakType = rep(NA_character_, n_peaks),
  distToGeneStart = rep(NA_integer_, n_peaks),
  nearestGene = rep(NA_character_, n_peaks),
  distToTSS = rep(NA_integer_, n_peaks),
  nearestTSS = rep(NA_character_, n_peaks),
  GC = rep(NA_real_, n_peaks),
  N = rep(NA_real_, n_peaks),
  row.names = all_peak_names
)

full_metadata[rownames(peak_metadata), colnames(peak_metadata)] <- peak_metadata

# Add to ATAC feature metadata
atac[["ATAC"]] <- AddMetaData(
  object = atac[["ATAC"]],
  metadata = full_metadata
)

saveRDS(atac, "260608_signac_atac_peak_matrix_object_annotated.rds")

#Plots

# ----------------------------------------
# 7. Peaktype and GC content: distributions
# ----------------------------------------
all_peaks <- atac[["ATAC"]]@meta.features

# Add peak IDs from rownames
all_peaks <- all_peaks %>%
  tibble::rownames_to_column("peakID")

# Check required columns exist
stopifnot("peakType" %in% colnames(all_peaks))
stopifnot("GC" %in% colnames(all_peaks))

# Keep only peaks with annotation and GC
all_peaks <- all_peaks %>%
  dplyr::filter(!is.na(peakType), !is.na(GC))

# Set peakType order
all_peaks <- all_peaks %>%
  dplyr::mutate(
    peakType = factor(peakType, levels = c("Distal", "Intronic", "Exonic", "Promoter")))

# Define colors
peak_type_colors <- c(Promoter = "#8B0000",Exonic   = "#DAA520",Intronic = "#2E8B57",Distal   = "#4682B4")


#Peak type distribution

summary_data <- all_peaks %>%
  dplyr::distinct(peakID, .keep_all = TRUE) %>%
  dplyr::group_by(peakType) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(
    percentage = round((count / sum(count)) * 100, 1),
    label = paste0(percentage, "% (n = ", count, ")"))

bar_chart <- ggplot(
  summary_data,
  aes(x = count, y = peakType, fill = peakType)
) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(
    aes(label = label),
    hjust = -0.1,
    size = 3.5
  ) +
  scale_fill_manual(values = peak_type_colors, drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  coord_cartesian(xlim = c(0, max(summary_data$count) * 1.2)) +
  theme_classic() +
  labs(
    title = "Peak Type Distribution",
    x = "Count",
    y = "Peak Type"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

# GC content per peak type
gc_per_peakType <- all_peaks %>%
  dplyr::select(peakID, GC, peakType) %>%
  dplyr::distinct(peakID, .keep_all = TRUE) %>%
  dplyr::mutate(
    peakType = factor(
      peakType,
      levels = c("Distal", "Intronic", "Exonic", "Promoter")
    )
  )

gc_per_peakType.plot <- ggplot(
  gc_per_peakType,
  aes(x = GC, y = peakType, fill = peakType, colour = peakType)
) +
  geom_point(
    position = position_jitter(height = 0.15),
    size = 0.25,
    alpha = 0.6
  ) +
  geom_boxplot(
    aes(x = GC, y = as.numeric(peakType) + 0.25),
    outlier.shape = NA,
    alpha = 0.3,
    width = 0.1,
    colour = "black"
  ) +
  labs(
    title = "",
    x = "GC Content",
    y = "Peak Type"
  ) +
  scale_fill_manual(values = peak_type_colors, drop = FALSE) +
  scale_color_manual(values = peak_type_colors, drop = FALSE) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# Combine plots
combined_plot <- bar_chart + gc_per_peakType.plot
combined_plot

# Save PDF
plot_name <- "chapter3_allPeaks_number_and_GCcontent.pdf"

pdf(file = plot_name, width = 8.968750, height = 4.666667 )
combined_plot
dev.off()



gc_median_per_peakType <- all_peaks %>%
  dplyr::filter(!is.na(peakType), !is.na(GC)) %>%
  dplyr::group_by(peakType) %>%
  dplyr::summarise(
    median_GC = median(GC, na.rm = TRUE),
    n_peaks = dplyr::n(),
    .groups = "drop"
  )

gc_median_per_peakType

