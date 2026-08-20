#!/usr/bin/env Rscript
# ===========================================
# Script Title: VTA vs SN ATAC Peak Calling
# Author: Luna Zea Redondo
# Date: 2023-11-13
#
# Description:
#   This script performs ATAC-seq peak calling to compare VTA and SN dopaminergic neurons.
#   Peak calling is performed using Signac and MACS2, starting from fragment files
#   generated during the ArchR preprocessing step (in script 3.4).
#
#   The workflow includes:
#     - Assignment of cells to VTA or SN based on RNA-derived metadata
#     - Peak calling stratified by brain region (VTA vs SN)
#     - Construction of region-specific peak matrices
#     - Generation and export of independent peak sets for VTA and SN neurons
#
#   Two alternative strategies for peak calling are explored:
#     1. Region-based peak calling using a single object (active)
#     2. Independent peak calling per subset using modified fragment files (standby)
#
#   The resulting peak sets are intended for downstream ArchR integration
#   and differential accessibility analyses.
# ===========================================

# ========== 0. Environment Setup ==========
rm(list = ls(all.names = TRUE))
gc()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

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


# ========== 2. Define paths and input data ==========
data_dir <- "/data/pombo/Luna/MultiomeCocaineTreatments/data_dir/"
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/final_peaksets/231113_VTAvsSN/231113_signac_peakCall"
setwd(dir)

load("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR/230626.ATAC_DN_condition_backup.rds") #contains wnn.seu, fragments_dir, etc. 

samples <- c("m30_cocaine_R1", "h1_cocaine_R1", "h1_cocaine_R2",   #h1_saline is removed
             "h4_saline_R1", "h4_cocaine_R1", "h4_cocaine_R2",
             "h8_saline_R1", "h8_cocaine_R1", "h8_cocaine_R2",
             "h24_saline_R1", "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3",
             "d14_saline_R1", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3")

# Load DN-VTA metadata
VTA_DN_metadata <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002_VTA_RNAonly_metadata_1570cells.tsv")
wnn.seu$region <- ifelse(colnames(wnn.seu) %in% VTA_DN_metadata$cellNames, "VTA", "SN")

pathToMacs2 <- "/home/lzeared/profiles/r_multiome230913/bin/macs2"

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# ========== 3. Assign Brain Region Labels ==========
# Identify VTA and SN cells using RNA metadata
# Annotate Seurat object with region information
wnn.seu@assays[["peaks"]] <- NULL
Fragments(object = wnn.seu) <- fragment_list 
DefaultAssay(wnn.seu) <- "ATAC"
VTA_vs_SN_peaks.list <- CallPeaks(
  object = wnn.seu,
  group.by = "region",
  combine.peaks = FALSE, 
  macs2.path = pathToMacs2,
  effective.genome.size = 1.87e9)
names(VTA_vs_SN_peaks.list) <- c("VTA", "SN")

# findOverlaps(VTA_vs_SN_peaks.list$VTA, VTA_vs_SN_peaks.list$SN)
# findOverlaps(VTA_vs_SN_peaks.list$SN, VTA_vs_SN_peaks.list$VTA)


# ========== 4. region-based peak calling (active strategy) ==========
# Remove existing peak assay
# Assign fragment list
# Call peaks per region using MACS2
# Generate separate peak matrices for VTA and SN

#Divide wnn.seu object now:
wnn.seu.VTA.231113 <- subset(wnn.seu, subset = region == "VTA")
wnn.seu.SN.231113 <- subset(wnn.seu, subset = region == "SN")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
wnn.seu.VTA.peaks.231113 <- keepStandardChromosomes(VTA_vs_SN_peaks.list$VTA, pruning.mode = "coarse")
wnn.seu.SN.peaks.231113 <- keepStandardChromosomes(VTA_vs_SN_peaks.list$SN, pruning.mode = "coarse")

# quantify counts in each peak
VTA_macs2_counts <- FeatureMatrix(
  fragments = Fragments(wnn.seu.VTA.231113),
  features = wnn.seu.VTA.peaks.231113,
  cells = colnames(wnn.seu.VTA.231113)
)

SN_macs2_counts <- FeatureMatrix(
  fragments = Fragments(wnn.seu.SN.231113),
  features = wnn.seu.SN.peaks.231113,
  cells = colnames(wnn.seu.SN.231113)
)

#Add peakset to each object
wnn.seu.VTA.231113[["peaks"]] <- CreateChromatinAssay(
  counts = VTA_macs2_counts,
  fragments = fragment_list,
  annotation = annotation
)

wnn.seu.SN.231113[["peaks"]] <- CreateChromatinAssay(
  counts = SN_macs2_counts,
  fragments = fragment_list,
  annotation = annotation
)

# ========== 5. Save Region-Specific Peak Sets ==========
# Export VTA peak GRanges
# Export SN peak GRanges
# Sanity checks
VTApeaks <- as.data.frame(wnn.seu.VTA.231113@assays[["peaks"]]@ranges)
saveRDS(VTApeaks, "231113.VTApeaks.rds")
SNpeaks <- as.data.frame(wnn.seu.SN.231113@assays[["peaks"]]@ranges)
saveRDS(SNpeaks, "231113.SNpeaks.rds")

#Sanity checks:
identical(VTApeaks, SNpeaks)
identical(wnn.seu.SN.231113, wnn.seu.VTA.231113)

#save.image("231113_peakCall_VTAvsSN_DNs_SIGNAC.rds")


# ========== 6. Alternative Peak Calling Strategies (Standby) ==========
# Fragment file subsetting per region
# Independent peak calling per subset

wnn.seu.VTA.231113 <- subset(wnn.seu, subset = region == "VTA")
wnn.seu.SN.231113 <- subset(wnn.seu, subset = region == "SN")

#Annotation and fragment files:
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# original fragment_list is loaded in line 52

#VTA fragments:
VTA_fragments_dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/final_peaksets/231113_VTAvsSN/231113_signac_peakCall/VTA_modified_fragments"
for (sample in samples) {
  # Analysis of sample
  print(glue("Processing fragments of sample: {sample}"))
  # Step 1: Set the file paths
  input_file <- glue("{data_dir}{sample}.atac_fragments.tsv.gz")
  temp_file <- glue("{VTA_fragments_dir}/{sample}.temp.tsv")
  output_file <- glue("{VTA_fragments_dir}/{sample}.atac_fragments.tsv")

  # Step 2: Extract atac fragments header and saved as temporary file:
  header_lines <- as.data.frame(readLines(gzfile(input_file), n = 52))
  write.table(header_lines, "temp.header.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Step 3: Read and modify the fragment file to keep only DN cells
  input_data <- fread(input_file, header = FALSE)  %>%
    mutate(V4 = glue("{sample}#{V4}")) %>%
    filter(V4 %in% colnames(wnn.seu.VTA.231113))

  # Step 4: Write the fragment file as temporary file
  write.table(input_data, file = temp_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Step 5: Concatenate header and DN fragments
  system(glue("cat temp.header.tsv {temp_file}> {output_file}"))

  # Step 6: Compress the output file
  system(glue("bgzip {output_file}"))

  # Step 7: generate index:
  system(glue("tabix -p bed {output_file}.gz"))

  # Step 8: remove temporary files
  system(glue("rm {temp_file}"))
}


#SN fragments:
SN_fragments_dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/final_peaksets/231113_VTAvsSN/231113_signac_peakCall/SN_modified_fragments"
for (sample in samples) {
  # Analysis of sample
  print(glue("Processing fragments of sample: {sample}"))
  # Step 1: Set the file paths
  input_file <- glue("{data_dir}{sample}.atac_fragments.tsv.gz")
  temp_file <- glue("{SN_fragments_dir}/{sample}.temp.tsv")
  output_file <- glue("{SN_fragments_dir}/{sample}.atac_fragments.tsv")

  # Step 2: Extract atac fragments header and saved as temporary file:
  header_lines <- as.data.frame(readLines(gzfile(input_file), n = 52))
  write.table(header_lines, "temp.header.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Step 3: Read and modify the fragment file to keep only DN cells
  input_data <- fread(input_file, header = FALSE)  %>%
    mutate(V4 = glue("{sample}#{V4}")) %>%
    filter(V4 %in% colnames(wnn.seu.SN.231113))

  # Step 4: Write the fragment file as temporary file
  write.table(input_data, file = temp_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Step 5: Concatenate header and DN fragments
  system(glue("cat temp.header.tsv {temp_file}> {output_file}"))

  # Step 6: Compress the output file
  system(glue("bgzip {output_file}"))

  # Step 7: generate index:
  system(glue("tabix -p bed {output_file}.gz"))

  # Step 8: remove temporary files
  system(glue("rm {temp_file}"))
}




## VTA-DN cells:

#Read fragment list
VTA_fragment_list <- list()
for (sample in samples) {
  frag.sample.path = glue("{VTA_fragments_dir}/{sample}.atac_fragments.tsv.gz")
  sample.cells <- wnn.seu.VTA.231113[, grepl(sample, colnames(wnn.seu.VTA.231113))]
  frag.sample <- CreateFragmentObject(path = frag.sample.path, cells = colnames(sample.cells), verbose = FALSE, tolerance = 0.5, validate.fragments = FALSE)
  VTA_fragment_list <- append(VTA_fragment_list, frag.sample)
}
names(VTA_fragment_list) <- samples  

# call peaks using MACS2
wnn.seu.VTA.231113@assays[["peaks"]] <- NULL

DefaultAssay(wnn.seu.VTA.231113) <- "ATAC"
wnn.seu.VTA.peaks.231113 <- CallPeaks(
  object = VTA_fragment_list,
  group.by = NULL,
  macs2.path = pathToMacs2,
  effective.genome.size = 1.87e9)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
wnn.seu.VTA.peaks.231018 <- keepStandardChromosomes(wnn.seu.VTA.peaks.231018, pruning.mode = "coarse")

# quantify counts in each peak
VTA_macs2_counts <- FeatureMatrix(
  fragments = Fragments(wnn.seu.VTA.231018),
  features = wnn.seu.VTA.peaks.231018,
  cells = colnames(wnn.seu.VTA.231018)
)

wnn.seu.VTA.231018[["peaks"]] <- CreateChromatinAssay(
  counts = VTA_macs2_counts,
  fragments = fragment_list,
  annotation = annotation
)



## SN-DN cells:

# call peaks using MACS2
DefaultAssay(wnn.seu.SN.231018) <- "peaks"
wnn.seu.SN.231018@assays[["peaks"]] <- NULL
wnn.seu.SN.peaks.231018 <- CallPeaks(
  object = wnn.seu.SN.231018,
  group.by = NULL,
  macs2.path = pathToMacs2,
  effective.genome.size = 1.87e9
)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
wnn.seu.SN.peaks.231018 <- keepStandardChromosomes(wnn.seu.SN.peaks.231018, pruning.mode = "coarse")

# quantify counts in each peak
SN_macs2_counts <- FeatureMatrix(
  fragments = Fragments(wnn.seu.SN.231018),
  features = wnn.seu.SN.peaks.231018,
  cells = colnames(wnn.seu.SN.231018)
)

#save.image("231019.temporary.rds")
wnn.seu.SN.231018[["peaks"]] <- CreateChromatinAssay(
  counts = SN_macs2_counts,
  fragments = fragment_list,
  annotation = annotation
)


# ========== 7. Output ==========
# RDS files with region-specific peak sets
# Objects for downstream ArchR and Signac analyses

#Save the peak files independently to be added to ArchR:
VTApeaks <- as.data.frame(wnn.seu.VTA.231018@assays[["peaks"]]@ranges)
saveRDS(VTApeaks, "231019.VTApeaks.rds")
SNpeaks <- as.data.frame(wnn.seu.SN.231018@assays[["peaks"]]@ranges)
saveRDS(SNpeaks, "231019.SNpeaks.rds")

identical(wnn.seu.SN.231018, wnn.seu.VTA.231018) #sanity check
#save.image("231018_peakCall_VTAvsSN_DNs_SIGNAC.rds")



# End of script
