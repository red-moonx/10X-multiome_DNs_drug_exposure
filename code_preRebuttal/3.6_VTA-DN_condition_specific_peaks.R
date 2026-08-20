# ===========================================
# Script Title: condition-Specific ATAC Peaks (VTA-DNs)
# Author: Luna Zea Redondo
# Date: 2023-11-23
#
# Description:
#   This script identifies and evaluates condition-specific ATAC-seq peaks
#   in VTA dopaminergic neurons.
#
#   The workflow includes:
#     - Condition-specific peak calling using Signac + MACS2
#     - Annotation of peak sets with genomic context and motifs
#     - Comparison of peak sets across cocaine timepoints
#     - Construction of a pooled and condition-exclusive peak catalogue
#     - Peak length, GC content, and RPKM-based accessibility analysis
#
#   This script focuses on peak-level characterization and comparison.
#   Downstream differential accessibility and regulatory analyses
#   are performed in separate scripts.
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
library(ggplot2)             
library(glue)               
library(tibble)              
library(readr)               
library(patchwork)           
library(ComplexHeatmap)      
library(scater)              
library(limma)               
library(UpSetR)              


# ========== 2. Define paths and input Files ==========

source("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/signac/231010/ValidationUtils.R")

data_dir <- "/data/pombo/Luna/MultiomeCocaineTreatments/data_dir/"
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/6_231107_updated_analyses/2_ATAC_alternative_processing/2_231123_condition_specific_peaks"
setwd(dir)

load("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR/230626.ATAC_DN_condition_backup.rds")


# ========== 3. Load Preprocessed Objects ==========
# WNN Seurat object
# ArchR VTA project
# RNA-based VTA DN metadata

samples <- c("m30_cocaine_R1", "h1_cocaine_R1", "h1_cocaine_R2",   #h1_saline is removed
             "h4_saline_R1", "h4_cocaine_R1", "h4_cocaine_R2",
             "h8_saline_R1", "h8_cocaine_R1", "h8_cocaine_R2",
             "h24_saline_R1", "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3",
             "d14_saline_R1", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3")

#Load DN-VTA metadata
VTA_DN_metadata <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002_VTA_RNAonly_metadata_1570cells.tsv")
wnn.seu$region <- ifelse(colnames(wnn.seu) %in% VTA_DN_metadata$cellNames, "VTA", "SN")

pathToMacs2 <- "/home/lzeared/profiles/r_multiome230913/bin/macs2"

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))


# ========== 4. Subset cells by brain region ==========
wnn.seu.VTA.231123 <- subset(wnn.seu, subset = region == "VTA")
wnn.seu.SN.231123 <- subset(wnn.seu, subset = region == "SN")


# ========== 5. Condition-Specific Peak Calling ==========

VTA_fragments_dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/final_peaksets/231113_VTAvsSN/231113_signac_peakCall/VTA_modified_fragments"

VTA_fragment_list <- list()
for (sample in samples) {
  frag.sample.path = glue("{VTA_fragments_dir}/{sample}.atac_fragments.tsv.gz")
  sample.cells <- wnn.seu.VTA.231123[, grepl(sample, colnames(wnn.seu.VTA.231123))]
  frag.sample <- CreateFragmentObject(path = frag.sample.path, cells = colnames(sample.cells), verbose = FALSE, tolerance = 0.5, validate.fragments = FALSE)
  VTA_fragment_list <- append(VTA_fragment_list, frag.sample)
}
names(VTA_fragment_list) <- samples  

# call peaks using MACS2
DefaultAssay(wnn.seu.VTA.231123) <- "peaks"
Fragments(wnn.seu.VTA.231123) <- NULL
Fragments(wnn.seu.VTA.231123) <- VTA_fragment_list

wnn.seu.VTA.peaks.231123 <- CallPeaks(
  object = wnn.seu.VTA.231123,
  group.by = "simpleIdent",
  macs2.path = pathToMacs2,
  effective.genome.size = 1.87e9, 
  combine.peaks = FALSE)

VTA.condition.specific.peaks <- wnn.seu.VTA.peaks.231123
names(VTA.condition.specific.peaks) <- unique(wnn.seu.VTA.231123$simpleIdent)

# Load the VTA peakset for the 103K total peaks
proj5.VTA <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/final_peaksets/231018_VTAvsSN/231019_ArchR_DARs_and_motifs/231102.proj5_signac_VTA.annotated.rds")
proj5.VTA.peaks <- getPeakSet(proj5.VTA)

# Calculate condition-specific peaks (in this case, condition is equivalent to timepoint)

for (condition in conditions) {
  message(glue("Analyzing: {condition}")) 
  peaks.condition <- VTA.condition.specific.peaks[[condition]]
  peaks.condition <- keepStandardChromosomes(peaks.condition, pruning.mode = "coarse")
  
  #Add peakset to proj.condition
  proj5.VTA.condition <- proj5.VTA[proj5.VTA$condition == condition]
  proj5.VTA.condition
  
  #General peak annotation:
  proj_condition_macs2_counts <- FeatureMatrix(
    fragments = VTA_fragment_list,
    features = peaks.condition,
    cells = proj5.VTA.condition$cellNames
  )  
  
  #General peak annotation:
  ##########################
  geneAnnotation = getGeneAnnotation(proj5.VTA.condition)
  genomeAnnotation = getGenomeAnnotation(proj5.VTA.condition)
  promoterRegion = c(1000, 500) #I changed this bit
  
  #Validate:
  peaks.condition <- .validGRanges(peaks.condition)
  peakSummits <- GenomicRanges::resize(peaks.condition,1,"center")
  geneAnnotation$genes <- .validGRanges(geneAnnotation$genes)
  geneAnnotation$exons <- .validGRanges(geneAnnotation$exons)
  geneAnnotation$TSS <- .validGRanges(geneAnnotation$TSS)
  BSgenome <- eval(parse(text = genomeAnnotation$genome))
  BSgenome <- validBSgenome(BSgenome)
  
  #Annotate:
  peaks.condition <- .validGRanges(peaks.condition)
  peakSummits <- GenomicRanges::resize(peaks.condition,1,"center")
  geneAnnotation$genes <- .validGRanges(geneAnnotation$genes)
  geneAnnotation$exons <- .validGRanges(geneAnnotation$exons)
  geneAnnotation$TSS <- .validGRanges(geneAnnotation$TSS)
  BSgenome <- validBSgenome(BSgenome)
  
  #First Lets Get Distance to Nearest Gene Start
  message("Annotating Peaks : Nearest Gene")
  distPeaks <- distanceToNearest(peakSummits, GenomicRanges::resize(geneAnnotation$genes, 1, "start"), ignore.strand = TRUE)
  mcols(peaks.condition)$distToGeneStart <- mcols(distPeaks)$distance
  mcols(peaks.condition)$nearestGene <- mcols(geneAnnotation$genes)$symbol[subjectHits(distPeaks)]
  message("Annotating Peaks : Gene")
  promoters <- extendGR(GenomicRanges::resize(geneAnnotation$genes, 1, "start"), upstream = promoterRegion[1], downstream = promoterRegion[2])
  op <- overlapsAny(peakSummits, promoters, ignore.strand = TRUE)
  og <- overlapsAny(peakSummits, geneAnnotation$genes, ignore.strand = TRUE)
  oe <- overlapsAny(peakSummits, geneAnnotation$exons, ignore.strand = TRUE)
  type <- rep("Distal", length(peaks.condition))
  type[which(og & oe)] <- "Exonic"
  type[which(og & !oe)] <- "Intronic"
  type[which(op)] <- "Promoter"
  mcols(peaks.condition)$peakType <- type
  
  distTSS <- distanceToNearest(peakSummits, GenomicRanges::resize(geneAnnotation$TSS, 1, "start"), ignore.strand = TRUE)
  mcols(peaks.condition)$distToTSS <- mcols(distTSS)$distance
  if("symbol" %in% colnames(mcols(geneAnnotation$TSS))){
    mcols(peaks.condition)$nearestTSS <- mcols(geneAnnotation$TSS)$symbol[subjectHits(distTSS)]
  }else if("tx_name" %in% colnames(mcols(geneAnnotation$TSS))){
    mcols(peaks.condition)$nearestTSS <- mcols(geneAnnotation$TSS)$tx_name[subjectHits(distTSS)]
  }
  #Get NucleoTide Content
  message("Annotating Peaks : GC")
  nucFreq <- BSgenome::alphabetFrequency(getSeq(BSgenome, peaks.condition))
  mcols(peaks.condition)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
  mcols(peaks.condition)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
  
  
  ## 3. Add annotated peakset to proj5_signac:
  ############################################
  proj5.VTA.condition <- addPeakSet(ArchRProj =  proj5.VTA.condition, peakSet = peaks.condition, force = TRUE)
  proj5.VTA.condition <- addPeakMatrix(proj5.VTA.condition, force = TRUE, threads = 1)
  proj5.VTA.condition <- addMotifAnnotations(ArchRProj = proj5.VTA.condition, motifSet = "cisbp", name = "Motif", force = TRUE)
  
  #Backup proj5_signac object:
  saveRDS(proj5.VTA.condition, glue("231205.proj5.VTA.{condition}.annotated.rds"))
  saveRDS(peaks.condition, glue("231205.projVTA_{condition}_peaks_metadata.rds"))
}


# ========== 6. Peakset overlap analysis ==========
proj5.VTA <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/final_peaksets/231018_VTAvsSN/231019_ArchR_DARs_and_motifs/231102.proj5_signac_VTA.annotated.rds")

pool_peaks <- data.frame(getPeakSet(proj5.VTA) ) %>% dplyr::select(seqnames, start, end)
write_tsv(pool_peaks, "231205.pool.peaks.bed", col_names = FALSE)

saline_peaks <- data.frame(readRDS("231205.projVTA_saline_peaks_metadata.rds")) %>% dplyr::select(seqnames, start, end)
write_tsv(saline_peaks, "231205.saline.peaks.bed", col_names = FALSE)

h1coc_peaks <- data.frame(readRDS("231205.projVTA_h1_cocaine_peaks_metadata.rds")) %>% dplyr::select(seqnames, start, end)
write_tsv(h1coc_peaks, "231205.h1_cocaine.peaks.bed", col_names = FALSE)

h4coc_peaks <- data.frame(readRDS("231205.projVTA_h4_cocaine_peaks_metadata.rds")) %>% dplyr::select(seqnames, start, end)
write_tsv(h4coc_peaks, "231205.h4_cocaine.peaks.bed", col_names = FALSE)

h8coc_peaks <- data.frame(readRDS("231205.projVTA_h8_cocaine_peaks_metadata.rds")) %>% dplyr::select(seqnames, start, end)
write_tsv(h8coc_peaks, "231205.h8_cocaine.peaks.bed", col_names = FALSE)

h24coc_peaks <- data.frame(readRDS("231205.projVTA_h24_cocaine_peaks_metadata.rds")) %>% dplyr::select(seqnames, start, end)
write_tsv(h24coc_peaks, "231205.h24_cocaine.peaks.bed", col_names = FALSE)

d14coc_peaks <- data.frame(readRDS("231205.projVTA_d14_cocaine_peaks_metadata.rds")) %>% dplyr::select(seqnames, start, end)
write_tsv(d14coc_peaks, "231205.d14_cocaine.peaks.bed", col_names = FALSE)

#Locally: run intervene. 
#GUIX_PROFILE="/fast/home/l/lzeared/profiles/r_multiome230310"
#. "$GUIX_PROFILE/etc/profile"

#In the cluster, getting the "only" peaks
#bedtools intersect -v -a 231205.saline.peaks.bed -b 231205.pool.peaks.bed > saline_only.bed

#bedtools intersect -v -a 231205.h1_cocaine.peaks.bed -b 231205.pool.peaks.bed | wc -l
#[lzeared@max-login6.mdc-berlin.net:/fast/AG_Pombo/luna/2023_vta_multiome/6_231107_updated_analyses/2_ATAC_alternative_processing/2_231123_condition_specific_peaks] $ bedtools intersect -v -a 231205.h4_cocaine.peaks.bed -b 231205.pool.peaks.bed > h4_cocaine_only.bed
#[lzeared@max-login6.mdc-berlin.net:/fast/AG_Pombo/luna/2023_vta_multiome/6_231107_updated_analyses/2_ATAC_alternative_processing/2_231123_condition_specific_peaks] $ bedtools intersect -v -a 231205.h8_cocaine.peaks.bed -b 231205.pool.peaks.bed > h8_cocaine_only.bed
#[lzeared@max-login6.mdc-berlin.net:/fast/AG_Pombo/luna/2023_vta_multiome/6_231107_updated_analyses/2_ATAC_alternative_processing/2_231123_condition_specific_peaks] $ bedtools intersect -v -a 231205.h24_cocaine.peaks.bed -b 231205.pool.peaks.bed > h24_cocaine_only.bed
#[lzeared@max-login6.mdc-berlin.net:/fast/AG_Pombo/luna/2023_vta_multiome/6_231107_updated_analyses/2_ATAC_alternative_processing/2_231123_condition_specific_peaks] $ bedtools intersect -v -a 231205.d14_cocaine.peaks.bed -b 231205.pool.peaks.bed > d14_cocaine_only.bed


#Load "only" peaks and add peakID
sal <- read_tsv("saline_only.bed", col_names = FALSE) %>% dplyr::mutate(peakID = glue("{X1}:{X2}-{X3}"))
h1coc <- read_tsv("h1_cocaine_only.bed", col_names = FALSE) %>% dplyr::mutate(peakID = glue("{X1}:{X2}-{X3}"))
h4coc <- read_tsv("h4_cocaine_only.bed", col_names = FALSE) %>% dplyr::mutate(peakID = glue("{X1}:{X2}-{X3}"))
h8coc <- read_tsv("h8_cocaine_only.bed", col_names = FALSE) %>% dplyr::mutate(peakID = glue("{X1}:{X2}-{X3}"))
h24coc <- read_tsv("h24_cocaine_only.bed", col_names = FALSE) %>% dplyr::mutate(peakID = glue("{X1}:{X2}-{X3}"))
d14coc <- read_tsv("d14_cocaine_only.bed", col_names = FALSE) %>% dplyr::mutate(peakID = glue("{X1}:{X2}-{X3}"))

#Load whole peakset (annotated), add condition ID and filter the "only"

pool_peaks <- data.frame(getPeakSet(proj5.VTA) ) %>%
  dplyr::mutate(peakID = glue("{seqnames}:{start}-{end}"), 
                call_in = "pool") %>% 
  dplyr::select(-idx)

saline_peaks <- data.frame(readRDS("231205.projVTA_saline_peaks_metadata.rds")) %>%
  dplyr::mutate(peakID = glue("{seqnames}:{start}-{end}"), 
                call_in = "saline_only") %>% 
  dplyr::filter(peakID %in% sal$peakID) %>% 
  dplyr::select(colnames(pool_peaks))

h1coc_peaks <- data.frame(readRDS("231205.projVTA_h1_cocaine_peaks_metadata.rds")) %>%
  dplyr::mutate(peakID = glue("{seqnames}:{start}-{end}"), 
                call_in = "h1coc_only") %>% 
  dplyr::filter(peakID %in% h1coc$peakID) %>% 
  dplyr::select(colnames(pool_peaks))

h4coc_peaks <- data.frame(readRDS("231205.projVTA_h4_cocaine_peaks_metadata.rds")) %>%
  dplyr::mutate(peakID = glue("{seqnames}:{start}-{end}"), 
                call_in = "h4coc_only") %>% 
  dplyr::filter(peakID %in% h4coc$peakID) %>% 
  dplyr::select(colnames(pool_peaks))

h8coc_peaks <- data.frame(readRDpool_peaksS("231205.projVTA_h8_cocaine_peaks_metadata.rds")) %>%
  dplyr::mutate(peakID = glue("{seqnames}:{start}-{end}"), 
                call_in = "h8coc_only") %>% 
  dplyr::filter(peakID %in% h8coc$peakID) %>% 
  dplyr::select(colnames(pool_peaks))

h24coc_peaks <- data.frame(readRDS("231205.projVTA_h24_cocaine_peaks_metadata.rds")) %>%
  dplyr::mutate(peakID = glue("{seqnames}:{start}-{end}"), 
                call_in = "h24coc_only") %>% 
  dplyr::filter(peakID %in% h24coc$peakID) %>% 
  dplyr::select(colnames(pool_peaks))

d14coc_peaks <- data.frame(readRDS("231205.projVTA_d14_cocaine_peaks_metadata.rds")) %>%
  dplyr::mutate(peakID = glue("{seqnames}:{start}-{end}"), 
                call_in = "d14coc_only") %>% 
  dplyr::filter(peakID %in% d14coc$peakID) %>% 
  dplyr::select(colnames(pool_peaks))


# ========== 7. Construct extended peak catalogue ==========

extended_peakset <- rbind(pool_peaks, saline_peaks, h1coc_peaks, h4coc_peaks, h8coc_peaks, h24coc_peaks, d14coc_peaks)
extended_peakset.gr <- sort(GRanges(extended_peakset))
extended_peakset.gr$length <- width(extended_peakset.gr)

extended_peakset.gr.collapsed <- GenomicRanges::reduce(extended_peakset.gr)

#Repeat annotations for the extended dataset:

geneAnnotation = getGeneAnnotation(proj5.VTA)
genomeAnnotation = getGenomeAnnotation(proj5.VTA)
promoterRegion = c(2000, 500) #Extended from 100 to 500

#Validate:
peaks_extended <- .validGRanges(extended_peakset.gr.collapsed)
peakSummits <- GenomicRanges::resize(peaks_extended,1,"center")
geneAnnotation$genes <- .validGRanges(geneAnnotation$genes)
geneAnnotation$exons <- .validGRanges(geneAnnotation$exons)
geneAnnotation$TSS <- .validGRanges(geneAnnotation$TSS)
BSgenome <- eval(parse(text = genomeAnnotation$genome))
BSgenome <- validBSgenome(BSgenome)

#Annotate:
peaks_extended <- .validGRanges(peaks_extended)
peakSummits <- GenomicRanges::resize(peaks_extended,1,"center")
geneAnnotation$genes <- .validGRanges(geneAnnotation$genes)
geneAnnotation$exons <- .validGRanges(geneAnnotation$exons)
geneAnnotation$TSS <- .validGRanges(geneAnnotation$TSS)
BSgenome <- validBSgenome(BSgenome)

#First Lets Get Distance to Nearest Gene Start
message("Annotating Peaks : Nearest Gene")
distPeaks <- distanceToNearest(peakSummits, GenomicRanges::resize(geneAnnotation$genes, 1, "start"), ignore.strand = TRUE)
mcols(peaks_extended)$distToGeneStart <- mcols(distPeaks)$distance
mcols(peaks_extended)$nearestGene <- mcols(geneAnnotation$genes)$symbol[subjectHits(distPeaks)]
message("Annotating Peaks : Gene")
promoters <- extendGR(GenomicRanges::resize(geneAnnotation$genes, 1, "start"), upstream = promoterRegion[1], downstream = promoterRegion[2])
op <- overlapsAny(peakSummits, promoters, ignore.strand = TRUE)
og <- overlapsAny(peakSummits, geneAnnotation$genes, ignore.strand = TRUE)
oe <- overlapsAny(peakSummits, geneAnnotation$exons, ignore.strand = TRUE)
type <- rep("Distal", length(peaks_extended))
type[which(og & oe)] <- "Exonic"
type[which(og & !oe)] <- "Intronic"
type[which(op)] <- "Promoter"
mcols(peaks_extended)$peakType <- type

distTSS <- distanceToNearest(peakSummits, GenomicRanges::resize(geneAnnotation$TSS, 1, "start"), ignore.strand = TRUE)
mcols(peaks_extended)$distToTSS <- mcols(distTSS)$distance
if("symbol" %in% colnames(mcols(geneAnnotation$TSS))){
  mcols(peaks_extended)$nearestTSS <- mcols(geneAnnotation$TSS)$symbol[subjectHits(distTSS)]
}else if("tx_name" %in% colnames(mcols(geneAnnotation$TSS))){
  mcols(peaks_extended)$nearestTSS <- mcols(geneAnnotation$TSS)$tx_name[subjectHits(distTSS)]
}
#Get NucleoTide Content
message("Annotating Peaks : GC")
nucFreq <- BSgenome::alphabetFrequency(getSeq(BSgenome, peaks_extended))
mcols(peaks_extended)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
mcols(peaks_extended)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)


#Add to proj
proj6.VTA <- proj5.VTA
proj6.VTA <- addPeakSet(ArchRProj =  proj6.VTA, peakSet = peaks_extended, force = TRUE)
proj6.VTA <- addPeakMatrix(proj6.VTA, force = TRUE, threads = 1)
proj6.VTA <- addMotifAnnotations(ArchRProj = proj6.VTA, motifSet = "cisbp", name = "Motif", force = TRUE)

# Backup proj5_signac object:
# saveRDS(proj6.VTA, glue("231218.proj6.VTA.annotated.rds"))
# saveRDS(extended_peakset.gr.collapsed, glue("231218.proj6.VTA_extended_peaks_metadata.rds"))

#Add call_in information (if overlap several calls: misc)
peaks_extended.df <- data.frame(peaks_extended) %>%
  dplyr::mutate(peakID = glue("{seqnames}:{start}-{end}"))

peaks_extended.df.annot <- left_join(peaks_extended.df, extended_peakset[, 13:14], by = "peakID")
is.na(peaks_extended.df.annot$call_in)

# ========== 8. Peak feature distributions ==========

# Peak type 
df.peaktype <- peaks_extended.df.annot %>% dplyr::select(peakType, call_in) %>%
  group_by(call_in, peakType) %>%
  summarize(count = n())

df.peaktype$call_in <- factor(df.peaktype$call_in, levels = c("pool", "saline_only", "h1coc_only", "h4coc_only", "h8coc_only", "h24coc_only", "d14coc_only"))
plot.peaktype <- ggplot(df.peaktype) +
  geom_bar(aes(y = count, x = call_in, fill = peakType), 
           position="fill", stat="identity")

# Peak length distributions
callin_colors <- c("gray80", "black", "#617641", "#C48208", "#326186", "#AE430A", "#564686")
callin_names <- c("pool", "saline_only", "h1coc_only", "h4coc_only", "h8coc_only", "h24coc_only", "d14coc_only")
names(callin_colors) <- callin_names

peaks_extended.df.annot$call_in <- factor(peaks_extended.df.annot$call_in, levels = callin_names)
plot.length <- ggplot(peaks_extended.df.annot, aes(x=width, color=call_in, fill=call_in)) +
  geom_density(alpha=0.3)+
  #geom_histogram(position="identity", alpha=0.5) +
  scale_color_manual(values=callin_colors) +
  scale_fill_manual(values=callin_colors) +
  theme_minimal() +
  scale_x_log10() +
  facet_wrap(~ call_in)


# ========== 9. Accessibility quantification ==========
# GroupSE per condition
# RPKM normalization
# Accessibility distributions across conditions

proj6.VTA.conditionSE <- getGroupSE(proj6.VTA, groupBy = "condition", divideN = FALSE, useMatrix = "PeakMatrix")
proj6.VTA.conditionSE.df <- data.frame(proj6.VTA.conditionSE@assays@data$PeakMatrix)
peakID <- paste0(proj6.VTA.conditionSE@elementMetadata@listData$seqnames, ":",
                     proj6.VTA.conditionSE@elementMetadata@listData$start, "-",
                     proj6.VTA.conditionSE@elementMetadata@listData$end)
proj6.VTA.conditionSE.df2 <- cbind(peakID, proj6.VTA.conditionSE.df)


#Add annotations
extended_peakset_tojoin <- peaks_extended.df.annot %>% 
  dplyr::mutate(length = end-start) %>% 
  dplyr::select(peakID, length, call_in)

proj6.VTA.conditionSE.df.complete <- left_join(proj6.VTA.conditionSE.df2, extended_peakset_tojoin, by = "peakID")


#RPKMs:
proj6.VTA.conditionSE.df.complete.rpkm <- proj6.VTA.conditionSE.df.complete
for (col in c("d14_cocaine", "h1_cocaine", "h24_cocaine", "h4_cocaine", "h8_cocaine", "saline")) {
  proj6.VTA.conditionSE.df.complete.rpkm[[paste0(col, "_rpkm")]] <- (proj6.VTA.conditionSE.df.complete.rpkm[[col]] / proj6.VTA.conditionSE.df.complete.rpkm$length) * 1e6 / sum(proj6.VTA.conditionSE.df.complete.rpkm[[col]] / proj6.VTA.conditionSE.df.complete.rpkm$length)
}

#write_tsv(proj6.VTA.conditionSE.df.complete.rpkm, "231219_proj6.VTA.conditionSE.df.complete.rpkm.tsv")

#Conditions
df3 <- melt(proj6.VTA.conditionSE.df.complete.rpkm[9:15])
condition.rpkm.length <- ggplot(df3, aes(x=value, color=call_in, fill=call_in)) +
  geom_density(alpha=0.4)+
  #geom_histogram(position="identity", alpha=0.5) +
  scale_color_manual(values=callin_colors) +
  scale_fill_manual(values=callin_colors) +
  theme_minimal() +
  scale_x_log10() +
  facet_wrap(~ variable)

# ========== 10. Quality Control and Filtering ==========
# Low-accessibility peak identification
# Summary statistics
# random_peaks_example <- proj6.VTA.conditionSE.df.complete.rpkm %>%
#   dplyr::filter(call_in != "NA") %>% 
#   dplyr::select(peakID, call_in) %>% 
#   group_by(call_in) %>%
#   sample_n(size = 5, replace = FALSE)

#write_tsv(random_peaks_example, "231205.random_peaks_example.tsv")
#Check peaks with RPKM < 2 in all 
less_2rpkm <- proj6.VTA.conditionSE.df.complete.rpkm %>% dplyr::mutate(less_2rpkm = ifelse(d14_cocaine_rpkm <= 2 & h1_cocaine_rpkm <= 2 & h24_cocaine_rpkm <= 2 & h4_cocaine_rpkm <= 2 & h8_cocaine_rpkm <=2 & saline_rpkm <= 2, "yes", "no"))
less_2rpkm.complete <- left_join(less_2rpkm, extended_peakset, by = "peakID")


# ========== 11. Save Outputs ==========
# save.image("231219.condition_specific_and_RPKMs.rds")

# End of script
