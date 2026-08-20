# =========================================================
# PROJECT UTILITIES & CONFIGURATION
# =========================================================

# 1. Custom Operators ----
`%notin%` <- Negate(`%in%`)

# 2. Sample Color Palettes ----

# Cocaine R1
colors_CR1 <- c("#9F2365", "#617641", "#C48208","#326186","#AE430A", "#564686")
names(colors_CR1) <- c("m30_cocaine_R1", "h1_cocaine_R1", "h4_cocaine_R1", "h8_cocaine_R1", "h24_cocaine_R1", "d14_cocaine_R1")

# Cocaine R2
colors_CR2 <- c("#6C8448", "#D78F09", "#376C95", "#C14B0B", "#6754A0")
names(colors_CR2) <- c("h1_cocaine_R2", "h4_cocaine_R2", "h8_cocaine_R2", "h24_cocaine_R2", "d14_cocaine_R2")

# Cocaine R3
colors_CR3 <- c("#D4520C", "#7D6CB2")
names(colors_CR3) <- c("h24_cocaine_R3", "d14_cocaine_R3")

# Saline R1
colors_SR1 <- c("#B2C596", "#F9CB76", "#88B2D3", "#F7A578", "#A194C7")
names(colors_SR1) <- c("h1_saline_R1", "h4_saline_R1", "h8_saline_R1", "h24_saline_R1", "d14_saline_R1")

# Combined Global Palette
sample_colors <- c(colors_SR1, colors_CR1, colors_CR2, colors_CR3)


# Simplified class
simpleIdent_colors <- c("gray", "#9F2365", "#617641", "#C48208","#326186","#AE430A", "#564686")
names(simpleIdent_colors) <- c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")

#Timepoints
timepoints <- c("m30", "h1", "h4", "h8", "h24", "d14")
timepoints_colors <- colors_CR1
names(timepoints_colors) <- timepoints


sample_order <- c(
  "m30_cocaine_R1",
  "h1_saline_R1",
  "h1_cocaine_R1", "h1_cocaine_R2",
  "h4_saline_R1",
  "h4_cocaine_R1", "h4_cocaine_R2",
  "h8_saline_R1",
  "h8_cocaine_R1", "h8_cocaine_R2",
  "h24_saline_R1",
  "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3",
  "d14_saline_R1",
  "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3"
)

batches <- read_csv("/fast/AG_Pombo/luna/2026_rebuttal/experimental_batches.csv")


#Color palettes
library(RColorBrewer)
SamplePrep <- unique(batches$sample_prep)
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
SamplePrep_colors <- getPalette(length(SamplePrep))
names(SamplePrep) <- SamplePrep_colors

#Batch 2: FACS
FACSmachine <- unique(batches$FACS)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
FACS_colors <- getPalette(length(FACSmachine))
names(FACS_colors) <- FACSmachine

#Batch 3: GEXlibrary
GEXlibrary <- unique(batches$GEXlibrary)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
GEXlibrary_colors <- getPalette(length(GEXlibrary))
names(GEXlibrary_colors) <- GEXlibrary

#Batch 4: Sequencing
Sequencing <- unique(batches$Sequencing)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
sequencing_colors <- getPalette(length(Sequencing))
names(sequencing_colors) <- Sequencing



all_contrasts <- c("h1_cocaine-saline",
                   "h4_cocaine-saline",
                   "h8_cocaine-saline",
                   "h24_cocaine-saline",
                   "d14_cocaine-saline", 
                   "h4_cocaine-h1_cocaine",
                   "h8_cocaine-h4_cocaine", 
                   "h24_cocaine-h8_cocaine", 
                   "d14_cocaine-h24_cocaine")





annotate_seurat_peaks_archr_style <- function(
    seurat_obj,
    assay = "ATAC",
    geneAnnotation,
    genomeAnnotation,
    promoterRegion = c(1000, 500)
) {
  DefaultAssay(seurat_obj) <- assay
  
  peaks <- granges(seurat_obj[[assay]])
  names(peaks) <- rownames(seurat_obj[[assay]])
  
  # Keep all feature names for later
  all_peak_names <- names(peaks)
  
  # ArchR-style validation
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
  
  # Create full metadata table for all original Seurat peaks.
  # This protects you if keepStandardChromosomes removed chrM, scaffolds, etc.
  full_metadata <- data.frame(
    peakType = NA_character_,
    distToGeneStart = NA_integer_,
    nearestGene = NA_character_,
    distToTSS = NA_integer_,
    nearestTSS = NA_character_,
    GC = NA_real_,
    N = NA_real_,
    row.names = all_peak_names
  )
  
  full_metadata[rownames(peak_metadata), colnames(peak_metadata)] <- peak_metadata
  
  seurat_obj[[assay]] <- AddMetaData(
    object = seurat_obj[[assay]],
    metadata = full_metadata
  )
  
  return(seurat_obj)
}
