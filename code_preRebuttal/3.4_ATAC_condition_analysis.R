# ===========================================
# Script Title: condition-based ATAC Analysis (VTA-DNs)
# Author: Luna Zea Redondo
# Date: 2023-06-07
# Description:
#   This script contains code to analyze time series ATAC-seq data. It includes:
#     - Marker peak and gene detection
#     - Motif enrichment analysis
#     - Differential peak testing (Wilcoxon)
#     - Volcano and MA plots per contrast
#     - Integration of Seurat UMAPs into ArchR
#     - Fragment QC per condition
#     - Upset plots of differential peaks
#     - Preparation of Pando input data
#
#   Note: These steps represent a first overview and have been further refined in subsequent scripts.

# ===========================================


# ========== 0. Setup ==========
rm(list = ls(all.names = TRUE))
gc()
.libPaths(c("~/profiles/r_multiome230310/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)
`%notin%` <- Negate(`%in%`)

# ========== 1. Load Required Libraries ==========
library(Seurat)
library(Signac)
library(ArchR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(glue)
library(parallel)
library(ComplexHeatmap)
library(Cairo)
library(ggrepel)
library(ggExtra)
library(PupillometryR)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVARmotifs)
library(presto)
library(enrichR)
library(DOSE)
library(RColorBrewer)
library(Pando)
library(scales)

# ========== 2. Define Paths and Metadata ==========
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR"
setwd(dir)

addArchRThreads(threads = 16)
addArchRGenome("mm10")

# Sample Colors
colors_CR1 <- c("#9F2365", "#617641", "#C48208","#326186","#AE430A", "#564686")
names(colors_CR1) <- c("m30_cocaine_R1", "h1_cocaine_R1", "h4_cocaine_R1", "h8_cocaine_R1", "h24_cocaine_R1", "d14_cocaine_R1")

colors_CR2 <- c("#6C8448", "#D78F09", "#376C95", "#C14B0B", "#6754A0")
names(colors_CR2) <- c("h1_cocaine_R2", "h4_cocaine_R2", "h8_cocaine_R2", "h24_cocaine_R2", "d14_cocaine_R2")

colors_CR3 <- c("#D4520C", "#7D6CB2")
names(colors_CR3) <- c("h24_cocaine_R3", "d14_cocaine_R3")

colors_SR1 <- c("#B2C596", "#F9CB76", "#88B2D3", "#F7A578", "#A194C7")
names(colors_SR1) <- c("h1_saline_R1", "h4_saline_R1", "h8_saline_R1", "h24_saline_R1", "d14_saline_R1")

sample_colors <- c(colors_SR1, colors_CR1, colors_CR2, colors_CR3)
samples <- names(sample_colors)

condition_names <- c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")
condition_colors <- c("black", "#9F2365", "#617641", "#C48208", "#326186", "#AE430A", "#564686")
names(condition_colors) <- condition_names

all_contrasts <- c(
  "h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline", 
  "h24_cocaine-saline", "d14_cocaine-saline",
  "h4_cocaine-h1_cocaine", "h8_cocaine-h4_cocaine",
  "h24_cocaine-h8_cocaine", "d14_cocaine-h24_cocaine"
)

# ========== 3. Load Preprocessed Objects ==========
seurat_dataframe_complete <- readRDS("seurat_dataframe_complete.rds")
RNAreduction <- readRDS("RNAreduction.rds")
ATACreduction <- readRDS("ATACreduction.rds")
WNNreduction <- readRDS("WNNreduction.rds")
wnn.seu <- readRDS("wnn.seu.rds")
proj4 <- readRDS("peakCalling/230602_ATAC_proj4.condition.rds")

# ========== 4. Integrate Seurat Reductions into ArchR ==========
proj4 <- proj4[proj4$cellNames %in% rownames(seurat_dataframe_complete),]
proj4_seurat <- proj4
proj4_seurat$wsnn_res.0.4 <- wnn.seu$wsnn_res.0.4
proj4_seurat@cellColData$wsnn_res.0.4 <- as.character(proj4_seurat@cellColData$wsnn_res.0.4)

proj4_seurat@embeddings$seuratUMAP_RNA <- SimpleList(df = RNAreduction, params = list())
proj4_seurat@embeddings$seuratUMAP_ATAC <- SimpleList(df = ATACreduction, params = list())
proj4_seurat@embeddings$seuratUMAP_WNN <- SimpleList(df = WNNreduction, params = list())

# ========== 5. QC Metrics per Sample ==========
ATACmetrics <- as.data.frame(proj4_seurat@cellColData)
ATACmetrics$Sample <- factor(ATACmetrics$Sample, levels = samples)

ATACmetrics.TSS <- ggplot(ATACmetrics, aes(x=Sample, y=TSSEnrichment, fill = Sample, colour = Sample))+
  geom_flat_violin(position = position_nudge(x = .25), adjust = 2, trim = TRUE) +
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(aes(x = as.numeric(Sample)+0.25, y = TSSEnrichment), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_y_log10() +
  labs(title = "DNs: TSS Enrichment") +
  scale_fill_manual(values = sample_colors) +
  scale_color_manual(values = sample_colors) +
  theme_minimal() + theme(legend.position = "none")

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

ATACmetrics.TSS/ATACmetrics.nFrags


# ========== 6. Marker gene and peak detection ==========

# 6.1. Marker genes (based on ATAC data)
proj4_seurat@cellColData$wsnn_res.0.4 <- as.character(proj4_seurat@cellColData$wsnn_res.0.4)
#proj4_seurat@cellColData$wnnCluster <- as.character(proj4_seurat@cellColData$wnnCluster)

#Gene markers:
markersGE <- getMarkerFeatures(
  ArchRProj = proj4_seurat, 
  useMatrix = "GeneExpressionMatrix", 
  bias = c("TSSEnrichment", "log10(nFrags)"),
  groupBy = "condition",
  testMethod = "wilcoxon", 
  threads=1
)
markersGE <- markersGE[, c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")]
markerList <- getMarkers(markersGE, cutOff = "FDR <= 0.3 & Log2FC >= 0.5")
all_gene_markers <- as.data.frame(as.data.frame(markerList))

markerList.top5 <- as.data.frame(markerList) %>%
  group_by(group_name) %>% 
  slice_min(FDR, n= 5, with_ties = FALSE) %>% 
  ungroup() %>% 
  select(name) 

heatmapGE <- plotMarkerHeatmap(
  seMarker = markersGE, 
  cutOff = "FDR <= 0.3 & Log2FC >= 0.5", 
  transpose = FALSE, 
  clusterCols = FALSE, 
  labelRows = FALSE
)

heatmapGE_CH <- ComplexHeatmap::draw(heatmapGE, heatmap_legend_side = "bot", annotation_legend_side = "bot")
png(file="markerGenes.heatmapGE.png", units = "px", width=500,height=1000)
draw(heatmapGE_CH)
dev.off()

# Which functions? EnrichR analysis:
dbs <- listEnrichrDbs()
dbs <- c("TF_Perturbations_Followed_by_Expression",
         "BioPlanet_2019", 
         "GO_Molecular_Function_2021")

enrichR_mastertable <-data.frame()
conditions <-   c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")

for (j in conditions) {
  print(glue("Calculating enrichR results for {j}"))
  genes_toenrich <- all_gene_markers %>% filter(group_name == j) %>% pull(name) %>%  as.character()
  enriched <- enrichr(genes_toenrich, dbs)
  enriched.df <- rbindlist(enriched, idcol = TRUE)
  enriched.df <- enriched.df %>% mutate(groupID = j)
  enrichR_mastertable <- rbind(enrichR_mastertable, enriched.df)
}
#write_tsv(enrichR_mastertable, glue("221111_proj4_markerGS_enrichR_mastertable.tsv"))

enrichR_mastertable_toplot <- enrichR_mastertable %>%
  group_by(groupID, .id) %>% 
  slice_min(Adjusted.P.value, n=2, with_ties = FALSE) %>% 
  mutate(Overlap = parse_ratio(Overlap))

enrichR_mastertable_toplot$groupID <- factor(enrichR_mastertable_toplot$groupID, levels = conditions) 
enrichR_mastertable_toplot <- enrichR_mastertable_toplot %>% arrange(.id, groupID)
enrichR_mastertable_toplot$Term <- factor(enrichR_mastertable_toplot$Term, levels = unique(enrichR_mastertable_toplot$Term))

enrichRplot <- ggplot(enrichR_mastertable_toplot,
                      aes(x = groupID, y = fct_rev(Term))) + 
  geom_point(aes(size = Overlap, color = Adjusted.P.value)) +
  theme_classic() +
  theme_bw(base_size = 10) +
  ylab(NULL) +
  xlab(NULL) +
  facet_wrap(~ .id, scales = "free", ncol = 1) +
  scale_colour_gradient2(limits=c(0, 0.1), low="#cb181d", mid= "#fb6a4a", high= "#fcbba1",
                         na.value = "grey50",
                         midpoint = 0.05)
enrichRplot
#ggsave(filename = glue("markerGenes.enrichRplot.png"), enrichRplot, units = "px", device = "png",width=3500,height=2000,dpi = 300)


# 6.2. Marker peaks

proj4_seurat_2reps <- proj4_seurat[proj4_seurat$condition != "m30_cocaine"]
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj4_seurat_2reps,
  useMatrix = "PeakMatrix",
  groupBy = "condition",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  threads=1
)

#markersPeaks <- markersPeaks[, c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")]
markersPeaks <- markersPeaks[, c("saline", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")]

markerList.peaks <- getMarkers(markersPeaks, cutOff = "Pval <= 0.05 & Log2FC >= 0.5")
markerList.df <- as.data.frame(markerList.peaks)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "Pval <= 0.05 & Log2FC >= 0.5",
  transpose = FALSE, 
  clusterCols = TRUE
)
#ggsave(filename = glue("markerPeaks_all_heatmapPeaks.png"), heatmapPeaks, units = "px", device = "png",width=3500,height=2000,dpi = 300)

#Add annotation
proj4_seurat <- addMotifAnnotations(ArchRProj = proj4_seurat, motifSet = "cisbp", name = "Motif")

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj4_seurat,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 1"
)

enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, cutOff = 1, transpose = TRUE)

proj4_seurat$condition_simple <- ifelse(proj4_seurat$Sample %in% c("h1_saline_R1", "h4_saline_R1", "h8_saline_R1"), "saline_early", 
                                        ifelse(proj4_seurat$Sample %in% c("h24_saline_R1", "d14_saline_R1"), "saline_late", proj4_seurat$condition))

proj4_seurat <- addMotifAnnotations(ArchRProj = proj4_seurat, motifSet = "cisbp", name = "Motif", force = TRUE)
proj4_seurat <- addDeviationsMatrix(
  ArchRProj = proj4_seurat,
  peakAnnotation = "Motif",
  force = TRUE,
  threads=1
)

# ========== 7. Pairwise Differential Peak Testing ==========
# Loop over all contrasts to compute:
# - Differential peaks
# - Volcano and MA plots
# - Motif enrichment (peakAnnoEnrichment)

condition_colors <- c("black", "#9F2365", "#617641", "#C48208", "#326186", "#AE430A", "#564686")
condition_names <- c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")
names(condition_colors) <- condition_names
all_contrasts <- c("h1_cocaine-saline",
                   "h4_cocaine-saline",
                   "h8_cocaine-saline",
                   "h24_cocaine-saline",
                   "d14_cocaine-saline", 
                   "h4_cocaine-h1_cocaine",
                   "h8_cocaine-h4_cocaine", 
                   "h24_cocaine-h8_cocaine", 
                   "d14_cocaine-h24_cocaine")

#dir.create("peaks_and_motifs2")
DiffPeaks_complete_results <-data.frame()
EnrichedMotifs_complete_results <- data.frame()

for (contrast in all_contrasts) {
  print(glue("Analyzing: {contrast}"))
  group1=str_split(contrast, "-")[[1]][1]
  group2=str_split(contrast, "-")[[1]][2]
  contrast_label=contrast
  
  # 1) Calculate diff peaks:
  markerTest <- getMarkerFeatures(
    ArchRProj = proj4_seurat, 
    useMatrix = "PeakMatrix",
    groupBy = "condition",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = group1,
    bgdGroups = group2,
  )
  markerList <- getMarkers(markerTest, cutOff = "FDR <=1")
  
  markerList.df <- as.data.frame(markerTest@elementMetadata) %>% mutate(
    peakID = glue("{seqnames}:{start}-{end}"),
    contrast=contrast_label,
    useGroups = group1,
    bdgGroups = group2,
    Log2FC = markerTest@assays@data@listData[["Log2FC"]]$x,
    Pval = markerTest@assays@data@listData[["Pval"]]$x,
    FDR = markerTest@assays@data@listData[["FDR"]]$x, 
    Mean = markerTest@assays@data@listData[["Mean"]]$x, 
    MeanBGD = markerTest@assays@data@listData[["MeanBGD"]]$x, 
    LM = log2((Mean + MeanBGD)/2 + 1),
    Differential = ifelse(Pval <= 0.05 & Log2FC >= 0.5, "Upregulated (pval <= 0.05 and Log2FC >= 0.5)",
                          ifelse(Pval <= 0.05 & Log2FC <= -0.5, "Downregulated (pval <= 0.05 and Log2FC <= -0.5)", "No significant")))
  
  DiffPeaks_complete_results <- rbind(DiffPeaks_complete_results, markerList.df)
  # 2)Plot volcano
  #Colors
  volcano_colors <- c(condition_colors[group1], condition_colors[group2], "gray")
  names(volcano_colors) <- c("Upregulated (pval <= 0.05 and Log2FC >= 0.5)", "Downregulated (pval <= 0.05 and Log2FC <= -0.5)", "No significant")
  
  #Plot
  volcano_plot <- ggplot(markerList.df, aes(x = Log2FC, y = -log10(Pval))) +
    geom_point(aes(color = Differential)) +
    scale_color_manual(values = volcano_colors) +
    theme_classic(base_size = 12) + theme(legend.position = "bottom") +
    geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") + 
    geom_vline(xintercept = -0.5, colour="#990000", linetype="dashed") + 
    geom_vline(xintercept = 0.5, colour="#990000", linetype="dashed") +
    labs(title = glue("Differential Expression Analysis: {group1} vs {group2}"),
         subtitle = "pval <0.05; |logFC| >= 0.5", 
         x= "Log2 (FC)", 
         y="-Log10 (p value)") +
    guides(fill=guide_legend(title="Differential")) 
  
  volcano_plot
  
  #Plot MA
  MA_plot <- ggplot(markerList.df, aes(x = LM, y = Log2FC)) +
    geom_point(alpha=0.6, aes(color = Differential)) +
    geom_point(alpha=0.6,aes(color = Differential), data=markerList.df[markerList.df$Differential != "No significant",]) +
    scale_color_manual(values = volcano_colors) +
    theme_classic(base_size = 12) + theme(legend.position = "bottom") +
    labs(title = glue("Differential Expression Analysis: {group1} vs {group2}"),
         subtitle = "pval <0.05; |logFC| >= 0.5", 
         x= "Log2 Mean", 
         y="Log2 Fold Change") +
    guides(fill=guide_legend(title="Differential")) 
  
  
  # 3) Motifs
  # a) Motifs up:
  motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj4_seurat,
    peakAnnotation = "Motif",
    cutOff = glue("Pval <= 0.05 & Log2FC >= 0.5"))
  
  # dataframe
  dfUp <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1]) %>% mutate(
    contrast=contrast_label, 
    group = group1
  )
  dfUp <- dfUp[order(dfUp$mlog10Padj, decreasing = TRUE),]
  dfUp$rank <- seq_len(nrow(dfUp))
  #plot
  ggUp <- ggplot(dfUp, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point(size = 2) +
    ggrepel::geom_label_repel(
      data = dfUp[rev(seq_len(10)), ], aes(x = rank, y = mlog10Padj, label = TF), 
      size = 3,
      nudge_x = 2,
      color = "black", 
      max.overlapsload("230626.ATAC_DN_condition_backup.rds") = 10
    ) + theme_ArchR() + 
    labs(title = glue("TF-Motifs enriched in {group1}")) +
    theme(plot.title = element_text(size = 15)) +
    ylab("-log10(P-adj) Motif Enrichment") + 
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 
  
  ggUp
  
  
  # b) Motifs down:
  motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj4_seurat,
    peakAnnotation = "Motif",
    cutOff = glue("Pval <= 0.05 & Log2FC <= -0.5"))
  
  # dataframe
  dfDo <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1]) %>% mutate(
    contrast=contrast_label, 
    group = group2
  )
  dfDo <- dfDo[order(dfDo$mlog10Padj, decreasing = TRUE),]
  dfDo$rank <- seq_len(nrow(dfDo))
  #plot
  ggDo <- ggplot(dfDo, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point(size = 2) +
    ggrepel::geom_label_repel(
      data = dfDo[rev(seq_len(10)), ], aes(x = rank, y = mlog10Padj, label = TF), 
      size = 3,
      nudge_x = 2,
      color = "black", 
      max.overlaps = 10
    ) + theme_ArchR() + 
    labs(title = glue("TF-Motifs enriched in {group2}")) +
    theme(plot.title = element_text(size = 15)) +
    ylab("-log10(P-adj) Motif Enrichment") + 
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 
  
  ggDo
  EnrichedMotifs_complete_results <- rbind(EnrichedMotifs_complete_results, dfUp, dfDo)
  
  contrast_results <- volcano_plot + ggUp + ggDo
  contrast_results2 <- MA_plot + ggUp + ggDo
  
  ggsave(filename = glue("peaks_and_motifs2/{contrast_label}_peakANDmotifs_results.png"), contrast_results, units = "px", device = "png", width=6000,height=2000,dpi = 300)
  ggsave(filename = glue("peaks_and_motifs2/{contrast_label}_peakANDmotifs_results_MAplot.png"), contrast_results2, units = "px", device = "png", width=6000,height=2000,dpi = 300)
  
}
#save tables
# write_tsv(DiffPeaks_complete_results, glue("peaks_and_motifs2/230607_DiffPeaks_complete_results.tsv"))
# write_tsv(EnrichedMotifs_complete_results, glue("peaks_and_motifs2/230607_EnrichedMotifs_complete_results.tsv"))

# ========== 8. Visualization; upset plots (overlap summary) ==========

# ATAC genomic tracks per condition:
getGroupBW(proj4_seurat, groupBy = "condition", threads = 1)


# Peaks overview: genomic distribution and DARs

peaks_differential <- DiffPeaks_complete_results %>%
  filter(Differential != "No significant") %>% 
  select(peakID) %>%
  distinct()

peakset.df <- data.frame(getPeakSet(proj4_seurat)) %>% mutate(
  peakID = glue("{seqnames}:{start}-{end}"), 
  differential = ifelse(peakID %in% peaks_differential$peakID, "yes", "no")) %>% 
  group_by(peakType, differential) %>%
  dplyr::summarise(n = n()) %>% 
  mutate(freq =(n / sum(n)) *100)

differetial_class <- c("no", "yes")
diff_colors <- c("black", "#961D4E")
peakset.gg <- ggplot(data=peakset.df, aes(x=peakType, y=freq, fill = differential)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=round(freq, digits = 1)), vjust=1.6, color="white", size=5)+
  scale_fill_manual(values = diff_colors) + 
  scale_color_manual(values = diff_colors) +
  theme_minimal() +
  guides(fill=guide_legend(title="Differential \n peaks")) +
  labs(title = "Genomic distribution of accesible regions", y = "Percentage (%) of peaks", x = "") +
  theme(title = element_text(color="black", 
                             size=14, angle=0),
        axis.title.x =  element_text(color="black", 
                                     size=14, angle=0),
        axis.title.y =  element_text(color="black", 
                                     size=14, angle=90),
        axis.text.x = element_text(color="black", 
                                   size=12, angle=0),
        axis.text.y = element_text(color="black", 
                                   size=12, angle=0), 
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))


# Heatmap peaks (230607)
############################
#Get average accessibility per peak per sample
peaks_SE <- getGroupSE(proj4_seurat, 
                       useMatrix = "PeakMatrix",
                       groupBy = "condition",
                       logFile = createLogFile("getGroupSE"))

#Add peakID
peaks.df <-as.data.frame(peaks_SE@assays@data@listData$PeakMatrix) 
rownames_peaks <- as.data.frame(peaks_SE@elementMetadata@listData) %>% 
  mutate(peakID = glue("{seqnames}:{start}-{end}"))
peaks.df <- cbind(peaks.df, "peakID" = rownames_peaks$peakID)

#Combine Diff peaks with peak info
differential_results <- left_join(DiffPeaks_complete_results, peaks.df, by = "peakID")

#Plot ATAC:
saline_contrasts <-c("h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline", "h24_cocaine-saline", "d14_cocaine-saline")
time_contrasts <- c("h1_cocaine-saline", "h4_cocaine-h1_cocaine", "h8_cocaine-h4_cocaine", "h24_cocaine-h8_cocaine", "d14_cocaine-h24_cocaine")

heatmap_peaks <- differential_results %>%
  filter(contrast %in% saline_contrasts,
         Differential != "No significant")

#Extract DEGs from the big table and create metadata
duplicated_peaks <- differential_results %>% 
  filter(Differential != "No significant", 
         contrast %in% saline_contrasts)
duplicated_peaks <- duplicated_peaks[duplicated(duplicated_peaks$peakID), "peakID"]

peak_withIDs <- differential_results %>% 
  filter(Differential != "No significant", 
         contrast %in% saline_contrasts) %>% distinct(peakID, .keep_all = TRUE) %>% 
  mutate(groupID = ifelse(peakID %in% duplicated_peaks, "No unique", 
                          ifelse(Log2FC < 0, bdgGroups, useGroups))) %>% 
  arrange(factor(groupID, levels = c(condition_names, "No unique")))

#Prepare matrix for heatmap (in the order I want to be plotted)
mat <- as.matrix(peak_withIDs[, c("saline", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")])
mat_scaled = t(scale(t(mat)))

# clusters <- hclust(dist(mat_scaled), method = 'average')
# plot(clusters)

#samples.coldata
condition <- c("saline", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")
condition <- factor(condition, levels = c("saline", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine"))

cluster <- factor(peak_withIDs$groupID, levels = c(condition)) #heatmap cluster
heatmap_color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)  

ha = HeatmapAnnotation(Condition = condition, show_legend = FALSE, 
                       col = list(Condition = condition_colors), show_annotation_name = FALSE)

row_ha = rowAnnotation(Cluster = cluster, show_legend = FALSE, 
                       col = list(Cluster = c(condition_colors, "No unique" = "gray")), show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"), labels = NULL)


htmp <- Heatmap(mat_scaled, cluster_columns = FALSE, cluster_rows = FALSE, name = "Z-score", column_split =condition, row_split =cluster, row_title = NULL, show_row_names = FALSE,
                show_column_names = FALSE, top_annotation = ha, right_annotation = row_ha, col = heatmap_color, heatmap_legend_param = list(direction = "horizontal"))
htmp_save <- draw(htmp, heatmap_legend_side = "bottom")

# dev.off()
# jpeg(glue("peaks_and_motifs2/DARs_vs_saline.jpeg"), bg = "transparent", width=600, height=700, units = "px")
# print(htmp_save)
# dev.off()


# UpsetPlot (230701)

#Make a list of diff peaks: contrast and direction 
upset.df <- heatmap_peaks %>%  mutate(
  direction = substr(Differential,1,nchar(Differential)-42),
  contrast_direction =  glue("{contrast}_{direction}")
)

upset.df$contrast_direction <- factor(upset.df$contrast_direction, levels = c("h1_cocaine-saline_Up",  "h1_cocaine-saline_Downr", "h4_cocaine-saline_Up",  "h4_cocaine-saline_Downr", 
                                                                              "h8_cocaine-saline_Up",  "h8_cocaine-saline_Downr", "h24_cocaine-saline_Up",  "h24_cocaine-saline_Downr", 
                                                                              "d14_cocaine-saline_Up",  "d14_cocaine-saline_Downr"))
#All combinations
peaks_list <- split(upset.df$peakID, upset.df$contrast_direction)
m1 = make_comb_mat(peaks_list)

comb_size(m1)
UpSet(m1)

#More than 100 elements (combs)
m2 = m1[comb_size(m1) > 100]
ht <- UpSet(m2, set_order = as.integer(1:10),  comb_order = order(comb_size(m2), decreasing = TRUE))

# Extract combinations:
m2_sorted <- sort(comb_size(m2), decreasing = TRUE) 


# ========== 9. Motif Matrix ==========

# Extract motif matrix
peakset.complete <- data.frame(getPeakSet(proj4_seurat)) %>% mutate(
  peakID = glue("{seqnames}:{start}-{end}"), 
  differential = ifelse(peakID %in% peaks_differential$peakID, "yes", "no")) 

motifMatrix <- getMatches(proj4_seurat, "Motif")
motifMatrix <- motifMatrix@assays@data$matches  
rownames(motifMatrix) <- peakset.complete$peakID


df <- as.data.frame(as.matrix(motifMatrix)) # %>% rownames_to_column("peakID")
df2 <- t(df)
motif.complete <- reshape2::melt(df2) %>% filter(value != "FALSE") %>% select(Var1, Var2)
colnames(motif.complete) <- c("TF", "peakID")

motif.final <- motif.complete %>% 
  group_by(peakID) %>% 
  dplyr::summarise(TF = paste(unique(TF), collapse = ', ')) 


# Join motif matrix to peakset
 
peakset.complete$peakID <- factor(peakset.complete$peakID)
peakset.complete2 <- left_join(peakset.complete, motif.final, by = "peakID")

write_tsv(peakset.complete2, "230712.condition_peakSet_with_motifs.txt")
#save.image("230607.condition_ATACanalysis_backup.rds")


# ========== 10. Join archR peak matrix to wnn ==========

wnn.seu <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR/wnn.seu.rds")

#Prepare the peak assay to add to wnn. 
archR_peaks <- getMatrixFromProject(proj4_seurat, useMatrix = "PeakMatrix")
archR_peakMatrix <- archR_peaks@assays@data$PeakMatrix
rownames(archR_peakMatrix) <- peakset.complete$peakID

#Process using signac code. 
data_dir <- "/data/pombo/Luna/MultiomeCocaineTreatments/data_dir/"

samples <- c("m30_cocaine_R1", "h1_saline_R1", "h1_cocaine_R1", "h1_cocaine_R2",
             "h4_saline_R1", "h4_cocaine_R1", "h4_cocaine_R2",
             "h8_saline_R1", "h8_cocaine_R1", "h8_cocaine_R2",
             "h24_saline_R1", "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3",
             "d14_saline_R1", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3")

fragments_dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR/modified_fragments"

for (sample in samples[15:18]) {
  # Analysis of sample
  print(glue("Processing fragments of sample: {sample}"))
  # Step 1: Set the file paths
  input_file <- glue("{data_dir}{sample}.atac_fragments.tsv.gz")
  temp_file <- glue("{fragments_dir}/{sample}.temp.tsv")
  output_file <- glue("{fragments_dir}/{sample}.atac_fragments.tsv")
  
  # Step 2: Extract atac fragments header and saved as temporary file:
  header_lines <- as.data.frame(readLines(gzfile(input_file), n = 52))
  write.table(header_lines, "temp.header.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Step 3: Read and modify the fragment file to keep only DN cells
  input_data <- fread(input_file, header = FALSE)  %>% 
    mutate(V4 = glue("{sample}#{V4}")) %>% 
    filter(V4 %in% colnames(wnn.seu))
  
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


fragment_list <- list()
for (sample in samples) {
  frag.sample.path = glue("{fragments_dir}/{sample}.atac_fragments.tsv.gz")
  sample.cells <- atac_counts[, grepl(sample, colnames(atac_counts))]
  frag.sample <- CreateFragmentObject(path = frag.sample.path, cells = colnames(sample.cells), verbose = FALSE, tolerance = 0.5, validate.fragments = FALSE)
  fragment_list <- append(fragment_list, frag.sample)
}
names(fragment_list) <- samples  
atac_counts <- archR_peakMatrix


# ========== 11. Save output ==========
#save.image("230626.ATAC_DN_condition_backup.rds")









