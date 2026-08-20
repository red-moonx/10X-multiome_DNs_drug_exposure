# ===========================================
# Script Title: DN-ATAC Core Analysis (combining ArchR + Signac)
# Author: Luna Zea Redondo
# Date: 2023-06-02

# Description:
#   This script processes dopaminergic neuron (DN) ATAC-seq data using a combination 
#   of ArchR and Signac. It performs the following key steps:
#
#   (1) Quality control of all cells and DN-only cells. 
#   (2) Filtering of low-quality DN cells
#   (3) Dimensionality reduction using Iterative LSI, clustering and UMAP visualization.
#   (4) Integration of ATAC (from ArchR) and RNA (from Seurat) modalities into a 
#       unified Seurat WNN object (`wnn.seu`). PCA and LSI are used for RNA/ATAC embeddings,
#       followed by WNN graph construction, clustering, and visualization.
#   (5) Scoring and visualization of SN-specific gene signatures 
#   (6) Analysis of ATAC-specific metrics (TSS enrichment, fragment count) across 
#       WNN clusters, followed by export of UMAP embeddings and metadata.
#   (7) The script prepares the final integrated data for downstream use and saves all reductions and key objects.
#
#   NOTE: Peak calling is run separately in:
#         230602_peakCalling_both_ways.R
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

# ========== Set Working Directory ==========
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR"
setwd(dir)

# ========== ArchR Setup ==========
addArchRThreads(threads = 16)
addArchRGenome("mm10")


# ========== 0) Load Data ==========
# Load ArchR project with all cells
# Subset to DNs using matching barcodes from RNA and store as proj_DNs
load("230601.projb2.rds")
DNs.RNA.seu <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/230421.DN.conos.sub.rds")


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
  select(Sample, log10.nFrags., TSSEnrichment) %>% 
  dplyr::rename(log10nFrags = log10.nFrags.)

allplots <- list()
for (sample in samples) {
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




# ========== 2) ATAC QC — DN-only ==========
# QC exploration for DNs, filtering (TSS ≥ 4, nFrags ≥ 1000)
proj <- proj_DNs
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
#ggsave(filename = "atac_qc_onlyDNs.png", atac_qc_all_plot, units = "px", device = "png", width=6000,height=2000,dpi = 300)

# #Save
dev.size()
width_original = 6.760417
height_original= 5.083333

pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
plot_name <- "chapter3_DNs_atac_allQC"
file_name <- glue("{pdf_dir}/{plot_name}.pdf" )

pdf(file_name, width = width_original, height =height_original)
atac_qc_all_plot
dev.off()

df2 <- as.data.frame(df) %>% rownames_to_column("cellNames") %>% 
  separate(cellNames, c("Sample", "CB"), sep = "#") %>% 
  select(Sample, log10.nFrags., TSSEnrichment) %>% 
  dplyr::rename(log10nFrags = log10.nFrags.)

allplots <- list()
for (sample in samples) {
  df3 <- df2 %>% filter(Sample == sample)
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

table_passing <- df2 %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::summarise(total = n(), 
                   passingQC = sum(log10nFrags > 3 & TSSEnrichment > 4), 
                   percentage = round((passingQC/total)*100))

#View(table_passing)

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
DNs_fromRNA <- DNs.RNA.seu@meta.data %>% dplyr::select(orig.ident, simpleIdent, seurat_clusters, sample_prep, FACS, GEXlibrary, Sequencing) %>% rownames_to_column("cellNames")
proj.metadata <- left_join(proj.metadata,DNs_fromRNA, by = "cellNames") 


# ========== 3) Compare QC metrics after filtering (only DNs) ==========
ATACmetrics <- as.data.frame(proj@cellColData)
ATACmetrics$Sample <- factor(ATACmetrics$Sample, levels = samples)
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

DNs_fromRNA <-  DNs.RNA.seu@meta.data %>%
  tibble::rownames_to_column("cellNames") %>% 
  dplyr::select(cellNames, orig.ident, simpleIdent, seurat_clusters, sample_prep, FACS, GEXlibrary, Sequencing) %>%
  dplyr::filter(cellNames %in% proj2$cellNames) %>%
  tibble::column_to_rownames("cellNames") 

proj <- addCellColData(ArchRProj = proj,
                       data = as.character(DNs_fromRNA$seurat_clusters),
                       cells = as.character(rownames(DNs_fromRNA)),
                       name = "seurat_clusters")


# ========== 4) Dimensionality Reduction & visual inspection ==========
# Add IterativeLSI (default method), Clusters, and UMAP to filtered DN ArchRProj
# Plot UMAPs by Sample, Cluster, Seurat cluster, TSS, nFrags, etc.
proj2 <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

proj2 <- addClusters(
  input = proj2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

proj2 <- addUMAP(
  ArchRProj = proj2,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
)
hex <- hue_pal()(6)
ggplot_pal <- c("0" = hex[1], "1"= hex[2], "2"= hex[3], "3"= hex[4], "4"= hex[5], "5"= hex[6])
p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP", pal = sample_colors, size = 1.4, labelMeans = FALSE) + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) + no_legend() 
p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", size = 1.4,   labelSize = 5) + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) + no_legend() 
p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "seurat_clusters", embedding = "UMAP", size = 1.4, pal =  ggplot_pal, labelMeans = FALSE) + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) + no_legend() 
#p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "seurat_clusters", embedding = "UMAP", pal = renamed_clusters_colors, size = 0.3) + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

p4 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", size = 1.4) + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
p5 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "nFrags", embedding = "UMAP", size = 1.4) + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
p6 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "NucleosomeRatio", embedding = "UMAP", size = 1.4) + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

umaps <- (p1+p2+p3) / (p4+p5+p6) + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
ggsave(filename = "umaps.png", umaps, units = "px", device = "png", width=3000,height=2000,dpi = 300)

DimPlot(DNs.RNA.seu) + p3

atac.projMetadata <- as.data.frame(proj2@cellColData)
sample_perCluster <- atac.projMetadata %>% 
  dplyr::group_by(Clusters, Sample) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n) * 100) %>% group_by(Clusters) %>% 
  dplyr::mutate(label = glue("n = {sum(n)}"))

sample_perCluster$Sample <- factor(sample_perCluster$Sample, levels = samples)
samplePerCluster_plot = ggplot() +
  geom_bar(data = sample_perCluster, aes(x = Clusters, y = freq, fill = Sample), stat="identity") +
  scale_fill_manual(values = sample_colors) + 
  theme_classic() + theme(legend.position = "right") + 
  labs(title = "Contribution of each sample to each DN-ATAC-cluster", x= "", y= "%cells") +
  geom_text(aes(Clusters, 100 + 2, label = label, fill = NULL), size = 5,  data = sample_perCluster)
ggsave(filename = "sample_perCluster_ATAC.png", samplePerCluster_plot, units = "px", device = "png", width=3000,height=2000,dpi = 300)

RNAgroup_perCluster <-atac.projMetadata %>% 
  dplyr::group_by(Clusters, seurat_clusters) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n) * 100) %>% group_by(Clusters) %>% 
  dplyr::mutate(label = glue("n = {sum(n)}"))

RNAgroup_perCluster_plot = ggplot() +
  geom_bar(data = RNAgroup_perCluster, aes(x = Clusters, y = freq, fill = seurat_clusters), stat="identity") +
  scale_fill_manual(values = ggplot_pal) + 
  theme_classic() + theme(legend.position = "right") + 
  labs(title = "Contribution of each RNA-group to each DN-ATAC-cluster", x= "", y= "%cells") +
  geom_text(aes(Clusters, 100 + 2, label = label, fill = NULL), size = 5,  data = sample_perCluster)
ggsave(filename = "RNAGroup_perCluster_ATAC.png", RNAgroup_perCluster_plot, units = "px", device = "png", width=3000,height=2000,dpi = 300)

plots_perCluster <- samplePerCluster_plot / RNAgroup_perCluster_plot
ggsave(filename = "plots_perCluster.ATAC.png", plots_perCluster, units = "px", device = "png", width=7000,height=4000,dpi = 300)


# ========== 5) RNA-ATAC WNN Integration ==========
# Before this a first peak calling was performed, using ArchR default method. 
# While it is not the final peak calling, the proj containing the peaks resulting from one of the first attempts
# was used for the proceeding with RNA-ATAC integration 
# more details regardng peak calling can be retrieved from script 3.3_first_attempts_peakCalling.R

proj3 <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/4_230601_peak_calling/archR/peakCalling/230602_ATAC_proj4.condition.rds")


# Convert ATAC UMAPs and matrices to Signac format
# Add ATAC assay to `wnn.seu`
# Run WNN using pca + lsi
# Visualize clusterings (RNA, ATAC, WNN)
##Only common cells
wnn.seu <- subset(DNs.RNA.seu, cells = proj3$cellNames)
wnn.seu@reductions$pca <- NULL
wnn.seu@reductions$umap <- NULL

#ATAC, convert from archr to signac
se <- getMatrixFromProject(proj3, "PeakMatrix", threads = 1)
se <- as(se, "SingleCellExperiment")
names(assays(se)) <- "counts"
rownames(se) <- paste0(seqnames(se), ":", start(se),"-",end(se))
seuratObj <- as.Seurat(se, data = NULL)
rD <- getReducedDims(proj3)
seuratObj[["umap.atac"]] <- CreateDimReducObject(embeddings = rD, key = "atacUMAP_", assay = DefaultAssay(seuratObj))
seuratObj

atac_counts <- seuratObj@assays[["originalexp"]]@counts

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

chrom_assay <- CreateChromatinAssay(counts  = atac_counts, sep = c(":", "-"))
#Add ATAC layer to wnn.seu
wnn.seu[["ATAC"]] <- chrom_assay

# Default reductions (to be replaced later)
# RNA analysis
DefaultAssay(wnn.seu) <- "RNA"
wnn.seu <- SCTransform(wnn.seu, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(wnn.seu) <- "ATAC"
wnn.seu <- RunTFIDF(wnn.seu)
wnn.seu <- FindTopFeatures(wnn.seu, min.cutoff = 'q0')
wnn.seu <- RunSVD(wnn.seu)
wnn.seu <- RunUMAP(wnn.seu, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


#This bit of "renaming" si what it is not working:
# #Rename: RNA
# wnn.seu@reductions$pca <- seu@reductions$pca
# wnn.seu@reductions$pca@assay.used <- "RNA"
# wnn.seu@reductions$umap.rna <- seu@reductions$umap
# wnn.seu@reductions$umap.rna@key <- "rnaUMAP_"


# #Rename: ATAC
# wnn.seu@reductions$lsi@cell.embeddings <- proj3@reducedDims$IterativeLSI@listData$matSVD
# wnn.seu@reductions$umap.atac@cell.embeddings <- as.matrix(proj3@embeddings$UMAP$df)

#save.image("221108.TEMPORAL.rds") Here the wnn.seu object is ready to go through wnn analysis.
DefaultAssay(wnn.seu) <- "RNA"
wnn.seu <- FindMultiModalNeighbors(wnn.seu, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:30))
wnn.seu <- RunUMAP(wnn.seu, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

wnn.seu <- FindClusters(wnn.seu, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.4)

p1 <- DimPlot(wnn.seu, reduction = "umap.rna", group.by = "RNA_snn_res.0.2", label = TRUE, label.size = 2.5, repel = TRUE, cols = ggplot_pal) + ggtitle("RNA-RNAclusters") + no_legend()
p2 <- DimPlot(wnn.seu, reduction = "umap.atac", group.by = "RNA_snn_res.0.2", label = TRUE, label.size = 2.5, repel = TRUE, cols = ggplot_pal) + ggtitle("ATAC-RNAclusters") + no_legend()
p3 <- DimPlot(wnn.seu, reduction = "wnn.umap", group.by = "RNA_snn_res.0.2", label = TRUE, label.size = 2.5, repel = TRUE, cols = ggplot_pal) + ggtitle("WNN-RNAclusters")

p4 <- DimPlot(wnn.seu, reduction = "umap.rna", group.by = "wsnn_res.0.4", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA-WNNClusters") + no_legend()
p5 <- DimPlot(wnn.seu, reduction = "umap.atac", group.by = "wsnn_res.0.4", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC-WNNClusters") + no_legend()
p6 <- DimPlot(wnn.seu, reduction = "wnn.umap", group.by = "wsnn_res.0.4", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN-WNNClusters")

C2cells <- as.data.frame(wnn.seu@meta.data) %>% dplyr::filter(RNA_snn_res.0.2 == "2") 
p7 <- DimPlot(wnn.seu, reduction = "umap.rna", group.by = "wsnn_res.0.4", label = TRUE, label.size = 2.5, repel = TRUE, cells.highlight	= rownames(C2cells)) + ggtitle("RNA-RNAcluster2")
p8 <- DimPlot(wnn.seu, reduction = "umap.atac", group.by = "wsnn_res.0.4", label = TRUE, label.size = 2.5, repel = TRUE, cells.highlight	= rownames(C2cells)) + ggtitle("ATAC-RNAcluster2")
p9 <- DimPlot(wnn.seu, reduction = "wnn.umap", group.by = "wsnn_res.0.4", label = TRUE, label.size = 2.5, repel = TRUE, cells.highlight	= rownames(C2cells)) + ggtitle("WNN-RNAcluster2")

#p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
(p4 + p5 + p6 & theme(plot.title = element_text(hjust = 0.5))) / (p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))) /(p7 + p8 + p9 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))

clusters_num <- 0:8
clusters_num_cols <- hue_pal()(10)
names(clusters_num_cols) <- clusters_num

weight.RNA.cluster <- VlnPlot(wnn.seu, features = "SCT.weight", group.by = 'wsnn_res.0.4', sort = TRUE, pt.size = 0.1, cols = clusters_num_cols) + NoLegend()
weight.ATAC.cluster <- VlnPlot(wnn.seu, features = "ATAC.weight", group.by = 'wsnn_res.0.4', sort = TRUE, pt.size = 0.1, cols = clusters_num_cols) + NoLegend()
weightCluster <- weight.RNA.cluster + weight.ATAC.cluster
weight.RNA.Group <- VlnPlot(wnn.seu, features = "SCT.weight", group.by = 'RNA_snn_res.0.2', sort = TRUE, pt.size = 0.1, cols = ggplot_pal) + NoLegend()
weight.ATAC.Group <- VlnPlot(wnn.seu, features = "ATAC.weight", group.by = 'RNA_snn_res.0.2', sort = TRUE, pt.size = 0.1, cols = ggplot_pal) + NoLegend()
weightGroup <- weight.RNA.Group + weight.ATAC.Group 

weight.RNA.cluster + weight.ATAC.cluster + weight.RNA.Group + weight.ATAC.Group + plot_layout(ncol = 2, nrow = 2)
weightCluster + weightGroup + plot_layout(ncol = 2, nrow = 2)

#renamed_clusters_colors2 <- c(renamed_clusters_colors, "AB" = "#D6B41C")

#Distribution of old RNAgroups per wnnCluster 
# wnnCluster_distribution <- as.data.frame(wnn.seu@meta.data) %>% 
#    group_by(wsnn_res.0.4, RNA_snn_res.0.2) %>% 
#    dplyr::summarise(n = n()) %>% 
#    mutate(freq = n / sum(n) * 100) %>% dplyr::group_by(wsnn_res.0.4, RNA_snn_res.0.2) %>% 
#    mutate(label = glue("n = {sum(n)}"))
#  wnnCluster_perGroup <- as.data.frame(wnnCluster_distribution)
#  
# wnn_plot1 = ggplot() +
#    geom_bar(data = wnnCluster_perGroup, aes(x = wsnn_res.0.4, y = n, fill = RNA_snn_res.0.2), stat="identity") +
#    scale_fill_manual(values = ggplot_pal) + 
#    theme_classic() + theme(legend.position = "right") + 
#    labs(title = "Contribution of each RNAgroup to each WNNcluster", x= "", y= "%cells") +
#    geom_text(aes(wsnn_res.0.4, 100 + 2, label = label, fill = NULL), size = 5,  data = wnnCluster_perGroup)
# wnn_plot1

#Samples
wnnCluster_distribution2 <- as.data.frame(wnn.seu@meta.data) %>% 
  select(wsnn_res.0.4, orig.ident) %>% 
  group_by(wsnn_res.0.4, orig.ident) %>% 
  dplyr::summarise(n = n()) %>% 
  mutate(freq = n / sum(n) * 100) %>% #group_by(Group) %>% 
  mutate(label = glue("n = {sum(n)}"))

wnnCluster_perSample <- as.data.frame(wnnCluster_distribution2)
wnnCluster_perSample$orig.ident <- factor(wnnCluster_perSample$orig.ident, levels = samples)
#wnnCluster_perSample$wsnn_res.0.4 <- factor(wnnCluster_perSample$wsnn_res.0.4, levels = my_levels)
wnn_plot2 = ggplot() +
  geom_bar(data = wnnCluster_perSample, aes(x = wsnn_res.0.4, y = freq, fill = orig.ident), stat="identity") +
  scale_fill_manual(values = sample_colors) + 
  theme_classic() + theme(legend.position = "bottom") + 
  labs(title = "Contribution of each sample to each WNNcluster", x= "", y= "%cells") +
  geom_text(aes(wsnn_res.0.4, 100 + 2, label = label, fill = NULL), size = 5,  data = wnnCluster_perSample)
wnn_plot2

#wnn_plot1/wnn_plot2


# #Set colors for plotting
# renamed_clusters_colors2 <- c(renamed_clusters_colors, "B3" = "#52C2C6", "AB" = "#D6B41C")
# renamed_clusters_colors2 <- renamed_clusters_colors2[my_levels]

p3 <- DimPlot(wnn.seu, reduction = "wnn.umap", group.by = "RNA_snn_res.0.2", label = TRUE, label.size = 2.5, repel = TRUE, cols = ggplot_pal) + ggtitle("WNN-RNAGroups")
p10 <- DimPlot(wnn.seu, reduction = "wnn.umap", group.by = "wsnn_res.0.4", label = TRUE, label.size = 2.5, repel = TRUE, cols = clusters_num_cols) + ggtitle("WNN-WNNClusters")
p3 + p10


# ========== 6) Metadata and Reductions Sync with ArchR ==========
# Add Seurat metadata (clusters, reductions) back into ArchRProj
# Save UMAP embeddings for RNA, ATAC, WNN in ArchR slot
seurat_dataframe <- as.data.frame(wnn.seu@meta.data) %>% rownames_to_column("cellNames") %>% 
  dplyr::select(cellNames, wsnn_res.0.4)

# Get RNA reduction:
RNAreduction <- as.data.frame(wnn.seu@reductions[["umap.rna"]]@cell.embeddings) %>% rownames_to_column("cellNames") 
colnames(RNAreduction) <- c("cellNames","custom#rnaUMAP_1","custom#rnaUMAP_2")

# Get ATAC reduction:
ATACreduction <- as.data.frame(wnn.seu@reductions[["umap.atac"]]@cell.embeddings) %>% rownames_to_column("cellNames") 
colnames(ATACreduction) <- c("cellNames", "custom#atacUMAP_1","custom#atacUMAP_2")

# Get WNN reduction:
WNNreduction <- as.data.frame(wnn.seu@reductions[["wnn.umap"]]@cell.embeddings) %>% rownames_to_column("cellNames") 
colnames(WNNreduction) <- c("cellNames","custom#wnnUMAP_1","custom#wnnUMAP_2")

seurat_dataframe_complete <- left_join(seurat_dataframe, RNAreduction, by = "cellNames") %>%
  left_join(., ATACreduction, by = "cellNames") %>% 
  left_join(., WNNreduction, by = "cellNames") %>% column_to_rownames("cellNames")

proj3 <- proj3[proj3$cellNames %in% rownames(seurat_dataframe_complete),] # sanity check
proj3_seurat <- proj3

# Add metadata: all this metadata will be added to PROJ4 (most likely?)
proj3_seurat$wsnn_res.0.4 <- wnn.seu$wsnn_res.0.4
proj3_seurat@cellColData$wsnn_res.0.4 <- as.character(proj3_seurat@cellColData$wsnn_res.0.4)
proj3_seurat$wsnn_res.0.4 <- wnn.seu$wsnn_res.0.4
proj3_seurat@cellColData$wsnn_res.0.4 <- as.character(proj3_seurat@cellColData$wsnn_res.0.4)
proj3_seurat@cellColData$wnnCluster_simple <- ifelse(proj3_seurat@cellColData$wsnn_res.0.4 )


#SIMPLIFIED
proj3_seurat@embeddings$seuratUMAP_RNA <- SimpleList(df = RNAreduction, params = list())
proj3_seurat@embeddings$seuratUMAP_ATAC <- SimpleList(df = ATACreduction, params = list())
proj3_seurat@embeddings$seuratUMAP_WNN <- SimpleList(df = WNNreduction, params = list())

#Sanity check
p1 <- plotEmbedding(ArchRProj = proj3_seurat, colorBy = "cellColData", name = "wnnCluster", embedding = "seuratUMAP_RNA", size = 0.3, pal = renamed_clusters_colors2) + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
p2 <- plotEmbedding(ArchRProj = proj3_seurat, colorBy = "cellColData", name = "wnnCluster", embedding = "seuratUMAP_ATAC", size = 0.3, pal = renamed_clusters_colors2) + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
p3 <- plotEmbedding(ArchRProj = proj3_seurat, colorBy = "cellColData", name = "wnnCluster", embedding = "seuratUMAP_WNN", size = 0.3, pal = renamed_clusters_colors2) + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
p1+p2+p3



# ========== 7) Explore SN Identity ==========
# Score SN-related gene sets (different sources tested)
# Plot feature scores and genes in WNN UMAP
SN_vedran <- list(c("Sox6","Aldh1a7","Ndnf","Serpine2","Rbp4","Fgf20"))
wnn.seu <- AddModuleScore(object = wnn.seu, features = SN_vedran, name = 'SN_Vedran_list')
p11 <- FeaturePlot(object = wnn.seu, features = 'SN_Vedran_list1', reduction="wnn.umap")
#FeaturePlot(object = wnn.seu, features = c("Slc17a6", "Aldh1a1"), reduction="wnn.umap")

SN_specific <- c("Sox6","Aldh1a7","Ndnf","Serpine2","Rbp4","Fgf20")
FeaturePlot(object = wnn.seu, features = SN_specific, reduction="wnn.umap", ncol = 3)
FeaturePlot(object = wnn.seu, features = c("Otx2", "Neurod6"), reduction="wnn.umap")
# So far it is not very clear which are SN cells. 

SN_kramer <- list(c("Aldh1a7", "Igf1", "Bsn", "Ntsr1", "Kcns3", "Anxa1"))
wnn.seu <- AddModuleScore(object = wnn.seu, features = SN_kramer, name = 'SN_kramer_list')
p12 <-FeaturePlot(object = wnn.seu, features = "SN_kramer_list1", reduction="wnn.umap")

SN_kramer <- list(c("Aldh1a7", "Igf1", "Bsn", "Ntsr1", "Kcns3", "Anxa1"))
wnn.seu <- AddModuleScore(object = wnn.seu, features = SN_kramer, name = 'SN_kramer_list')
FeaturePlot(object = wnn.seu, features = "SN_kramer_list1", reduction="wnn.umap")

SN_kramer <- c("Aldh1a7", "Igf1", "Bsn", "Ntsr1", "Kcns3", "Anxa1")
FeaturePlot(object = wnn.seu, features = SN_kramer, reduction="wnn.umap", ncol = 3)

SN_lamanno <- list(c("Sox6", "Aldh1a1", "Aldh1a7", "Anxa1", "Ndnf"))
wnn.seu <- AddModuleScore(object = wnn.seu, features = SN_lamanno, name = 'lamanno_list')
p13 <- FeaturePlot(object = wnn.seu, features = "lamanno_list1", reduction="wnn.umap")

SN_lamanno <- c("Sox6", "Aldh1a1", "Aldh1a7", "Anxa1", "Ndnf")
FeaturePlot(object = wnn.seu, features = SN_lamanno, reduction="wnn.umap", ncol=3)

p10 + p11 + p12 + p13 + plot_layout(ncol=4, nrow=1)



# ========== 8) Marker Gene & Peak Detection ==========
# Run `getMarkerFeatures()` for GeneExpressionMatrix and PeakMatrix
# Plot heatmaps with top markers per cluster

proj3_seurat@cellColData$wsnn_res.0.4 <- as.character(proj3_seurat@cellColData$wsnn_res.0.4)
proj3_seurat@cellColData$wnnCluster <- as.character(proj3_seurat@cellColData$wnnCluster)

#Gene markers:
markersGE <- getMarkerFeatures(
  ArchRProj = proj3_seurat, 
  useMatrix = "GeneExpressionMatrix", 
  bias = c("TSSEnrichment", "log10(nFrags)"),
  groupBy = "wsnn_res.0.4",
  testMethod = "wilcoxon", 
  threads=1
)

markerList <- getMarkers(markersGE, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

markerList.top5 <- as.data.frame(markerList) %>%
  group_by(group_name) %>% 
  slice_min(FDR, n= 4, with_ties = FALSE) %>% 
  ungroup() %>% 
  select(name) 

heatmapGE <- plotMarkerHeatmap(
  seMarker = markersGE, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  transpose = FALSE, 
  clusterCols = FALSE, 
  labelMarkers = markerList.top5$name
)
heatmapGE


#Peaks markers
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj3_seurat, 
  useMatrix = "PeakMatrix", 
  groupBy = "wsnn_res.0.4",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  threads=1
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList.peaks.top2 <- as.data.frame(markerList) %>%
  group_by(group_name) %>% 
  slice_min(FDR, n= 2, with_ties = FALSE) %>% 
  ungroup() %>% 
  dplyr::mutate(name = glue("{seqnames}:{start}-{end}")) %>% 
  select(name) 

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.3 & Log2FC >= 0.5",
  transpose = FALSE, 
  labelMarkers = markerList.peaks.top2$name,
  clusterCols = FALSE)


# ========== 9) ATAC Metric Distribution by WNN Cluster ==========
# Violin plots of TSS enrichment and nFrags by WNN cluster
ATACmetrics <- as.data.frame(proj3_seurat@cellColData)
ATACmetrics$wsnn_res.0.4 <- factor(ATACmetrics$wsnn_res.0.4)
ATACmetrics.TSS <- ggplot(ATACmetrics, aes(x=wsnn_res.0.4, y=TSSEnrichment, fill = wsnn_res.0.4, colour = wsnn_res.0.4))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  geom_boxplot(aes(x = as.numeric(wsnn_res.0.4)+0.25, y = TSSEnrichment),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_y_log10() +
  labs(title = "DNs: TSS enrichment", x= "", y= "") +
  # scale_fill_manual(values = renamed_clusters_colors2) + 
  # scale_color_manual(values = renamed_clusters_colors2) +
  theme_minimal() +
  theme(legend.position = "none")

ATACmetrics.nFrags <- ggplot(ATACmetrics, aes(x=wsnn_res.0.4, y=nFrags, fill = wsnn_res.0.4, colour = wsnn_res.0.4))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  geom_boxplot(aes(x = as.numeric(wsnn_res.0.4)+0.25, y = nFrags),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_y_log10() +
  labs(title = "DNs: Number of unique fragments", x= "", y= "") +
  # scale_fill_manual(values = renamed_clusters_colors2) + 
  # scale_color_manual(values = renamed_clusters_colors2) +
  theme_minimal() +
  theme(legend.position = "none") 

ATACmetrics.TSS/ATACmetrics.nFrags

# ========== 10) Export Dimensionality Reductions ========== IMPORTANT; used in further steps
# Save RNA/ATAC/WNN UMAP matrices and combined metadata as .rds
saveRDS(seurat_dataframe_complete, "seurat_dataframe_complete.rds")
saveRDS(RNAreduction, "RNAreduction.rds")
saveRDS(ATACreduction, "ATACreduction.rds")
saveRDS(WNNreduction, "WNNreduction.rds")
saveRDS(wnn.seu, "wnn.seu.rds")


p1 <- DimPlot(wnn.seu, reduction = "wnn.umap", group.by = "wsnn_res.0.4", label = FALSE, label.size = 2.5, repel = TRUE) +
  theme(legend.position="bottom", legend.direction="horizontal") +
  guides(color = guide_legend(nrow = 3, override.aes = list(size = 2))) 

p2 <- DimPlot(wnn.seu, reduction = "wnn.umap", group.by = "orig.ident", label = FALSE, label.size = 2.5, repel = TRUE, cols = sample_colors) +
  theme(legend.position="bottom", legend.direction="horizontal") +
  guides(color = guide_legend(nrow = 6, override.aes = list(size = 2))) 

p1+ p2


# ========== 11) Backup ==========
#save.image("230607_final_basic_DN_ATAC_backup.rds")



