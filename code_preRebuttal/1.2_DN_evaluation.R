# ===========================================
# Script Title: 1.2 DN Identification
# Author: Luna Zea
# Date: 2025-06-05 (for github)
# Description:
#   This script the integrated conos object obtained in 1.1, identify and 
#   characterize the the DN population, evaluate batch effects, filter poor-quality
#   cells and generate visualization plots.
# ===========================================

# ========== Load Required Libraries ==========
library(dplyr)
library(tibble)
library(ggplot2)
library(ggside)
library(ggrepel)
library(Seurat)
library(Signac)
library(harmony)
library(glue)
library(EnsDb.Mmusculus.v79)
library(ggExtra)
library(tidyverse)
library(PupillometryR)
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Cairo)
library(muscat)
library(scCustomize)
library(DESeq2)
library(RColorBrewer)
library(ComplexHeatmap)
library(GenomicRanges)
library(Pando)
library(patchwork)
library(ggalluvial)
library(ggVennDiagram)
library(conos)
library(SingleCellExperiment)

# ========== Set Environment ==========
rm(list = ls(all.names = TRUE))
gc()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

# ========== Load Data ==========
# seu.merged
seu.merged <- ("./data/250530.integrated_RNA_all.rds")


# ========== Extract DN population ==========

#Evaluate DN enrichment
vedran_markers <- read_tsv("../data/Cell_Markers_Update.txt")
dopaminess_set1 <- vedran_markers %>% dplyr::filter(cell_type == "Dopaminergic") %>% dplyr::select(gene_name) %>% pull()
dopaminess_set2 <- c(dopaminess_set1[dopaminess_set1 %notin% c("Foxa2")], "Ddc", "Pbx1", "Lmo3", "Calb1") #Foxa2 lowly detected (technical, snRNAseq)

seu.merged.dopaminess.test1 <- seu.merged
seu.merged.dopaminess.test1 <- AddModuleScore(seu.merged.dopaminess.test1, features = dopaminess_set1, name = "dopaminess_set1_")
f1 <- FeaturePlot(seu.merged.dopaminess.test1, features = "dopaminess_set1_1", label = TRUE, repel = TRUE)
v1 <- VlnPlot(seu.merged.dopaminess.test1, features = "dopaminess_set1_1")
f1 + v1+plot_layout(ncol=2)
rm(seu.merged.dopaminess.test1)

seu.merged.dopaminess.test2 <- seu.merged
seu.merged.dopaminess.test2 <- AddModuleScore(seu.merged.dopaminess.test2, features = dopaminess_set2, name = "dopaminess_set2_")
f2 <- FeaturePlot(seu.merged.dopaminess.test2, features = "dopaminess_set2_1", label = TRUE, repel = TRUE)
v2 <- VlnPlot(seu.merged.dopaminess.test2, features = "dopaminess_set2_1")
rm(seu.merged.dopaminess.test2)

#Plot DN score distributions:
f1+v1+f2+v2 +plot_layout(ncol=2)

#Subset DNs (cluster 12) and subcluster
seu.merged <- AddModuleScore(seu.merged, features = dopaminess_set2, name = "dopaminess_set2_")
DN.conos.seu <- subset(x = seu.merged, subset = conos.RNA.clusters == "12")

DN.conos.seu <- NormalizeData(DN.conos.seu)
DN.conos.seu <- ScaleData(DN.conos.seu, verbose = FALSE)
DN.conos.seu <- FindVariableFeatures(DN.conos.seu, selection.method = "vst", nfeatures = 2000)
DN.conos.seu <- RunPCA(DN.conos.seu, npcs = 30, verbose = FALSE)
DN.conos.seu <- RunUMAP(DN.conos.seu, reduction = "pca", dims = 1:30)
DN.conos.seu <- FindNeighbors(DN.conos.seu, reduction = "pca", dims = 1:30)
DN.conos.seu <- FindClusters(DN.conos.seu, resolution = 0.2)

DimPlot(DN.conos.seu)


# ========== Bacth evaluation ==========

batches <- read.csv("../data/220815_info_experimental_batches.csv")

DN.conos.seu.metadata <- DN.conos.seu@meta.data %>% mutate(
  simpleIdent = substr(orig.ident, 1, (str_length(orig.ident)-3)), 
  simpleIdent = if_else(str_detect(simpleIdent, "saline"), "saline", simpleIdent))

DN.conos.seu.metadata.complete <- DN.conos.seu.metadata %>% mutate(sample=orig.ident) %>% 
  separate(sample, c("timepoint", "treatment", "replicate"), sep = "_") %>% rownames_to_column("cellNames")

DN.conos.seu.metadata.complete <- left_join(DN.conos.seu.metadata.complete, batches, by = "orig.ident") %>% column_to_rownames("cellNames")
DN.conos.seu@meta.data <- DN.conos.seu.metadata.complete


#Visualize batches in UMAP plot:
#Set color schemes:
colors_CR1 <- c("#9F2365", "#617641", "#C48208","#326186","#AE430A", "#564686")
names(colors_CR1) <- c("m30_cocaine_R1", "h1_cocaine_R1", "h4_cocaine_R1", "h8_cocaine_R1", "h24_cocaine_R1", "d14_cocaine_R1")
colors_CR2 <- c("#6C8448", "#D78F09", "#376C95", "#C14B0B", "#6754A0")
names(colors_CR2) <- c("h1_cocaine_R2", "h4_cocaine_R2", "h8_cocaine_R2", "h24_cocaine_R2", "d14_cocaine_R2")
colors_CR3 <- c("#D4520C", "#7D6CB2")
names(colors_CR3) <- c("h24_cocaine_R3", "d14_cocaine_R3")
colors_SR1 <- c("#B2C596", "#F9CB76", "#88B2D3", "#F7A578", "#A194C7")
names(colors_SR1) <- c("h1_saline_R1", "h4_saline_R1", "h8_saline_R1", "h24_saline_R1", "d14_saline_R1")
sample_colors <- c(colors_SR1,colors_CR1,colors_CR2,colors_CR3)
#sample_colors <- sample_colors[samples]

# Simplified class
simpleIdent_colors <- c("gray", "#9F2365", "#617641", "#C48208","#326186","#AE430A", "#6728B8")
names(simpleIdent_colors) <- c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")

#Timepoints
timepoints <- c("m30", "h1", "h4", "h8", "h24", "d14")
timepoints_colors <- colors_CR1
names(timepoints_colors) <- timepoints

# ========== Bias inspection across samples ==========

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

Idents(DN.conos.seu) <- "orig.ident"
p1 <- print(DimPlot(DN.conos.seu, reduction = "umap", cols=sample_colors) + ggtitle("DNs: Samples") + NoLegend())
p1

Idents(DN.conos.seu) <- "simpleIdent"
p2 <- print(DimPlot(DN.conos.seu, reduction = "umap", cols = simpleIdent_colors) + ggtitle("DNs: Treatments and saline"))
p2

Idents(DN.conos.seu) <- "timepoint"
p3 <- print(DimPlot(DN.conos.seu, reduction = "umap", cols = timepoints_colors) + ggtitle("DNs: Time points"))
p3
p2+p3

Idents(DN.conos.seu) <- "sample_prep"
p4 <- print(DimPlot(DN.conos.seu, reduction = "umap", cols = SamplePrep_colors)  + ggtitle("DNs: Sample prep day"))
p4

Idents(DN.conos.seu) <- "FACS"
p5 <- print(DimPlot(DN.conos.seu, reduction = "umap", cols = FACS_colors)  + ggtitle("DNs: FACS machine"))
p5
p4+p5

Idents(DN.conos.seu) <- "GEXlibrary"
p6 <- print(DimPlot(DN.conos.seu, reduction = "umap", cols = GEXlibrary_colors)  + ggtitle("GEX library prep day"))
p6

Idents(DN.conos.seu) <- "Sequencing"
p7 <- print(DimPlot(DN.conos.seu, reduction = "umap", cols = sequencing_colors)  + ggtitle("DNs: Sequencing pool"))
p7
p6+p7

# Plot all sample, clusters and contribution to each cluster:
Idents(DN.conos.seu) <- "seurat_clusters"
p8 <- print(DimPlot(DN.conos.seu, reduction = "umap", label = TRUE) + ggtitle("DNs: Subclusters")) + NoLegend()
p8

(p1+p8)/(p2+p3) 
(p4+p5)/(p6+p7) 

DN.conos.seu.processed.metadata <- DN.conos.seu@meta.data #reclustered etc. 
DN.conos.seu.processed.metadata$orig.ident <- factor(DN.conos.seu.processed.metadata$orig.ident, levels = samples)

nUMIs_plot = ggplot(DN.conos.seu.processed.metadata, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident, colour = orig.ident)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust =2, trim = TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(aes(x = as.integer(orig.ident)+0.25, y = nCount_RNA), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_y_log10() +
  scale_fill_manual(values=sample_colors) +
  scale_color_manual(values=sample_colors) +
  theme_classic() + theme(legend.position = "none") + 
  geom_hline(yintercept=500, lty = "dashed", size = 0.5) +
  geom_hline(yintercept=1000, lty = "dashed", size = 0.5, col = "red") +
  geom_hline(yintercept=1500, lty = "dashed", size = 0.5) +
  labs(title = "Number of transcripts (UMIs) per sample (only DNs)", x= "", y= "# UMIs") 
nUMIs_plot

nGenes_plot = ggplot(DN.conos.seu.processed.metadata, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident, colour = orig.ident)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust =2, trim = TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(aes(x = as.integer(orig.ident)+0.25, y = nFeature_RNA), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_y_log10() +
  scale_fill_manual(values=sample_colors) +
  scale_color_manual(values=sample_colors) +
  theme_classic() + theme(legend.position = "none") + 
  geom_hline(yintercept=500, lty = "dashed", size = 0.5) +
  geom_hline(yintercept=1000, lty = "dashed", size = 0.5, col = "red") +
  geom_hline(yintercept=1500, lty = "dashed", size = 0.5) +
  labs(title = "Number of features (genes) per sample (only DNs)", x= "", y= "# genes") 
nGenes_plot
nUMIs_plot / nGenes_plot


# ========== Filter and Recluster DN Population ==========
# DN.conos.sub contains only good quality DNs (nCount_RNA >= 1000 & nFeature_RNA >= 500)

DN.conos.seu <- subset(x = seu.merged, subset = conos.RNA.clusters == "12")
DN.conos.sub <- DN.conos.seu[, DN.conos.seu$nCount_RNA >= 1000 & DN.conos.seu$nFeature_RNA >= 500]

DN.conos.sub <- NormalizeData(DN.conos.sub)
DN.conos.sub <- ScaleData(DN.conos.sub, verbose = FALSE)
DN.conos.sub <- FindVariableFeatures(DN.conos.sub, selection.method = "vst", nfeatures = 1500) 
DN.conos.sub <- RunPCA(DN.conos.sub, npcs = 30, verbose = FALSE)
DN.conos.sub <- RunUMAP(DN.conos.sub, reduction = "pca", dims = 1:30)
DN.conos.sub <- FindNeighbors(DN.conos.sub, reduction = "pca", dims = 1:30)
DN.conos.sub <- FindClusters(DN.conos.sub, resolution = 0.2)
DimPlot(DN.conos.sub)

DN.conos.sub.metadata <- DN.conos.sub@meta.data %>% mutate(
  simpleIdent = substr(orig.ident, 1, (str_length(orig.ident)-3)), 
  simpleIdent = if_else(str_detect(simpleIdent, "saline"), "saline", simpleIdent))

DN.conos.sub.metadata.complete <- DN.conos.sub.metadata %>% mutate(sample=orig.ident) %>% 
  separate(sample, c("timepoint", "treatment", "replicate"), sep = "_") %>% rownames_to_column("cellNames")

DN.conos.sub.metadata.complete <- left_join(DN.conos.sub.metadata.complete, batches, by = "orig.ident") %>% column_to_rownames("cellNames")
DN.conos.sub@meta.data <- DN.conos.sub.metadata.complete


Idents(DN.conos.sub) <- "orig.ident"
p1 <- print(DimPlot(DN.conos.sub, reduction = "umap", cols=sample_colors) + ggtitle("DNs: Samples") + NoLegend())
p1

Idents(DN.conos.sub) <- "simpleIdent"
p2 <- print(DimPlot(DN.conos.sub, reduction = "umap", cols = simpleIdent_colors) + ggtitle("DNs: Treatments and saline"))
p2

Idents(DN.conos.sub) <- "timepoint"
p3 <- print(DimPlot(DN.conos.sub, reduction = "umap", cols = timepoints_colors) + ggtitle("DNs: Time points"))
p3
p2+p3

Idents(DN.conos.sub) <- "sample_prep"
p4 <- print(DimPlot(DN.conos.sub, reduction = "umap", cols = SamplePrep_colors)  + ggtitle("DNs: Sample prep day"))
p4

Idents(DN.conos.sub) <- "FACS"
p5 <- print(DimPlot(DN.conos.sub, reduction = "umap", cols = FACS_colors)  + ggtitle("DNs: FACS machine"))
p5
p4+p5

Idents(DN.conos.sub) <- "GEXlibrary"
p6 <- print(DimPlot(DN.conos.sub, reduction = "umap", cols = GEXlibrary_colors)  + ggtitle("GEX library prep day"))
p6

Idents(DN.conos.sub) <- "Sequencing"
p7 <- print(DimPlot(DN.conos.sub, reduction = "umap", cols = sequencing_colors)  + ggtitle("DNs: Sequencing pool"))
p7
p6+p7

# Plot all sample, clusters and contribution to each cluster:
Idents(DN.conos.sub) <- "seurat_clusters"
p8 <- print(DimPlot(DN.conos.sub, reduction = "umap", label = TRUE) + ggtitle("DNs: Subclusters")) + NoLegend()
p8

(p1+p8)/(p2+p3) 
(p4+p5)/(p6+p7) 

DN.conos.sub.processed.metadata <- DN.conos.sub@meta.data #reclustered etc. 
DN.conos.sub.processed.metadata$orig.ident <- factor(DN.conos.sub.processed.metadata$orig.ident, levels = samples)

nUMIs_plot = ggplot(DN.conos.sub.processed.metadata, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident, colour = orig.ident)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust =2, trim = TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(aes(x = as.integer(orig.ident)+0.25, y = nCount_RNA), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_y_log10() +
  scale_fill_manual(values=sample_colors) +
  scale_color_manual(values=sample_colors) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.x = element_blank()) + 
  labs(title = "", x= "", y= "# UMIs") 

nGenes_plot = ggplot(DN.conos.sub.processed.metadata, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident, colour = orig.ident)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust = 2, trim = TRUE) +
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(aes(x = as.integer(orig.ident) + 0.25, y = nFeature_RNA), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_y_log10(
    breaks = c(1000, 3000, 10000),               # Specify the breaks for ticks
    labels = c("1000", "3000", "10000"),         # Provide labels for the breaks
    limits = c(1000, 10000)                      # Set limits to include all breaks
  ) +
  scale_fill_manual(values = sample_colors) +
  scale_color_manual(values = sample_colors) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank()               # Remove x-axis tick labels
  ) +
  guides(fill = guide_legend(nrow = 2), colour = guide_legend(nrow = 2)) +
  labs(title = "", x = "", y = "# genes")


UMIs_and_features <- nUMIs_plot / nGenes_plot

#Save figure
# dev.size()
# width_original = 14.697917
# height_original= 7.708333
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter1_DNs_nUMIs_nFeatures"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# UMIs_and_features
# dev.off()


# ========== Evaluate Dopaminergic Signature ==========

DN.conos.sub.processed.metadata <- DN.conos.sub.processed.metadata %>% 
  mutate(DN_module_score = rowSums(select(., dopaminess_set2_1:dopaminess_set2_11))/ 11)

DN_module_plot_sample = ggplot(DN.conos.sub.processed.metadata, aes(x = orig.ident, y = DN_module_score, fill = orig.ident, colour = orig.ident)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust =2, trim = TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(aes(x = as.integer(orig.ident)+0.25, y = DN_module_score), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  #scale_y_log10() +
  scale_fill_manual(values=sample_colors) +
  scale_color_manual(values=sample_colors) +
  theme_classic() + theme(legend.position = "none") + 
  labs(title = "DN module score per sample", x= "", y= "DN module score") 
DN_module_plot_sample

DN_module_plot_RNAgroup = ggplot(DN.conos.sub.processed.metadata, aes(x = seurat_clusters, y = DN_module_score, fill = seurat_clusters, colour = seurat_clusters)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust =2, trim = TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(aes(x = as.integer(seurat_clusters)+0.25, y = DN_module_score), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  #scale_y_log10() +
  #scale_fill_manual(values=sample_colors) +
  #scale_color_manual(values=sample_colors) +
  theme_classic() + theme(legend.position = "none") + 
  labs(title = "DN module score per DN subcluster", x= "", y= "DN module score") 
DN_module_plot_RNAgroup

DN_module_plot_sample/DN_module_plot_RNAgroup

#Samples over/underrepresented?
sample_perCluster <- DN.conos.sub.processed.metadata %>% 
  dplyr::group_by(seurat_clusters, orig.ident) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n) * 100) %>% group_by(seurat_clusters) %>% 
  dplyr::mutate(label = glue("n = {sum(n)}"))


sample_perCluster <- as.data.frame(sample_perCluster)
samplePerCluster_plot = ggplot() +
  geom_bar(data = sample_perCluster, aes(x = seurat_clusters, y = freq, fill = orig.ident), stat="identity") +
  scale_fill_manual(values = sample_colors) + 
  theme_classic() + theme(legend.position = "right") + 
  labs(title = "Contribution of each sample to each DN-cluster", x= "", y= "%cells") +
  geom_text(aes(seurat_clusters, 100 + 2, label = label, fill = NULL), size = 5,  data = sample_perCluster)

# ========== Back-up files ==========

#write_tsv(DN.conos.sub.processed.metadata %>% rownames_to_column("cellNames"), ../results/DN.conos.sub.processed.metadata.tsv", col_names = NA)
#saveRDS(DN.conos.sub, "../results/DN.conos.sub.rds")




