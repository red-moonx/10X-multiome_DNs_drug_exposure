# ===========================================
# Script Title: 1.2 Evaluation of the DN specific cluster
# Author: Luna Zea Redondo
# Date: 2026-05-02
# Description:
#   This script explores the DN specific cluster and yield a final set of high quality DNs.
#   Includes: removal of low quality cells, SN putative DNs, and other sources of contamination 
# ===========================================


# ========== Set Environment ==========
rm(list = ls(all.names = TRUE))
gc()
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

# ========== Load Required Libraries ==========
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(PupillometryR)
library(patchwork)
library(glue)
library(readr)
library(tidyr)
library(pheatmap)

source("/fast/AG_Pombo/luna/2026_rebuttal/config.R", local = FALSE)

# ========== Set wd and load data ==========
setwd("/fast/AG_Pombo/luna/2026_rebuttal/2_DN-evaluation")
seurat_DNs <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/1_conos/seu.cluster13_putativeDNs.rds")


# ----------------------------
# 1. Remove low quality cells
# ----------------------------
message(paste("DN count before filtering:", ncol(seurat_DNs))) #2,761

#Filtering
seurat_DNs <- subset(seurat_DNs, subset = nCount_RNA >= 1000 & nFeature_RNA >= 500) 
message(paste("DN count after filtering:", ncol(seurat_DNs))) #2,597

#Plot: midbrain DNs
cols_to_keep <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "conos_clusters")
seurat_DNs@meta.data <- seurat_DNs@meta.data[, cols_to_keep, drop = FALSE]

seurat_DNs.processed.metadata <- seurat_DNs@meta.data
seurat_DNs.processed.metadata$orig.ident <- factor(
  seurat_DNs.processed.metadata$orig.ident,
  levels = sample_order
)

nUMIs_plot <- ggplot(seurat_DNs.processed.metadata, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident, colour = orig.ident)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust = 2, trim = TRUE) +
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.3, colour = "black", position = position_nudge(x = .25)) +
  scale_y_log10() +
  scale_fill_manual(values = sample_colors) +
  scale_color_manual(values = sample_colors) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_hline(yintercept=1000, lty = "dashed", size = 0.5, col = "red") +
  labs(title = "Number of transcripts (UMIs) per sample (only DNs)", x = "", y = "# UMIs")

nGenes_plot = ggplot(seurat_DNs.processed.metadata, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident, colour = orig.ident)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust =2, trim = TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(aes(x = as.integer(orig.ident)+0.25, y = nFeature_RNA), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  scale_y_log10() +
  scale_fill_manual(values=sample_colors) +
  scale_color_manual(values=sample_colors) +
  theme_classic() + theme(legend.position = "none") + 
  geom_hline(yintercept=500, lty = "dashed", size = 0.5, col = "red") +
  labs(title = "Number of features (genes) per sample (only DNs)", x= "", y= "# genes") 

nUMIs_plot / nGenes_plot


# -------------------------------------
# 2. Is there SN derived contamination?
# -------------------------------------

# 2.A. Embedding for midbrain DNs
# -------------------------------
seurat_DNs <- NormalizeData(seurat_DNs)
seurat_DNs <- FindVariableFeatures(seurat_DNs, selection.method = "vst", nfeatures = 2000)
seurat_DNs <- ScaleData(seurat_DNs)
seurat_DNs <- RunPCA(seurat_DNs, npcs = 20)

seurat_DNs <- FindNeighbors(seurat_DNs, dims = 1:20)
seurat_DNs <- FindClusters(seurat_DNs, resolution = 0.4)

seurat_DNs <- RunUMAP(seurat_DNs, dims = 1:20, min.dist = 0.4, n.neighbors = 50)
DimPlot(seurat_DNs)


# 2.B. Overlap SN markers and apply cutoff for SN-derived DNs 
# ------------------------------------------------------------
#based on the SN-DN list of genes in Szabo et al. 2024 (biorxiv)
SN_markers_szabo <- c("Sox6","Aldh1a7","Ndnf","Serpine2","Rbp4","Fgf20")

## 1. Calculate the score 
seurat_DNs <- AddModuleScore(
  object = seurat_DNs, 
  features = list(SN_markers_szabo = SN_markers_szabo), 
  name = "SN_markers_szabo"
)

# Plot:
SN_featurePlot <- FeaturePlot(seurat_DNs, features = "SN_markers_szabo1")

## 2. Classify based on your > 0 threshold
seurat_DNs$region <- ifelse(seurat_DNs$SN_markers_szabo1 > 0, "SN", "VTA")

## 3. Extract only the VTA
seurat_vta_clean <- subset(seurat_DNs, subset = region == "VTA")

## 4. Check numbers
message(paste("Filtered out", sum(seurat_DNs$region == "SN"), "SN cells."))
DimPlot(seurat_vta_clean)



# Plot (save for manuscript)
SN_featurePlot <- FeaturePlot(seurat_DNs, features = "SN_markers_szabo1")

SN_vedran1_plot <- ggplot(seurat_DNs@meta.data, aes(x = RNA_snn_res.0.4, y = SN_markers_szabo1, fill = RNA_snn_res.0.4, colour = RNA_snn_res.0.4)) +
  geom_flat_violin(position = position_nudge(x = 0.25, y = 0), adjust = 2, trim = TRUE) +
  geom_point(position = position_jitter(width = 0.15), size = 0.25) +
  geom_boxplot(aes(x = as.integer(as.factor(RNA_snn_res.0.4)) + 0.25, y = SN_markers_szabo1),
               outlier.shape = NA, alpha = 0.3, width = 0.1, colour = "BLACK") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +  # Add dashed red line
  # Uncomment and customize the following for manual color scales, if needed
  # scale_fill_manual(values = class_colors) +
  # scale_color_manual(values = class_colors) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "SN_markers_szabo1 values by RNA_snn_res.0.4", x = "", y = "SN module score") 

final_combined_plot <- SN_featurePlot + SN_vedran1_plot +
  plot_layout(nrow = 1, ncol=2, widths = c(1, 2))

#Save plot
# dev.size()
# width_original = 12.479167
# height_original= 4.458333
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2026_rebuttal/figures"
# plot_name <- "chapter1_SNcontamination"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# final_combined_plot
# dev.off()




# ----------------------------------------------------
# 3. Which is the cellular identity of each subcluster
# ----------------------------------------------------
# After testing several resolutions, we are keeping 0.04
backup <- seurat_vta_clean

# 3.a. Re-run embedding for VTA-DNs only
# ----------------------------------------
seurat_vta_clean <- NormalizeData(seurat_vta_clean) #1725 cells
seurat_vta_clean <- FindVariableFeatures(seurat_vta_clean, selection.method = "vst", nfeatures = 2000)
seurat_vta_clean <- ScaleData(seurat_vta_clean)
seurat_vta_clean <- RunPCA(seurat_vta_clean, npcs = 20)
seurat_vta_clean <- FindNeighbors(seurat_vta_clean, dims = 1:20)

#Run UMAP
seurat_vta_clean <- RunUMAP(
  seurat_vta_clean,
  dims = 1:15,
  min.dist = 0.6,
  n.neighbors = 60
)

# Cluster with resolution = 0.04
seurat_vta_clean <- FindClusters(
  seurat_vta_clean,
  resolution = 0.4,
  cluster.name = "res_0_4"
)

# Set identities
Idents(seurat_vta_clean) <- "res_0_4"
DimPlot(seurat_vta_clean)

# 3.b. Explore the subclusters identity and filter
# ------------------------------------------------
#1. Calculate marker genes
markers_res_0_4 <- FindAllMarkers(seurat_vta_clean, assay = "RNA", slot = "data", only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.25)

top10_res_0_4 <- markers_res_0_4 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 15, with_ties = FALSE) %>%
  ungroup()

#save results
#write.csv(top10_res_0_4, "clusters0to8_top15_markers_res_0_4.csv", row.names = FALSE)

# 2. Filtering
# Clusters 6 and 8 represent contamination (oligo/myelinating & serotonergic, respectively).
# Remove, recluster and calculate markers again. 
seu.VTA_DNs <- subset(
  seurat_vta_clean,
  idents = c("0", "1", "2", "3", "4", "5", "7")
)

#3. Re-run the full workflow on the cleaned VTA DN object
seu.VTA_DNs <- NormalizeData(seu.VTA_DNs)
seu.VTA_DNs <- FindVariableFeatures(seu.VTA_DNs, selection.method = "vst", nfeatures = 2000)
seu.VTA_DNs <- ScaleData(seu.VTA_DNs)
seu.VTA_DNs <- RunPCA(seu.VTA_DNs, npcs = 20)
seu.VTA_DNs <- FindNeighbors(seu.VTA_DNs, dims = 1:15)

seu.VTA_DNs <- RunUMAP(
  seu.VTA_DNs,
  dims = 1:8,
  min.dist = 0.8,
  spread = 1,       # just slightly up from 0.3
  n.neighbors = 50
)


seu.VTA_DNs <- FindClusters(
  seu.VTA_DNs,
  resolution = 0.2,
  cluster.name = "res_0_2"
)

Idents(seu.VTA_DNs) <- "res_0_2"
DimPlot(seu.VTA_DNs, label = TRUE, repel = TRUE)

###################
#Add some metadata
###################
# Define top 1% highest-UMI cells
umi_cutoff <- quantile(seu.VTA_DNs$nCount_RNA, 0.99)

high_umi_cells <- rownames(seu.VTA_DNs@meta.data)[
  seu.VTA_DNs$nCount_RNA > umi_cutoff
]

# Add metadata column
seu.VTA_DNs$high_UMI <- ifelse(
  colnames(seu.VTA_DNs) %in% high_umi_cells,
  "High_UMI",
  "Normal"
)

seu.VTA_DNs[["percent.mt"]] <- PercentageFeatureSet(
  seu.VTA_DNs,
  pattern = "^mt-"
)

seu.VTA_DNs$simpleIdent <- ifelse(
  grepl("saline", seu.VTA_DNs$orig.ident),
  "saline",
  sub("_R[0-9]+$", "", seu.VTA_DNs$orig.ident)
)


#4. Calculate marker genes for all VTA-DN clusters
markers <- FindAllMarkers(seu.VTA_DNs, assay = "RNA", slot = "data", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top 10 markers per cluster
top10_markers <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup() 

# Save
# write.csv(top10_markers, "seu.VTA_DNs_top10_markers_res0.2.csv", row.names = FALSE)


#5. Can we find DN-only and DN-combinatorial? According to Fitzgerald 2026
DA_only <- c("Sv2c", "Slc6a3", "Slc18a2", "Th", "Gch1", "Ddc", "Bdnf", "Vip")
combinatorial <- c("Reln", "Slc32a1", "Slc17a6", "Gad2", "Hcn1", "Hcn2", "Slc26a7")

FeaturePlot(seu.VTA_DNs, features = DA_only, ncol = 4, order = TRUE)
FeaturePlot(seu.VTA_DNs, features = combinatorial, ncol = 4, order = TRUE)


# Now that we have the final seurat object gfor VTA-DNs, lets prepare the clean metadata to evaluate batches
seu.VTA_DNs@meta.data <- seu.VTA_DNs@meta.data %>%
  dplyr::select(-RNA_snn_res.0.4,-res_0_4,-res_0_2)



#6. Last: check nCounts and nFeatures 
# Table
aggregate(cbind(nCount_RNA, nFeature_RNA) ~ seurat_clusters, 
          data = seu.VTA_DNs@meta.data, FUN = median)

# Violin plots
p1 <- ggplot(seu.VTA_DNs@meta.data, 
             aes(x = factor(seurat_clusters), 
                 y = nCount_RNA, 
                 fill = factor(seurat_clusters))) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "nCount_RNA") +
  theme(legend.position = "none")

p2 <- ggplot(seu.VTA_DNs@meta.data, 
             aes(x = factor(seurat_clusters), 
                 y = nFeature_RNA, 
                 fill = factor(seurat_clusters))) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "nFeature_RNA") +
  theme(legend.position = "none")

p3 <- ggplot(seu.VTA_DNs@meta.data, 
             aes(x = factor(seurat_clusters), 
                 y = percent.mt, 
                 fill = factor(seurat_clusters))) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  theme_classic() +
  labs(x = "Cluster", y = "% Mitochondrial") +
  theme(legend.position = "none")

p1 / p2 / p3


# --------------------
# 4. Batch evaluation
# --------------------
# 4.1. Add injection and experimental batches:
injection_batches <- read_csv("/fast/AG_Pombo/luna/2026_rebuttal/drug_treatment_batches.csv")

experimental_batches <- read_csv("/fast/AG_Pombo/luna/2026_rebuttal/experimental_batches.csv")
experimental_batches$sample_prep <- format(
  experimental_batches$sample_prep,
  "%d/%m/%Y"
)

experimental_batches$GEXlibrary <- format(
  experimental_batches$GEXlibrary,
  "%d/%m/%Y"
)
# Keep cell names during joins

md <- seu.VTA_DNs@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::left_join(injection_batches, by = "orig.ident") %>%
  dplyr::left_join(experimental_batches, by = "orig.ident") %>%
  tibble::column_to_rownames("cell")

seu.VTA_DNs@meta.data <- md

seu.VTA_DNs@meta.data <- seu.VTA_DNs@meta.data %>%
  mutate(Time_injection = format(
    as.POSIXct(Time_injection, format = "%H:%M:%S"),
    "%H:%M"
  ))


# 4.2. Make plot and save manuscript-ready figure:
Idents(seu.VTA_DNs) <- "orig.ident"
p1 <- DimPlot(seu.VTA_DNs, reduction = "umap", group.by = "orig.ident",  cols = sample_colors) +
  ggtitle("Samples") + theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

Idents(seu.VTA_DNs) <- "sample_prep"
p2 <- print(DimPlot(seu.VTA_DNs, reduction = "umap", cols = SamplePrep_colors)) +
  ggtitle("Sample preparation day") + theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

Idents(seu.VTA_DNs) <- "FACS"
p3 <- print(DimPlot(seu.VTA_DNs, reduction = "umap", cols = SamplePrep_colors))  + 
  ggtitle("FACS instrument") + theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

Idents(seu.VTA_DNs) <- "GEXlibrary"
p4 <- print(DimPlot(seu.VTA_DNs, reduction = "umap", cols = SamplePrep_colors)) +
  ggtitle("GEX library preparation day") + theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

Idents(seu.VTA_DNs) <- "Sequencing"
p5 <- print(DimPlot(seu.VTA_DNs, reduction = "umap", cols = c(SamplePrep_colors[2], SamplePrep_colors[1]))) +
  ggtitle("Sequencing pool") + theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")
            
experimental_batches <- print(p1 + p2 + p3 + p4 + p5) + plot_layout(nrow=1)
#ggsave(filename = "experimental_batches.pdf", plot = experimental_batches, device = "pdf", width = 22.85, height = 5.55, dpi = 300)

Idents(seu.VTA_DNs) <- "Cage_Number"
p1 <- print(DimPlot(seu.VTA_DNs, reduction = "umap", cols=SamplePrep_colors)) +
  ggtitle("Cage number") + theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

Idents(seu.VTA_DNs) <- "Age"
p2 <- print(DimPlot(seu.VTA_DNs, reduction = "umap", cols = SamplePrep_colors)) + 
  ggtitle("Age (weeks)") + theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

Idents(seu.VTA_DNs) <- "Time_injection"
p3 <- print(DimPlot(seu.VTA_DNs, reduction = "umap")) +
  ggtitle("Time of injection") + theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

Idents(seu.VTA_DNs) <- "Date_Dissection"
p4 <- print(DimPlot(seu.VTA_DNs, reduction = "umap", cols = SamplePrep_colors)) +
  ggtitle("Date of dissection") + theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

Idents(seu.VTA_DNs) <- "Dissection_window"
p5 <- print(DimPlot(seu.VTA_DNs, reduction = "umap", cols = SamplePrep_colors[c(1:5, 7)])) +
  ggtitle("Dissection window") + theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal")

injection_batches_plot <- print(p1 + p2 + p3 + p4 + p5) + plot_layout(nrow=1)
#ggsave(filename = "injection_batches.pdf", plot = injection_batches_plot, device = "pdf", width = 22.85, height = 5.55, dpi = 300)


# -------------------------------
# 5. Cell composition per cluster
# -------------------------------
seu.VTA_DNs_metadata <- seu.VTA_DNs@meta.data 
sample_perCluster <- seu.VTA_DNs_metadata %>% 
  dplyr::group_by(seurat_clusters, orig.ident) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n) * 100) %>% group_by(seurat_clusters) %>% 
  dplyr::mutate(label = glue("{sum(n)}"))
sample_perCluster <- as.data.frame(sample_perCluster)

sample_perCluster$orig.ident <- factor(sample_perCluster$orig.ident, levels = sample_order)
samplePerCluster_plot = ggplot() +
  geom_bar(data = sample_perCluster, aes(x = seurat_clusters, y = freq, fill = orig.ident), stat="identity") +
  scale_fill_manual(values = sample_colors) + 
  theme_classic() + theme(legend.position = "right") + 
  labs(title = "Contribution of each sample to each DN-cluster", x= "", y= "%cells") +
  geom_text(aes(seurat_clusters, 100 + 2, label = label, fill = NULL), size = 4,  data = sample_perCluster) 
 
# dev.size()
# width_original = 8.572917
# height_original= 5.604167
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2026_rebuttal/figures/"
# plot_name <- "chapter1_VTA_DNs_barPlot"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# samplePerCluster_plot
# dev.off()


seu.VTA_DNs$simpleIdent <- ifelse(
  grepl("saline", seu.VTA_DNs_metadata$orig.ident),
  "saline",
  sub("_R[0-9]+$", "", seu.VTA_DNs_metadata$orig.ident)
)

View(table(seu.VTA_DNs_metadata$seurat_clusters, seu.VTA_DNs_metadata$simpleIdent))

Idents(seu.VTA_DNs) <- "seurat_clusters"
DotPlot(
  seu.VTA_DNs,
  features = c("Slc18a2","Slc6a3","Ddc", "Th","Drd2", "Gad1", "Gad2", "Slc32a1", "Slc6a1", "Slc17a6")) + RotatedAxis()


# ----------------------------------------
# Counts and proportions per cluster/time
# ---------------------------------------
count_pivot <- seu.VTA_DNs@meta.data %>%
  count(seurat_clusters, simpleIdent) %>%
  pivot_wider(
    names_from = simpleIdent,
    values_from = n,
    values_fill = 0
  )

count_pivot

prop_pivot <- seu.VTA_DNs@meta.data %>%
  count(seurat_clusters, simpleIdent) %>%
  group_by(seurat_clusters) %>%
  mutate(prop = round(n / sum(n) * 100, 2)) %>%
  select(-n) %>%
  pivot_wider(
    names_from = simpleIdent,
    values_from = prop,
    values_fill = 0
  )

prop_pivot


# Chi-square test
# ----------------------------------------
timepoint_table <- table(
  seu.VTA_DNs$seurat_clusters,
  seu.VTA_DNs$simpleIdent
)

chisq_test <- chisq.test(timepoint_table)
chisq_test


# 1) Per-cell QC table
qc_per_cell <- seu.VTA_DNs@meta.data %>%
  rownames_to_column("cell") %>%
  dplyr::select(cell, orig.ident, nCount_RNA, nFeature_RNA)

# 2) Per-sample QC summary
qc_per_sample <- qc_per_cell %>%
  group_by(orig.ident) %>%
  summarise(
    n_cells = n(),
    median_UMI = median(nCount_RNA),
    median_genes = median(nFeature_RNA),
    min_UMI = min(nCount_RNA),
    max_UMI = max(nCount_RNA),
    min_genes = min(nFeature_RNA),
    max_genes = max(nFeature_RNA),
    IQR_UMI = IQR(nCount_RNA),
    IQR_genes = IQR(nFeature_RNA),
    mean_UMI = mean(nCount_RNA),
    mean_genes = mean(nFeature_RNA),
    .groups = "drop"
  )

qc_per_cell
qc_per_sample


# -------------------
# 6. Data Archiving
# -------------------
# seuVTA object (for muscat etc. )
# saveRDS(seu.VTA_DNs, "250507_seu.VTA_DNs.rds")

# write.csv(qc_per_cell, "qc_per_cell.csv", row.names = FALSE)
# write.csv(qc_per_sample, "qc_per_sample.csv", row.names = FALSE)
# load("260513_DN-evaluation.rds")

# --------------------------
# 6. Figures for manuscript:
# --------------------------

#1. Unbiased clustering p0:
p0 <- DimPlot(seu.VTA_DNs, label = FALSE) + NoAxes()

samplePerCluster_plot <- ggplot() +
  geom_bar(data = sample_perCluster,
           aes(x = seurat_clusters, y = freq, fill = orig.ident),
           stat = "identity") +
  geom_text(data = dplyr::distinct(sample_perCluster, seurat_clusters, label),
            aes(x = seurat_clusters, y = 102, label = label),
            size = 4) +
  scale_fill_manual(values = sample_colors) +
  theme_classic() +
  theme(legend.position = "right") +
  labs(title = "Contribution of each sample to each DN-cluster",
       x = "", y = "%cells")

# p0 + samplePerCluster_plot
# dev.size()
# width_original = 12.364583
# height_original= 5.447917
# 
# plot_name <- "chapter1_VTA_DNs_barPlot.pdf"
# 
# pdf(plot_name, width = width_original, height =height_original)
# p0 + samplePerCluster_plot
# dev.off()

# 2. Correlation with literature subtypes
genes_poulin <- c("Th", "Aldh1a1", "Sox6", "Slc17a6", "Otx2", "Slc32a1", "Vip", "Gad2", "Cck", "Calb1")
poulin_markers <- FeaturePlot(seu.VTA_DNs, features = genes_poulin, ncol = 5, order = TRUE, cols = c("lightgrey", "red")) & NoAxes() & NoLegend()

# pdf("chapter1_VTA_DNs_FeaturePlots.pdf", width = 9.604167 * 2, height = 5.447917 *2)
# poulin_markers
# dev.off()

#and conclusion:
dopalit <- c("Th", "Aldh1a1", "Sox6", "Slc17a6", "Otx2", "Slc32a1", "Vip", "Gad2", "Cck", "Calb1")
meta_dopa <- seu.VTA_DNs@meta.data
meta_dopa_mark_mat = meta_dopa %>% 
  cbind(
    GetAssayData(seu.VTA_DNs, "RNA")[dopalit,] %>%
      as.matrix() %>%
      t() %>%
      as.data.frame()) %>%
  # setnames(gids$gene_name)) %>%
  mutate(dopa_literature_type = case_when(
    Vip > 0 ~ "DN-6", #Vip+
    Aldh1a1 > 0 & Sox6 > 0 ~ "DN-1", #SN-1
    #Slc17a6 > 0 & Aldh1a1 == 0 ~ "#feca33", #SN-2
    Aldh1a1 == 0 & Sox6 > 0 ~ "DN-2",  #type2
    Aldh1a1 > 0 & Slc17a6 > 0 ~ "DN-5", #type3
    Aldh1a1 == 0 &  Slc17a6 > 0 ~ "DN-3",
    Slc32a1 > 0 | Gad2 > 0 ~ "DN-4",
    Cck | Calb1 ~ "DN-7",
    TRUE ~ as.character(NA)
  ))

literature.colors <- c("DN-1"="#b31700", "DN-2"="#f4822b", "DN-3"="#0075ba", "DN-4"="#55c1ff",
                       "DN-5"="#28B01D", "DN-6"="#97185d", "DN-7"="#0F8029", "NA"="gray90")

seu.VTA_DNs@meta.data <- meta_dopa_mark_mat
#How many cells:
table(seu.VTA_DNs$dopa_literature_type)
table(is.na(seu.VTA_DNs$dopa_literature_type))

#Make the plots:
Idents(seu.VTA_DNs) <- "dopa_literature_type"
literature.subtypes.plot <- DimPlot(seu.VTA_DNs, group.by = "dopa_literature_type", cols = literature.colors) + NoAxes() + theme_void()

#3. Distribution of "simpleIdents"
Idents(seu.VTA_DNs) <- "simpleIdent"
distribution_simpleIdents <-  DimPlot(seu.VTA_DNs, cols = simpleIdent_colors) + NoAxes()
# 
# pdf("chapter1_VTA_DNs_annotations.pdf", width = 12.36, height = 5.45)
# literature.subtypes.plot + distribution_simpleIdents
# dev.off()


