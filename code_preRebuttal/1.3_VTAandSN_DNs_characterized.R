# ===========================================
# Script Title: 1.3 VTA vs SN 
# Author: Luna Zea
# Date: 2025-06-05 (for github)
# Description:
#   This script first identifies and separates substantia nigra (SN) cells from 
#   ventral tegmental area (VTA) cells. It then characterizes in detial the VTA-DN population.
#   Since dissections targeted the VTA, the SN population likely represents partial contaminatio
#   and is not comprehensively profiled. For this reason, SN-specific results are not interpreted further.
# ===========================================

# ========== Set Environment ==========
rm(list = ls(all.names = TRUE))
gc()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

# ========== Load Required Libraries ==========
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(glue)
library(patchwork)
library(readr)

# ========== Set Working Directory ==========
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly"
setwd(dir)

# ========== Load Seurat Object ==========
#DNs.RNA.seu <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/230421.DN.conos.sub.rds") original path
DNs.RNA.seu <- readRDS("../data/DN.conos.sub.rds")
DNs.RNA.seu@meta.data[5:15] <- NULL

# ========== Define Sample Colors ==========
colors_CR1 <- c("#9F2365", "#617641", "#C48208","#326186","#AE430A", "#564686")
names(colors_CR1) <- c("m30_cocaine_R1", "h1_cocaine_R1", "h4_cocaine_R1", "h8_cocaine_R1", "h24_cocaine_R1", "d14_cocaine_R1")
colors_CR2 <- c("#6C8448", "#D78F09", "#376C95", "#C14B0B", "#6754A0")
names(colors_CR2) <- c("h1_cocaine_R2", "h4_cocaine_R2", "h8_cocaine_R2", "h24_cocaine_R2", "d14_cocaine_R2")
colors_CR3 <- c("#D4520C", "#7D6CB2")
names(colors_CR3) <- c("h24_cocaine_R3", "d14_cocaine_R3")
colors_SR1 <- c("#B2C596", "#F9CB76", "#88B2D3", "#F7A578", "#A194C7")
names(colors_SR1) <- c("h1_saline_R1", "h4_saline_R1", "h8_saline_R1", "h24_saline_R1", "d14_saline_R1")
sample_colors <- c(colors_SR1, colors_CR1, colors_CR2, colors_CR3)

# ========== 1. Increase Clustering Resolution ==========
resolutions <- c(3.2, 3.6, 4, 4.4)

DNs.RNA.seu <- FindNeighbors(DNs.RNA.seu, dims = 1:10)
DNs.RNA.seu <- FindClusters(DNs.RNA.seu, resolution = resolutions)

# ========== 2. Overlay Gene Lists and Sample Contribution ==========
#Contribution of each sample to each VTA-DN cluster
samples <- c("h1_saline_R1","h4_saline_R1", "h8_saline_R1", "h24_saline_R1", "d14_saline_R1", #salines
             "m30_cocaine_R1", "h1_cocaine_R1", "h1_cocaine_R2", "h4_cocaine_R1", "h4_cocaine_R2", "h8_cocaine_R1", "h8_cocaine_R2", #ETP
             "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3") #LTP

DNs.RNA.seu_metadata <- DNs.RNA.seu@meta.data 
sample_perCluster <- DNs.RNA.seu_metadata %>% 
  dplyr::group_by(RNA_snn_res.4.4, orig.ident) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n) * 100) %>% group_by(RNA_snn_res.4.4) %>% 
  dplyr::mutate(label = glue("n = {sum(n)}"))
sample_perCluster <- as.data.frame(sample_perCluster)

sample_perCluster$orig.ident <- factor(sample_perCluster$orig.ident, levels = samples)
samplePerCluster_plot = ggplot() +
  geom_bar(data = sample_perCluster, aes(x = RNA_snn_res.4.4, y = freq, fill = orig.ident), stat="identity") +
  scale_fill_manual(values = sample_colors) + 
  theme_classic() + theme(legend.position = "right") + 
  labs(title = "Contribution of each sample to each DN-cluster", x= "", y= "%cells") +
  geom_text(aes(RNA_snn_res.4.4, 100 + 2, label = label, fill = NULL), size = 2,  data = sample_perCluster) +
  theme(legend.position="bottom")
samplePerCluster_plot + no_legend()

#Plots
DimPlot(DNs.RNA.seu) + no_legend() + samplePerCluster_plot + no_legend() + plot_layout(width = c(1,4))

#Gene lists (not all have been explored in detail)
#SN genes:
SN_vedran <- c("Sox6","Aldh1a7","Ndnf","Serpine2","Rbp4","Fgf20")
SN_kramer <- c("Aldh1a7", "Igf1", "Bsn", "Ntsr1", "Kcns3", "Anxa1")
SN_lamanno <- c("Sox6", "Aldh1a1", "Aldh1a7", "Anxa1", "Ndnf")

#IEG lists
IEGs_104genes <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/IEGs_vedran.txt", col_names = FALSE) %>% pull(X1)
IEGs_29genes <- c("Egr2", "Cyr61", "Egr4", "Fos", "Arc", "Atf3", "Junb", "Fosb", "Dusp1", "Nr4a1",
                  "Dusp5", "Npas4", "Dusp6", "Elovl1", "Egr1", "Hdc", "Ier2", "Klf10", "Cartpt", "Cirbp",
                  "Drd2", "Enpp6", "Gad67", "Homer2", "Ier5", "Igsf1", "Irs2", "Jun", "Ngfr")
IEGs_cluster0_vedran <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230818.IEG_Cluster_Markers.tsv") %>% pull(gene_name)

ARG_list <- read_tsv("/fast/AG_Pombo/luna/8_defaultARCreference/4_refineDN_set/2_DEA_muscat_improved/ARG_list_tyssowski.tsv")
colnames(ARG_list) <- c("gene", "ARG_class")
ARG_list <- ARG_list %>%filter(ARG_class %in% c("rPRG", "dPRG", "SRG"))
ARG_colors <- c("rPRG" = "#FDBD0D", "dPRG" = "#B94227", "SRG" = "#087A87","no" = "gray")

rPRG_tissowsky <- ARG_list %>% dplyr::filter(ARG_class == "rPRG") %>% pull(gene)
dPRG_tissowsky <- ARG_list %>% dplyr::filter(ARG_class == "dPRG") %>% pull(gene)
SRG_tissowsky <- ARG_list %>% dplyr::filter(ARG_class == "SRG") %>% pull(gene)


#Add gene lists to metadata
gene_lists <- list(SN_vedran = SN_vedran, SN_kramer = SN_kramer, SN_lamanno = SN_lamanno, 
                   IEGs_104genes = IEGs_104genes, IEGs_29genes = IEGs_29genes, IEGs_cluster0_vedran = IEGs_cluster0_vedran, 
                   rPRG_tissowsky = rPRG_tissowsky, dPRG_tissowsky = dPRG_tissowsky, SRG_tissowsky = SRG_tissowsky)

signature.names = names(gene_lists)

DNs.RNA.seu <- AddModuleScore(DNs.RNA.seu, features = gene_lists, name = signature.names, search = TRUE)

# ========== 3. Classify VTA vs SN based on marker scores ==========
# is there SN contamination?
SN_markers_umap<- FeaturePlot(DNs.RNA.seu, features = paste0(signature.names[1:2], c(1:2)), label = TRUE, repel = TRUE, ncol = 1)
SN_markers_violin <- VlnPlot(DNs.RNA.seu, features = paste0(signature.names[1:2], c(1:2)), group.by = "RNA_snn_res.4.4", ncol = 1)

(SN_markers_umap | SN_markers_violin) + plot_layout(widths = c(1,4))

#Sn and VTA classification
DNs.RNA.seu$region <- ifelse(DNs.RNA.seu$SN_vedran1 > 0, "SN", "VTA")
DNs.RNA.seu.VTA <- subset(x = DNs.RNA.seu, subset = region == "VTA")

# ========== 4. Reclustering of VTA-only Cells ==========
DNs.RNA.seu.VTA <- FindNeighbors(DNs.RNA.seu.VTA, dims = 1:10)
DNs.RNA.seu.VTA <- FindClusters(DNs.RNA.seu.VTA, resolution = resolutions)
DNs.RNA.seu.VTA <- RunUMAP(DNs.RNA.seu.VTA, dims = 1:10)

#4.1 Sample per cluster
DNs.RNA.seu.VTA_metadata <- DNs.RNA.seu.VTA@meta.data 
sample_perCluster <- DNs.RNA.seu.VTA_metadata %>% 
  dplyr::group_by(RNA_snn_res.4.4, orig.ident) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n) * 100) %>% group_by(RNA_snn_res.4.4) %>% 
  dplyr::mutate(label = glue("{sum(n)}"))
sample_perCluster <- as.data.frame(sample_perCluster)

sample_perCluster$orig.ident <- factor(sample_perCluster$orig.ident, levels = samples)
samplePerCluster_plot = ggplot() +
  geom_bar(data = sample_perCluster, aes(x = RNA_snn_res.4.4, y = freq, fill = orig.ident), stat="identity") +
  scale_fill_manual(values = sample_colors) + 
  theme_classic() + theme(legend.position = "right") + 
  labs(title = "Contribution of each sample to each DN-cluster", x= "", y= "%cells") +
  geom_text(aes(RNA_snn_res.4.4, 100 + 2, label = label, fill = NULL), size = 4,  data = sample_perCluster) +
  theme(legend.position="bottom") +
  guides(fill = guide_legend(nrow = 3))
samplePerCluster_plot + no_legend()

#4.2. Save plot
# dev.size()
# width_original = 9.952381
# height_original= 9.247619
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter1_VTA_DNs_samples_LEGEND"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# samplePerCluster_plot
# dev.off()


# ========== 5. Literature-based subtype classification ==========
dopalit <- c("Th", "Aldh1a1", "Sox6", "Slc17a6", "Otx2", "Slc32a1", "Vip", "Gad2", "Cck", "Calb1")
FeaturePlot(DNs.RNA.seu.VTA, features = dopalit, ncol = 5)
meta_dopa <- DNs.RNA.seu.VTA@meta.data
meta_dopa_mark_mat = meta_dopa %>% 
  cbind(
    GetAssayData(DNs.RNA.seu.VTA, "data")[dopalit,] %>%
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

DNs.RNA.seu.VTA@meta.data <- meta_dopa_mark_mat
#How many cells:
table(DNs.RNA.seu.VTA$dopa_literature_type)
table(is.na(DNs.RNA.seu.VTA$dopa_literature_type))

#Make the plots:
literature.subtypes.plot <- DimPlot(DNs.RNA.seu.VTA, group.by = 'dopa_literature_type', cols = literature.colors) + no_legend()
Idents(DNs.RNA.seu.VTA) <- "dopa_literature_type"
literature.plot.split <- DimPlot(DNs.RNA.seu.VTA, split.by = "dopa_literature_type", cols = literature.colors, ncol = 4)

#Barplot: 
dopalit_perCluster <- meta_dopa_mark_mat %>% 
  dplyr::group_by(RNA_snn_res.4.4, dopa_literature_type) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n) * 100) %>% group_by(RNA_snn_res.4.4) %>% 
  dplyr::mutate(label = glue("{sum(n)}"))
dopalit_perCluster <- as.data.frame(dopalit_perCluster)

dopalit_perCluster_plot = ggplot() +
  geom_bar(data = dopalit_perCluster, aes(x = RNA_snn_res.4.4, y = freq, fill = dopa_literature_type), stat="identity") +
  scale_fill_manual(values = literature.colors) + 
  theme_classic() + theme(legend.position = "right") + 
  labs(title = "Contribution of DN-subtype to each DN-cluster", x= "", y= "%cells") +
  geom_text(aes(RNA_snn_res.4.4, 100 + 2, label = label, fill = NULL), size = 3,  data = dopalit_perCluster) +
  theme(legend.position="none")
dopalit_perCluster_plot

layout <- "
AABBBB
"
DNs.RNA.seu.VTA.plots <- literature.subtypes.plot + no_legend() + 
  dopalit_perCluster_plot + no_legend() + plot_layout(design = layout)
DNs.RNA.seu.VTA.plots


# ========== 6. Gene Signature Visualization ==========
IEGs_umap<- FeaturePlot(DNs.RNA.seu.VTA, features = paste0(signature.names[4:5], c(4:5)), label = FALSE, repel = TRUE, ncol = 1)
IEGs_violin <- VlnPlot(DNs.RNA.seu.VTA, features = paste0(signature.names[4:5], c(4:5)), group.by = "RNA_snn_res.4.4", ncol = 1)
IEGs_plot <- IEGs_umap | IEGs_violin

#5.2. IEG positive cluster
cluster0_umap <- FeaturePlot(DNs.RNA.seu.VTA, features = paste0(signature.names[6], c(6)), label = FALSE, repel = TRUE, ncol = 1)
cluster0_violin <- VlnPlot(DNs.RNA.seu.VTA, features = paste0(signature.names[6], c(6)), group.by = "RNA_snn_res.4.4", ncol = 1)
Idents(DNs.RNA.seu.VTA) <- "RNA_snn_res.4.4"

layout <- "
AACCC
BBCCC
"
cluster0 <- DimPlot(DNs.RNA.seu.VTA) + cluster0_umap + cluster0_violin + plot_layout(design = layout)

#5.3. Tissowsky plots
###############################
tissowsky_umap <- FeaturePlot(DNs.RNA.seu.VTA, features = paste0(signature.names[9], c(9)), label = FALSE, repel = TRUE, ncol = 1)
tissowsky_violin <- VlnPlot(DNs.RNA.seu.VTA, features = paste0(signature.names[9], c(9)), group.by = "RNA_snn_res.2.4", ncol = 1)
tissowsky_plot <- tissowsky_umap| tissowsky_violin + plot_layout(widths = c(1,3))


#6 Plot other genes: Pcdh9, Rbfox1, Vip 
#######################################
important_genes <- c("Rbfox1", "Vip", "Pcdh9", "Cartpt", "Nr4a1", "Bdnf")
important_genes_plot <- FeaturePlot(DNs.RNA.seu.VTA, features = important_genes, label = FALSE, repel = TRUE, ncol = 3)


# ========== 7. Save Figures and Metadata ==========
#load("231002.RNAonly.excluding_putative_SN_cells.rds")
#write_tsv(meta_dopa_mark_mat %>% rownames_to_column("cellNames"), "231002_VTA_RNAonly_metadata_1570cells.tsv")

#Sample contribution per cluster
DNs.RNA.seu.VTA_metadata <- DNs.RNA.seu.VTA@meta.data 
sample_perCluster <- DNs.RNA.seu.VTA_metadata %>% 
  dplyr::group_by(RNA_snn_res.4.4, orig.ident) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n) * 100) %>% group_by(RNA_snn_res.4.4) %>% 
  dplyr::mutate(label = glue("{sum(n)}"))
sample_perCluster <- as.data.frame(sample_perCluster)

sample_perCluster$orig.ident <- factor(sample_perCluster$orig.ident, levels = samples)
samplePerCluster_plot = ggplot() +
  geom_bar(data = sample_perCluster, aes(x = RNA_snn_res.4.4, y = freq, fill = orig.ident), stat="identity") +
  scale_fill_manual(values = sample_colors) + 
  theme_classic() + theme(legend.position = "right") + 
  labs(title = "Contribution of each sample to each DN-cluster", x= "", y= "%cells") +
  geom_text(aes(RNA_snn_res.4.4, 100 + 2, label = label, fill = NULL), size = 4,  data = sample_perCluster) +
  no_legend()

# Print the plot
# print(samplePerCluster_plot)
# dev.size()
# width_original = 11.156250
# height_original= 4.52
# 
# pdf_dir <- "/fast/AG_Pombo/luna/2025_pdf_files/"
# plot_name <- "chapter1_VTA_DNs_barPlot"
# file_name <- glue("{pdf_dir}/{plot_name}.pdf" )
# 
# pdf(file_name, width = width_original, height =height_original)
# samplePerCluster_plot
# dev.off()