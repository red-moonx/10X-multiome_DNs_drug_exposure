# ===========================================
# Script Title: 1.1 Integration 18 samples through Conos
# Author: Luna Zea Redondo
# Date: 2026-05-02
# Description:
#   This script performs the integration of the 18 different VTA samples at the  
#   RNA level and identifies and stores the DN specific cluster (13). 
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
library(conos)


# ========== Set wd and load data ==========

setwd("/fast/AG_Pombo/luna/2026_rebuttal/1_conos")

#Seurat list of 18 VTA samples
seu.list <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/2_cellBender/1_QCanalysis/1_HVGmethods_comparison/230410.combined.seu.obj.list.rds")


# --------------------------------------------------------
# 1. Compute consensus of highly variable genes (HVGs)
# --------------------------------------------------------
sample_list <- lapply(seu.list, function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})


# -----------------------------------
# 2. Run integration and embedding
# -----------------------------------

con <- Conos$new(sample_list, n.cores=4)
con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)

con$embedGraph(method="UMAP", min.dist=0.01, spread=15, min.prob.lower=1e-3)
#load("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/260502_backup.rds")


# -------------------------------------------------------
# 3. Run clustering and identify the DN enriched cluster
# -------------------------------------------------------
con$findCommunities(method=leiden.community, resolution=1)
integration_Th <- con$plotGraph(gene='Th', title='Th expression',  show.legend=TRUE)

#Cluster 13 seems to be the DN specific cluster
cells_in_13 <- which(con$clusters$leiden$groups == "13")
integration_Th <- con$plotGraph(gene='Th', title='Th expression',  show.legend=TRUE)
integration_Slc6a3 <- con$plotGraph(gene='Slc6a3', title='Slc6a3 expression',  show.legend=TRUE)
integration_Nr4a2 <- con$plotGraph(gene='Nr4a2', title='Nr4a2 expression',  show.legend=TRUE)
integration_Slc18a2 <- con$plotGraph(gene='Slc18a2', title='Slc18a2 expression',  show.legend=TRUE)
integration_Snca <- con$plotGraph(gene='Snca', title='Snca expression',  show.legend=TRUE)
integration_Foxa2 <- con$plotGraph(gene='Foxa2', title='Foxa2 expression',  show.legend=TRUE)
integration_Lmx1b <- con$plotGraph(gene='Lmx1b', title='Lmx1b expression',  show.legend=TRUE)
integration_Kcnj6 <- con$plotGraph(gene='Kcnj6', title='Kcnj6 expression',  show.legend=TRUE)
integration_Ddc <- con$plotGraph(gene='Ddc', title='Ddc expression',  show.legend=TRUE)
integration_Pbx1 <- con$plotGraph(gene='Pbx1', title='Pbx1 expression',  show.legend=TRUE)
integration_Lmo3 <- con$plotGraph(gene='Lmo3', title='Lmo3 expression',  show.legend=TRUE)
integration_Calb1 <- con$plotGraph(gene='Calb1', title='Calb1 expression',  show.legend=TRUE)

integration_Th + integration_Slc6a3 + integration_Nr4a2 + integration_Slc18a2 +
  integration_Snca + integration_Foxa2 + integration_Lmx1b + integration_Kcnj6 + 
  integration_Ddc + integration_Pbx1 + integration_Lmo3 + integration_Calb1 + plot_layout(ncol = 4)


# I am just getting the clusters and the embedding
con.embeddings <- con[["embeddings"]][["UMAP"]]
con.clusters <- con[["clusters"]][["leiden"]][["groups"]]
cell_names_13 <- names(con$clusters$leiden$groups)[cells_in_13] # 2761 DN cells



# ---------------------
# 4. Convert to Seurat 
# ---------------------
seurat_integrated <- merge(
  x = con$samples[[1]], 
  y = con$samples[-1], 
  add.cell.ids = names(con$samples)
)

# Clean the double prefix (e.g. "m30_cocaine_R1_m30_cocaine_R1#..." naming issue)
current_seurat_names <- colnames(seurat_integrated)

clean_names <- sapply(current_seurat_names, function(x) {
  parts <- strsplit(x, "#")[[1]]
  barcode <- parts[2]
  prefix_double <- parts[1]
  
  # Split prefix by underscores and take the second half
  prefix_parts <- strsplit(prefix_double, "_")[[1]]
  half_len <- length(prefix_parts) / 2
  single_prefix <- paste(prefix_parts[(half_len + 1):length(prefix_parts)], collapse = "_")
  
  return(paste0(single_prefix, "#", barcode))
})

# Apply the names using Seurat v5 syntax
seurat_integrated <- RenameCells(seurat_integrated, new.names = unname(clean_names))

# sanity check: are all cells in con found in serat object? 
common_cells <- intersect(rownames(con$embeddings$UMAP), colnames(seurat_integrated)) #yes: 100431 cells
seurat_integrated <- subset(seurat_integrated, cells = common_cells)
umap_coords <- con$embeddings$UMAP[common_cells, ]


# Inject Conos results: 
# UMAP coordinates:
seurat_integrated[["umap"]] <- CreateDimReducObject(
  embeddings = umap_coords, 
  key = "UMAP_", 
  assay = DefaultAssay(seurat_integrated)
)

# Leiden clusters calculated in Conos
seurat_integrated$conos_clusters <- con$clusters$leiden$groups[common_cells]
Idents(seurat_integrated) <- "conos_clusters"

# -----------------------
# 5. Extract cluster 13 
# -----------------------
Idents(seurat_integrated) <- "conos_clusters"
cluster13 <- subset(seurat_integrated, idents = "13") #2761 putative DNs


# -------------------
# 5. Data Archiving
# -------------------
# Save the entire session
# save.image("/fast/AG_Pombo/luna/2026_rebuttal/1_conos/260502.conos.RNA.backup.session.rds")

# Save specific files:
# saveRDS(seurat_integrated, "/fast/AG_Pombo/luna/2026_rebuttal/1_conos/seu.entireVTA-RNA.rds")
# saveRDS(cluster13, "/fast/AG_Pombo/luna/2026_rebuttal/1_conos/seu.cluster13_putativeDNs.rds")

# write.csv(cell_names_13, "/fast/AG_Pombo/luna/2026_rebuttal/1_conos/cellNames_cluster13.csv")



# ------------------------
# 6. Plots for manuscript:
# ------------------------

# =========================================================
# Marker lists
# =========================================================

DN_module_score_genes <- c("Th", "Slc6a3", "Nr4a2", "Slc18a2", "Snca", "Foxa2", "Lmx1b", "Kcnj6")
OPC_markers <- c("Pdgfra", "Cspg4", "Ptprz1", "Sox6", "Olig1", "Olig2", "Sox10", "Tnr", "Fgfr3", "Gpr17")
endothelial_markers <- c("Cldn5", "Flt1", "Kdr", "Pecam1", "Emcn", "Klf2", "Esam", "Robo4", "Ly6c1", "Vwf")
microglia_markers <- c("P2ry12", "Tmem119", "Cx3cr1", "Csf1r", "Hexb", "Fcrls", "Sall1", "Aif1", "Tyrobp", "Ctss")
astrocyte_markers <- c("Aldh1l1", "Aqp4", "Slc1a3", "Gfap", "Glul", "Fgfr3", "S100b", "Clu", "Hepacam", "Agt")
oligodendrocyte_markers <- c("Mbp", "Plp1", "Mog", "Mobp", "Mag", "Cnp", "Opalin", "Ermn", "Ugt8a", "Aspa")
pan_neuronal_markers <- c("Snap25", "Rbfox3", "Syt1", "Tubb3", "Elavl2", "Elavl3", "Elavl4", "Map2", "Dlg4", "Stmn2")

# =========================================================
# Add module scores
# =========================================================

marker_lists <- list(
  DN_Score = DN_module_score_genes,
  OPC_Score = OPC_markers,
  Endothelial_Score = endothelial_markers,
  Microglia_Score = microglia_markers,
  Astrocyte_Score = astrocyte_markers,
  Oligodendrocyte_Score = oligodendrocyte_markers,
  PanNeuronal_Score = pan_neuronal_markers
)

for(i in names(marker_lists)){
  seurat_integrated <- AddModuleScore(
    seurat_integrated,
    features = list(marker_lists[[i]]),
    name = i
  )
}

# =========================================================
# Plots
# =========================================================

clean_theme <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  legend.position = "none"
)

clean_theme <- theme_void() + theme(legend.position = "none")

p1 <- DimPlot(
  seurat_integrated,
  reduction = "umap", pt.size = 0.05,
  raster = FALSE
) +
  ggtitle("Clusters") +
  clean_theme


p2 <- FeaturePlot(
  seurat_integrated,
  features = "DN_Score1",
  cols = c("lightgrey", "red"), pt.size = 0.05,
  raster = FALSE
) +
  ggtitle("DN Module Score") +
  clean_theme

pop_plots <- wrap_plots(
  FeaturePlot(seurat_integrated, features = "OPC_Score1", cols = c("lightgrey", "red"), pt.size = 0.05, raster = FALSE) + ggtitle("OPCs") + clean_theme,
  FeaturePlot(seurat_integrated, features = "Endothelial_Score1", cols = c("lightgrey", "red"), pt.size = 0.05, raster = FALSE) + ggtitle("Endothelial") + clean_theme,
  FeaturePlot(seurat_integrated, features = "Microglia_Score1", cols = c("lightgrey", "red"), pt.size = 0.05, raster = FALSE) + ggtitle("Microglia") + clean_theme,
  FeaturePlot(seurat_integrated, features = "Astrocyte_Score1", cols = c("lightgrey", "red"), pt.size = 0.05, raster = FALSE) + ggtitle("Astrocytes") + clean_theme,
  FeaturePlot(seurat_integrated, features = "Oligodendrocyte_Score1", cols = c("lightgrey", "red"), pt.size = 0.05, raster = FALSE) + ggtitle("Oligodendrocytes") + clean_theme,
  FeaturePlot(seurat_integrated, features = "PanNeuronal_Score1", cols = c("lightgrey", "red"), pt.size = 0.05, raster = FALSE) + ggtitle("Pan-neuronal") + clean_theme,
  ncol = 3
)

clusters_and_DNs <- (p1 + p2)


# Save plot
# dev.size()
# width_original = 8.968750 
# height_original= 4.864583
# plot_name <- "chapter1_clusters_and_DNmodulescore.pdf"
# pdf(plot_name, width = width_original, height =height_original)
# clusters_and_DNs
# dev.off()


# Save plot
# dev.size()
# width_original = 8.968750 
# height_original= 4.864583
# plot_name <- "chapter1_population_plots.pdf"
# pdf(plot_name, width = width_original, height =height_original)
# pop_plots
# dev.off()


# =========================================================
# Save clusters + DN plot
# =========================================================

# tiff(
#   "chapter1_clusters_and_DNmodulescore.tiff",
#   width = 8.968750,
#   height = 4.864583,
#   units = "in",
#   res = 300,
#   compression = "lzw"
# )
# 
# print(clusters_and_DNs)
# 
# dev.off()

# =========================================================
# Save population plots
# =========================================================

# tiff(
#   "chapter1_population_plots.tiff",
#   width = 8.968750,
#   height = 4.864583,
#   units = "in",
#   res = 300,
#   compression = "lzw"
# )
# 
# print(pop_plots)
# 
# dev.off()
