# ===========================================
# Script Title: 1.4 VTA-specific DEA with Muscat
# Author: Luna Zea
# Date: 2025-06-05 (for github)
# Description:
#   This script performs differential gene expression analysis (DEA) using Muscat on 
#   VTA-specific DN cells. We run each vs. saline and each vs. preceding comparisons. 
#   It includes DEG calling, result filtering, and visualization through barplots,
#   upset plots, volcano/MA plots, and heatmap.
# ===========================================

# ========== Set Environment ==========
rm(list = ls(all.names = TRUE))
gc()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

# ========== Load Required Libraries ==========
library(muscat)
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(glue)
library(ggplot2)
library(patchwork)
library(reshape2)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(UpSetR)

# ========== Set Directory ==========
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly"
setwd(dir)

# ========== 1. Pre-processing ==========

load("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002.RNAonly.excluding_putative_SN_cells.rds")
DNs.RNA.seu #2411 cells; "region" separates VTA from SN cells

#0.1) Select time points and treatments of interest: we are excluding h1_saline
Idents(DNs.RNA.seu) <- "orig.ident"
samples_totest <- c("h4_saline_R1", "h8_saline_R1", "h24_saline_R1", "d14_saline_R1", #salines
                    "h1_cocaine_R1", "h1_cocaine_R2", "h4_cocaine_R1", "h4_cocaine_R2", "h8_cocaine_R1", "h8_cocaine_R2", #ETP
                    "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3") #LTP

sample_totest_colors <- sample_colors[names(sample_colors) %in% samples_totest]
sample_totest_colors <- sample_totest_colors[samples_totest]
DNs.RNA.seu.2reps <- subset(x = DNs.RNA.seu, idents = samples_totest)

#0.2) Convert the seurat object to sce: https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
DNs.RNA.sce <- as.SingleCellExperiment(DNs.RNA.seu.2reps)

#0.3) Add subpopulations:
DNs.RNA.sce$Analysis1 <- "All_DNs"
DNs.RNA.sce$Analysis2 <- DNs.RNA.sce$region


#0.4) Set subpopulations and colors
groups_colors <- c('All_DNs'="#989C30", 'VTA'="#9C3848", 'SN'="#47A8BD")


#0.5) Back up:
sce <- DNs.RNA.sce # Keep DNs.RNA.sce as unaltered copy

#Just for me, for easy access:
coldata.dataframe <- colData(sce)
table(coldata.dataframe$orig.ident)
table(coldata.dataframe$Analysis1)
table(coldata.dataframe$Analysis2)

sample.info <- as.data.frame(coldata.dataframe) %>%
  select(ident, simpleIdent, timepoint, treatment, replicate, sample_prep, FACS, GEXlibrary, Sequencing) %>% 
  distinct()
all_analysis <- c("Analysis1", "Analysis2")

saline_contrasts <-c("h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline", "h24_cocaine-saline", "d14_cocaine-saline")
time_contrasts <- c("h1_cocaine-saline", "h4_cocaine-h1_cocaine", "h8_cocaine-h4_cocaine", "h24_cocaine-h8_cocaine", "d14_cocaine-h24_cocaine")

#Analysis 1 and 2 refer to different groupings: 1 (all DNs, without VTA or SN distinction) and 2 (VTA and SN separated)
#Both could be run, we focused on Analysis 2, aiming for max specificity for the VTA DNs

for (analysis in all_analysis) {
  analysis <- all_analysis[2]
  
  # ========== 2. Initial set up, per analysis, and QC plots ==========
  print(analysis)
  sce <- DNs.RNA.sce # Keep DNs.RNA.sce as unaltered copy
  dir <- glue("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002_muscat_VTAvsSN_all_subtypes/{analysis}/")
  dir.create(dir)
  setwd(dir)
  
  #Plot samples info
  samples.coldata <- as.data.frame(coldata.dataframe) %>% dplyr::select(ident, simpleIdent, analysis)
  samples.coldata <- melt(dcast(samples.coldata, ident ~ samples.coldata[, analysis], value.var = analysis)) %>% filter(value > 0) 
  samples.coldata <- left_join(samples.coldata, sample.info, by = "ident")
  
  samples.coldata$ident <- factor(samples.coldata$ident, levels = samples_totest)
  samples.coldata$simpleIdent <- factor(samples.coldata$simpleIdent, levels = c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine"))
  
  samples.coldata.plot <- ggplot(samples.coldata) + 
    geom_bar(aes(y = value, x = simpleIdent, fill = ident), stat="identity") +
    geom_text(data=samples.coldata, aes(x = simpleIdent, y = as.numeric(value), label = value, size=4, fill = ident), position = position_stack(vjust = 0.5), show.legend = F) + 
    scale_fill_manual(values= sample_totest_colors) + theme_classic() + ggtitle("Number of cells per sample") + 
    theme(legend.position = "right", 
          axis.title.x=element_blank(), axis.title.y=element_blank()) +
    guides(fill=guide_legend(title="Sample")) +
    facet_wrap(~variable)
  samples.coldata.plot
  ggsave(filename = glue("{dir}/{analysis}_samples.coldata.plot.png"), samples.coldata.plot, units = "px", device = "png",width=5000,height=2400,dpi = 300)
  
  #Reduce number of genes and normalize
  sce <- sce[rowSums(counts(sce) > 0) > 0, ]
  dim(sce)
  nondetected <- (nrow(DNs.RNA.sce) - nrow(sce))
  # calculate per-cell quality control (QC) metrics
  qc <- perCellQCMetrics(sce)
  
  # remove lowly expressed genes
  #sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
  dim(sce)
  lowly_expressed <- (nrow(DNs.RNA.sce) - nondetected - nrow(sce))
  tested <- nrow(sce)
  # compute sum-factors & normalize
  sce <- computeLibraryFactors(sce)
  sce <- logNormCounts(sce)
  
  # plot genes info
  #Genes:
  gene_categories <- c("Non-detected", "Lowly-expressed", "Being tested")
  gene_numbers <- c(nondetected, lowly_expressed, tested)
  gene_info <- data.frame(gene_categories, gene_numbers)
  gene_info$gene_categories <- factor(gene_info$gene_categories, levels = gene_categories)
  genes_colors <- c("#1B1B1E", "#1B4B5D", "#F09E05")
  names(genes_colors) <- gene_categories
  
  genes_info.plot <- ggplot(gene_info) + 
    geom_bar(aes(y = gene_numbers, x = gene_categories, fill = gene_categories), stat="identity") +
    geom_text(data=gene_info, aes(x = gene_categories, y = as.numeric(gene_numbers), label = gene_numbers, size=5, fill = gene_categories), color = "white", position = position_stack(vjust = 0.5), show.legend = F) + 
    scale_fill_manual(values= genes_colors) + theme_classic() +
    ggtitle("Number of genes per category") +theme(legend.position = "right", 
                                                   axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
    guides(fill=guide_legend(title="Gene categories"))
  genes_info.plot
  ggsave(filename = glue("{dir}/{analysis}_genes_info.plot.png"), genes_info.plot, units = "px", device = "png",width=5000,height=2400,dpi = 300)
  
  
  
  # ========== 1. Run Muscat DEG Analysis ==========
  
  #Data preparation 
  (sce <- prepSCE(sce, 
                  kid = analysis, # subpopulation assignments
                  gid = "simpleIdent",  # group IDs (ctrl/stim)
                  sid = "orig.ident",   # sample IDs (ctrl/stim.1234)
                  drop = TRUE))  # drop all other colData columns
  
  nk <- length(kids <- levels(sce$cluster_id))
  ns <- length(sids <- levels(sce$sample_id))
  names(kids) <- kids; names(sids) <- sids
  
  
  #Cluster-samples sizes
  # Number of cells per cluster-sample
  t(table(sce$cluster_id, sce$sample_id))
  
  
  #Dimensional reduction
  #compute UMAP using 1st 20 PCs
  sce <- runUMAP(sce, pca = 30)
  # wrapper to prettify reduced dimension plots
  .plot_dr <- function(sce, dr, col)
    plotReducedDim(sce, dimred = dr, colour_by = col) +
    guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_minimal() + theme(aspect.ratio = 1)
  
  # downsample to max. 100 cells per cluster
  
  cs_by_k <- split(colnames(sce), sce$cluster_id)
  cs100 <- unlist(sapply(cs_by_k, function(u) 
    sample(u, min(length(u), 100))))
  
  # plot t-SNE & UMAP colored by cluster & group ID
  for (dr in c("UMAP"))
    for (col in c("cluster_id", "group_id"))
      .plot_dr(sce[, cs100], dr, col)
  DR <-.plot_dr(sce[, cs100], dr, col)
  
  
  #Aggregation of single-cell to pseudobulk data
  pb <- aggregateData(sce,
                      assay = "counts", fun = "sum",
                      by = c("cluster_id", "sample_id"))
  #one sheet per subpopulation
  assayNames(pb)
  
  #pseudobulks for 1st subpopulation
  t(head(assay(pb)))
  
  #Pseudo-bulk level MDS plot
  (pb_mds <- pbMDS(pb))
  three_groups_colors <- groups_colors
  #names(three_groups_colors) <- assayNames(pb)
  # use very distinctive shaping of groups & change cluster colors
  pb_mds <- pb_mds + 
    scale_shape_manual(values = c(15, 17, 18, 19, 20, 4)) +
    scale_color_manual(values = three_groups_colors)
  # change point size & alpha
  pb_mds$layers[[1]]$aes_params$size <- 5
  pb_mds$layers[[1]]$aes_params$alpha <- 1
  pb_mds
  

  #3) Differential state (DS) analysis)

  ###3.1) Aggregation of single-cell to pseudobulk data
  #####################################################
  pb <- aggregateData(sce,
                      assay = "counts", fun = "sum",
                      by = c("cluster_id", "sample_id"))
  # one sheet per subpopulation
  assayNames(pb)
  
  # pseudobulks for 1st subpopulation
  t(head(assay(pb)))
  
  ##Pseudo-bulk level MDS plot
  ####################################
  (pb_mds <- pbMDS(pb))
  pb_colors <- groups_colors[names(pb@assays@data@listData)]
  
  # use very distinctive shaping of groups & change cluster colors
  pb_mds <- pb_mds + 
    scale_shape_manual(values = c(15, 17, 18, 19, 20, 4)) +
    scale_color_manual(values = pb_colors)
  # change point size & alpha
  pb_mds$layers[[1]]$aes_params$size <- 5
  pb_mds$layers[[1]]$aes_params$alpha <- 1
  pb_mds 
  
  muscat_plots <- DR + pb_mds
  ggsave(filename = glue("{dir}/muscat_plots.png"), muscat_plots, units = "px", device = "png",width=5000,height=3000,dpi = 300)
  
  #Sample-level analysis: Pseudobulk 
  ##########################################
  # construct design & contrast matrix
  
  ei <- metadata(sce)$experiment_info
  mm <- model.matrix(~ 0 + ei$group_id)
  dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
  
  condition_colors <- c("black", "#9F2365", "#617641", "#C48208", "#326186", "#AE430A", "#564686")
  condition_names <- c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine")
  names(condition_colors) <- condition_names
  
  #Establish all possible contrasts:
  all_contrasts <- c("h1_cocaine-saline",
                     "h4_cocaine-saline",
                     "h8_cocaine-saline",
                     "h24_cocaine-saline",
                     "d14_cocaine-saline", 
                     "h4_cocaine-h1_cocaine",
                     "h8_cocaine-h4_cocaine", 
                     "h24_cocaine-h8_cocaine", 
                     "d14_cocaine-h24_cocaine")
  
  
  #I cannot work with m30 (one replicate)
  dir <- glue("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002_muscat_VTAvsSN_all_subtypes//{analysis}/")
  setwd(dir)
  dir.create("deg_results")
  dir.create("volcano_and_MA")
  
  DEG_complete_results <-data.frame()
  for (contrast in all_contrasts) {
    message(glue("Analyzing: {contrast}"))
    query=str_split(contrast, "-")[[1]][1]
    control=str_split(contrast, "-")[[1]][2]
    contrast_label = contrast
    
    message(glue("Query is: {query} and control is: {control}"))
    
    #Establish current contrast:
    contrast <- makeContrasts(contrast, levels = mm)
    
    # run DS analysis and save table
    res <- pbDS(pb, method = "edgeR", design = mm, contrast = contrast, filter = "none")
    DEG_contrast_results <- resDS(sce, res, bind = "row", cpm = TRUE, frq = TRUE) %>% mutate(
      query = query,
      control = control, 
      significant = ifelse(p_val < 0.05 & logFC > 0.5, "Upregulated (pval < 0.05 and logFC > 0.5)", 
                           ifelse(p_val < 0.05 & logFC < -0.5, "Downregulated (pval < 0.05 and logFC < -0.5)", "No significant"))) 
    
    DEG_contrast_results$percent5 <- ifelse(DEG_contrast_results[[glue("{query}.frq")]] >= 0.05 |  DEG_contrast_results[[glue("{query}.frq")]] >= 0.05, "Yes",  "No")
    DEG_contrast_results$significant_and_5per = ifelse(DEG_contrast_results$significant != "No significant" & DEG_contrast_results$percent5 == "Yes", DEG_contrast_results$significant, "No significant")
    DEG_complete_results <- rbind(DEG_complete_results,DEG_contrast_results)
    
  }
  
  write_tsv(DEG_complete_results, glue("deg_results/231002_DEG_complete_results_{analysis}.tsv"))
  
  
  # ========== 3. DEG Summary Plot (per cluster, per contrast) ==========
  # How many genes per contrast are DEG present in 10% of the cells
  # Plot samples info
  DEG.info <- DEG_complete_results %>% 
    dplyr::filter(significant != "No significant") %>% 
    dplyr::select(contrast, percent5,cluster_id)  
  
  DEG.info$contrast <- factor(DEG.info$contrast, levels = all_contrasts)
  DEG.info$percent5 <- factor(DEG.info$percent5, levels = c("Yes", "No"))
  percent5_colors <- c(pb_colors, "#BABFD1")
  names(percent5_colors) <- c(names(pb_colors), "No")
  
  percent5_colors <- percent5_colors[c("No",names(pb_colors))]
  DEG.info.melt <- melt(dcast(DEG.info, cluster_id + percent5 ~ contrast))
  
  DEG.info.melt$percent5 <- as.character(DEG.info.melt$percent5)
  DEG.info.melt$percent5 <- ifelse(DEG.info.melt$percent5 == "Yes", DEG.info.melt$cluster_id, DEG.info.melt$percent5) 
  
  #Temporal code (221011)
  DEG.info.melt.saline <- DEG.info.melt %>% filter(variable %in% saline_contrasts)
  DEG.info.melt.previous <- DEG.info.melt %>% filter(variable %in% time_contrasts)
  
  DEG.info.plot <- ggplot(DEG.info.melt) + 
    geom_bar(aes(y = value, x = variable, fill = percent5), stat="identity") +
    geom_text(data=DEG.info.melt, aes(x = variable, y = as.numeric(value), label = value, size=16, fill = percent5), col= "white", position = position_stack(vjust = 0.5), show.legend = F) + 
    scale_fill_manual(values= percent5_colors) + theme_classic() + ggtitle(glue("Number of DEGs per contrast:")) +
    theme(legend.position = "bottom", 
          axis.title.x=element_blank(), axis.title.y=element_blank()) +
    guides(fill=guide_legend(title="Detected in >5% of cells")) +
    facet_wrap(~cluster_id, ncol = 1)
  
  DEG.info.plot
  ggsave(filename = glue("{dir}/deg_results/DEG.info.plot.png"), DEG.info.plot, units = "px", device = "png",width=2100,height=1500,dpi = 300)
  
  a <- DEG.info.plot
  
  #4.2) Upset plot:
  upset_list <- list()
  for (subpop in kids) {
    print(subpop)
    clusterid_genes <- DEG_complete_results %>% filter(cluster_id == subpop, 
                                                       significant_and_5per != "No significant") %>% select(gene) %>% pull()
    upset_list[[subpop]] <- clusterid_genes
  }
  
  if(length(upset_list) >1 ) {
    upset_plot<- upset(fromList(upset_list),
                       order.by = "freq", text.scale = 2, point.size = 5, 
                       sets.bar.color= c("#47A8BD", "#9C3848")) 
    
    png(glue("{dir}/deg_results/upsetplot.png"))
    upset_plot
    dev.off()
    
  } else {
    print("Only one group tested")
  }
  
  # ========== 4. Volcano & MA Plots ==========
  #For violin plots
  seurat_clusters_selected <- levels(DNs.RNA.seu$seurat_clusters) 
  seu.temporal <- DNs.RNA.seu[, DNs.RNA.seu$seurat_clusters %in% seurat_clusters_selected]
  Idents(seu.temporal) <- "simpleIdent"
  Idents(seu.temporal) <- factor(Idents(seu.temporal), levels= c("saline", "m30_cocaine", "h1_cocaine", "h4_cocaine", "h8_cocaine", "h24_cocaine", "d14_cocaine"))
  
  for (comparison in all_contrasts) {
    print(glue("Generating volcano and MA plot for contrast: {comparison}"))
    plot.data <- DEG_complete_results %>% filter(contrast == comparison)
    plot.data.label <- plot.data %>%
      filter(significant_and_5per != "No significant", 
             cluster_id == "VTA") %>% 
      slice_max(abs(logFC), n=50)
    query=str_split(comparison, "-")[[1]][1]
    control=str_split(comparison, "-")[[1]][2]
    
    volcano_colors <- c(condition_colors[query], condition_colors[control], "gray")
    names(volcano_colors) <- c("Upregulated (pval < 0.05 and logFC > 0.5)", "Downregulated (pval < 0.05 and logFC < -0.5)", "No significant")
    volcano_plot <- ggplot(plot.data, aes(x = logFC, y = -log10(p_val))) +
      geom_point(aes(color = significant_and_5per)) +
      scale_color_manual(values = volcano_colors) +
      theme_classic(base_size = 12) + theme(legend.position = "bottom") +
      geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") + 
      geom_vline(xintercept = -0.5, colour="#990000", linetype="dashed") + 
      geom_vline(xintercept = 0.5, colour="#990000", linetype="dashed") + 
      geom_text_repel(
        data = plot.data.label, aes(label = gene),
        size = 3, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"), max.overlaps = 200) +
      labs(title = glue("Differential Expression Analysis: {query} vs {control}: top 50 VTA genes"),
           subtitle = "pval<0.05; |logFC|>0.5; detected in >5% cells", 
           x= "Log (FC)", 
           y="-Log10 (p value)") +
      theme(legend.position = "bottom", legend.title = element_blank()) +
      guides(fill=guide_legend(title="Detected in >5% of cells")) +
      facet_wrap(~ cluster_id)
    
    ggsave(filename = glue("volcano_and_MA/{comparison}.volcano.png"), volcano_plot, units = "px", device = "png",width=5000,height=3000,dpi = 300)
    
    MA_plot <- ggplot(plot.data, aes(x = logCPM, y = logFC)) +
      geom_point(aes(color = significant_and_5per)) +
      scale_color_manual(values = volcano_colors) +
      theme_classic(base_size = 12) + theme(legend.position = "bottom") +
      geom_hline(yintercept = c(-0.5, 0.5), colour="#990000", linetype="dashed") + 
      geom_hline(yintercept = 0, colour="#990000") + 
      geom_text_repel(
        data = plot.data.label, aes(label = gene),
        size = 3, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"), max.overlaps = 200) +
      labs(title = glue("Differential Expression Analysis: {query} vs {control}: top 50 VTA genes"),
           subtitle = "pval<0.05; |logFC|>0.5; detected in >5% cells", 
           x= "Log (CPM)", 
           y="Log (FC)") +
      theme(legend.position = "bottom", legend.title = element_blank()) + 
      facet_wrap(~ cluster_id)
    ggsave(filename = glue("volcano_and_MA/{comparison}.MA.png"), MA_plot, units = "px", device = "png",width=5000,height=2000,dpi = 300)
    
  }
  
  # ========== 5. Heatmap: All DEGs per Group ==========
  saline_contrasts <-c("h1_cocaine-saline", "h4_cocaine-saline", "h8_cocaine-saline", "h24_cocaine-saline", "d14_cocaine-saline")
  time_contrasts <- c("h1_cocaine-saline", "h4_cocaine-h1_cocaine", "h8_cocaine-h4_cocaine", "h24_cocaine-h8_cocaine", "d14_cocaine-h24_cocaine")
  
  #kids_clusters <- get(glue("kids_{analysis}"))
  dir <- glue("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230926.muscat/{analysis}")
  setwd(dir)
  
  groups <- names(table(coldata.dataframe$Analysis2))
  
  #each vs saline
  dirplots <- glue("{dir}/salines")  
  dir.create(dirplots)
  for (i in 1:length(groups)) {
    clusterID = groups[i]
    heatmap_genes <- DEG_complete_results %>%
      filter(contrast %in% saline_contrasts,
             significant != "No significant",  #Changed this bit to include also < 5%
             cluster_id == clusterID) %>% select("gene") %>% pull()
    
    if (length(heatmap_genes != 0)) {
      seu.temporal <- DNs.RNA.seu[, DNs.RNA.seu$region == clusterID]
      # seu.temporal@assays[["ATAC"]] <- NULL
      # seu.temporal@assays[["SCT"]] <- NULL
      # seu.temporal@assays[["peaks"]] <- NULL
      
      #Calculate average expression and format columns
      Idents(seu.temporal) <- "orig.ident"
      seu.temporal.avgexp = as.data.frame(AverageExpression(seu.temporal, group.by = 'orig.ident')) 
      names(seu.temporal.avgexp) <- substring(names(seu.temporal.avgexp), 5)
      
      DEG_perGroup_results <- DEG_complete_results %>% filter(cluster_id == clusterID, 
                                                              contrast %in% saline_contrasts)
      
      #Extract DEGs from the big table and create metadata
      DEG_withIDs <- DEG_perGroup_results %>% 
        filter(significant != "No significant") %>% distinct(gene, .keep_all = TRUE) %>% 
        mutate(groupID = ifelse(logFC < 0, control, query)) %>% 
        arrange(factor(groupID, levels = condition_names)) 
      
      seu.DEGs.avgexp.heatmap <- seu.temporal.avgexp %>% rownames_to_column("gene") %>% 
        filter(gene %in% DEG_withIDs$gene) %>% column_to_rownames("gene")
      
      seu.DEGs.avgexp.heatmap <- seu.DEGs.avgexp.heatmap[DEG_withIDs$gene, ]
      
      
      #Prepare matrix for heatmap (in the order I want to be plotted)
      mat <- as.matrix(seu.DEGs.avgexp.heatmap)
      mat_scaled = t(scale(t(mat)))[, levels(samples.coldata$ident)]
      
      condition <- samples.coldata %>% filter(variable == clusterID) %>% arrange(ident) %>% select(simpleIdent) %>% pull()
      cluster <- factor(DEG_withIDs$groupID, levels = condition_names) #heatmap cluster
      heatmap_color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)  
      
      ha = HeatmapAnnotation(Condition = condition, show_legend = FALSE, 
                             col = list(Condition = condition_colors), show_annotation_name = FALSE)
      row_ha = rowAnnotation(Cluster = cluster, show_legend = FALSE, 
                             col = list(Cluster = condition_colors, show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"), labels = NULL))
      
      
      htmp <- Heatmap(mat_scaled, cluster_columns = FALSE, cluster_rows = FALSE, name = "Z-score", column_split =condition, row_split =cluster, row_title = NULL, show_row_names = FALSE,
                      show_column_names = FALSE, top_annotation = ha, right_annotation = row_ha, col = heatmap_color, heatmap_legend_param = list(direction = "horizontal"))
      
      
      htmp_save <- draw(htmp, heatmap_legend_side = "bottom")
      
      dev.off()
      jpeg(glue("{dirplots}/{clusterID}.salines.DEGs.heatmap.jpeg"), bg = "transparent", width=600, height=700, units = "px")
      print(htmp_save)
      dev.off()
    }
  }
  

  # each vs preceding
  dirplots <- "/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230926.muscat/Analysis2/times"
  dir.create(dirplots)
  for (i in 1:length(groups)) {
    heatmap_genes <- DEG_complete_results %>%
      filter(contrast %in% time_contrasts,
             significant != "No significant",  #Changed this bit to include also < 5%
             cluster_id == clusterID) %>% select("gene") %>% pull()
    
    if (length(heatmap_genes != 0)) {
      seu.temporal <- DNs.RNA.seu[, DNs.RNA.seu$region == clusterID]
      # seu.temporal@assays[["ATAC"]] <- NULL
      # seu.temporal@assays[["SCT"]] <- NULL
      # seu.temporal@assays[["peaks"]] <- NULL
      
      #Calculate average expression and format columns
      Idents(seu.temporal) <- "orig.ident"
      seu.temporal.avgexp = as.data.frame(AverageExpression(seu.temporal, group.by = 'orig.ident')) 
      names(seu.temporal.avgexp) <- substring(names(seu.temporal.avgexp), 5)
      
      DEG_perGroup_results <- DEG_complete_results %>% filter(cluster_id == clusterID, 
                                                              contrast %in% time_contrasts)
      
      #Extract DEGs from the big table and create metadata
      DEG_withIDs <- DEG_perGroup_results %>% 
        filter(significant != "No significant") %>% distinct(gene, .keep_all = TRUE) %>% 
        mutate(groupID = ifelse(logFC < 0, control, query)) %>% 
        arrange(factor(groupID, levels = condition_names)) 
      
      seu.DEGs.avgexp.heatmap <- seu.temporal.avgexp %>% rownames_to_column("gene") %>% 
        filter(gene %in% DEG_withIDs$gene) %>% column_to_rownames("gene")
      
      seu.DEGs.avgexp.heatmap <- seu.DEGs.avgexp.heatmap[DEG_withIDs$gene, ]
      
      
      #Prepare matrix for heatmap (in the order I want to be plotted)
      mat <- as.matrix(seu.DEGs.avgexp.heatmap)
      mat_scaled = t(scale(t(mat)))[, levels(samples.coldata$ident)]
      
      condition <- samples.coldata %>% filter(variable == clusterID) %>% arrange(ident) %>% select(simpleIdent) %>% pull()
      cluster <- factor(DEG_withIDs$groupID, levels = condition_names) #heatmap cluster
      heatmap_color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)  
      
      ha = HeatmapAnnotation(Condition = condition, show_legend = FALSE, 
                             col = list(Condition = condition_colors), show_annotation_name = FALSE)
      row_ha = rowAnnotation(Cluster = cluster, show_legend = FALSE, 
                             col = list(Cluster = condition_colors), show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"), labels = NULL)
      
      
      htmp <- Heatmap(mat_scaled, cluster_columns = FALSE, cluster_rows = FALSE, name = "Z-score", column_split =condition, row_split =cluster, row_title = NULL, show_row_names = FALSE,
                      show_column_names = FALSE, top_annotation = ha, right_annotation = row_ha, col = heatmap_color, heatmap_legend_param = list(direction = "horizontal"))
      
      
      htmp_save <- draw(htmp, heatmap_legend_side = "bottom")
      
      dev.off()
      jpeg(glue("{dirplots}/{clusterID}.time.DEGs.heatmap.jpeg"), width=600, height=700, units = "px")
      print(htmp_save)
      dev.off()
    }
  }
  
  # ========== 6. EXTRA: ARGs heatmap ==========
  for (i in 1:length(groups)) {
    heatmap_genes <- DEG_complete_results %>%
      filter(contrast %in% time_contrasts,
             significant != "No significant",  #Changed this bit to include also < 5%
             cluster_id == clusterID) %>% select("gene") %>% pull()
    
    if (length(heatmap_genes != 0)) {
      seu.temporal <- DNs.RNA.seu[, DNs.RNA.seu$region == clusterID]
      # seu.temporal@assays[["ATAC"]] <- NULL
      # seu.temporal@assays[["SCT"]] <- NULL
      # seu.temporal@assays[["peaks"]] <- NULL
      
      #Calculate average expression and format columns
      Idents(seu.temporal) <- "orig.ident"
      seu.temporal.avgexp = as.data.frame(AverageExpression(seu.temporal, group.by = 'orig.ident')) 
      names(seu.temporal.avgexp) <- substring(names(seu.temporal.avgexp), 5)
      
      ARG_list <- read_tsv("/fast/AG_Pombo/luna/8_defaultARCreference/4_refineDN_set/2_DEA_muscat_improved/ARG_list_tyssowski.tsv")
      colnames(ARG_list) <- c("Gene", "ARG_class")
      
      ARG_list <- ARG_list %>%filter(ARG_class %in% c("rPRG", "dPRG", "SRG"))
      
      seu.temporal.avgexp.ARG <- seu.temporal.avgexp %>% 
        rownames_to_column("Gene") %>% 
        filter(Gene %in% ARG_list$Gene,)
      
      seu.temporal.avgexp.ARG <- left_join(seu.temporal.avgexp.ARG, ARG_list, by = "Gene")
      ARGclases <- c("rPRG", "dPRG", "SRG")
      ARG_colors <- c("#BA7A12", "#B94227", "#087A87")
      names(ARG_colors) <- ARGclases
      
      seu.temporal.avgexp.ARG.melt <- melt(seu.temporal.avgexp.ARG)
      
      seu.temporal.avgexp.ARG.melt$variable <- factor(seu.temporal.avgexp.ARG.melt$variable, levels = samples)
      seu.temporal.avgexp.ARG.melt$ARG_class <- factor(seu.temporal.avgexp.ARG.melt$ARG_class, levels = ARGclases)
      
      # ARG.RNAstats <- seu.temporal.avgexp.ARG.melt %>% group_by(ARG_class) %>%
      #   dplyr::mutate(MedianCounts = median(nCount_RNA)) %>%
      #   mutate(MedianGenes = median(nFeature_RNA)) %>%
      #   select(seurat_clusters, MedianCounts, MedianGenes) %>% distinct()
      
      ARGs_Gex_plot <- ggplot(seu.temporal.avgexp.ARG.melt, aes(x=ARG_class, y=value, fill = ARG_class, colour = ARG_class))+
        geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
        geom_point(position = position_jitter(width = .15), size = .25)+
        #geom_boxplot(aes(x = as.numeric(ARG_class)+0.25, y = value),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
        scale_y_log10() +
        labs(title = glue("230926: {clusterID} Gene expression per ARG class"), x= "", y= "Gene expression") +
        geom_hline(yintercept = 0.01, color = "red", linetype = "dashed") +
        scale_fill_manual(values = ARG_colors) + 
        scale_color_manual(values = ARG_colors) +
        theme_minimal() +
        theme(legend.position = "bottom") +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        facet_wrap(~variable)
      ARGs_Gex_plot
      ggsave(filename = glue("{dir}{clusterID}.ARGs_Gex_plot_v1.png"), ARGs_Gex_plot, units = "px", device = "png",width=4000,height=4000,dpi = 300)
      
      ARGs_Gex_sample_plot <- ggplot(seu.temporal.avgexp.ARG.melt, aes(x=variable, y=value, fill = variable, colour = variable))+
        geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
        geom_point(position = position_jitter(width = .15), size = .25)+
        #geom_boxplot(aes(x = as.numeric(variable)+0.25, y = value),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
        scale_y_log10() +
        labs(title = glue("230926: {clusterID} Gene expression of ARG genes per sample"), x= "", y= "Gene expression") +
        geom_hline(yintercept = 0.01, color = "red", linetype = "dashed") +
        scale_fill_manual(values = sample_colors) + 
        scale_color_manual(values = sample_colors) +
        #theme_minimal() +
        #theme(legend.position = "bottom") +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        facet_wrap(~ARG_class, ncol = 1)
      ARGs_Gex_sample_plot
      ggsave(filename = glue("{dir}{clusterID}.ARGs_Gex_plot_v2.png"), ARGs_Gex_sample_plot, units = "px", device = "png",width=4000,height=4000,dpi = 300)
      
      rm(sce)
    }
  }
}

# ========== 7. Backup: complete analysis ==========
#save.image("231003.VTAvsSN.muscat.rds")

