# ===========================================
# Script Title: 4.2_PRC_enzymes_expression.R
# Author: Luna Zea Redondo
# Date: 2024-06-11
# Description:
#   Evaluates the expression of PRC enzymes in mouse multiome/RNA datasets.
#   Main steps:
#     - curate PRC enzyme gene symbols
#     - map human symbols to human Ensembl IDs and mouse orthologs
#     - add RPKM / expression information
#     - merge with DEG information
#     - generate PRC2 / PRC1 / ncPRC1 expression plots
#     - generate MA plots for PRC genes across cocaine contrasts
# ===========================================

#!/usr/bin/env Rscript

# ========== Environment setup ==========
rm(list = ls(all.names = TRUE))
gc()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))
set.seed(1)

# ========== Libraries ==========
library(biomaRt)
library(dplyr)
library(HGNChelper)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ========== Project setup ==========
setwd("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/4_240611_PRCenzymes_expr")
DEGs <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/250205_DEG_complete_results_kmeans.corrected.4243.rds")


# ===========================================
# 1) Curate PRC enzyme gene symbols and map to mouse
# ===========================================

# 1A. Read PRC enzyme list and update gene symbols
loh_PRCenzymes <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/4_240611_PRCenzymes_expr/loh2022_PRCgenes.txt")
all_genes <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002_muscat_VTAvsSN_all_subtypes/Analysis2/deg_results/231002_DEG_complete_results_Analysis2.tsv") %>% 
  dplyr::select(gene) %>% distinct() %>% dplyr::rename(Gene = "gene")

gene_symbols <- loh_PRCenzymes$Gene

updated_genes <- checkGeneSymbols(gene_symbols)
updated_genes <- updated_genes$Suggested.Symbol

loh_PRCenzymes$Updated_gene <- updated_genes 
loh_PRCenzymes <- loh_PRCenzymes %>% 
  separate_rows(Updated_gene, sep = " /// ") %>% 
  dplyr::mutate(Updated_gene = ifelse(is.na(Updated_gene), Gene, Updated_gene))

# 1B. Map human genes to Ensembl IDs and mouse homologs
# Use biomaRt to get the mouse Ensembl IDs
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

human_ensembl_ids <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
                           filters = "hgnc_symbol",
                           values = loh_PRCenzymes$Updated_gene,
                           mart = human_mart) 

loh_PRCenzymes_ensemblID <- left_join(loh_PRCenzymes, human_ensembl_ids, by = c("Updated_gene" = "hgnc_symbol"))

mouse_gene_info <- getBM(
  attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"),
  filters = "ensembl_gene_id",
  values = human_ensembl_ids$ensembl_gene_id,
  mart = human_mart
)

loh_PRCenzymes_ensemblID <- left_join(loh_PRCenzymes_ensemblID, mouse_gene_info, by = "ensembl_gene_id") %>% 
  dplyr::mutate(mouse_gene = ifelse(is.na(ensembl_gene_id) & is.na(mmusculus_homolog_ensembl_gene), Updated_gene, NA))


# Do for all genes:
all_genes <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002_muscat_VTAvsSN_all_subtypes/Analysis2/deg_results/231002_DEG_complete_results_Analysis2.tsv") %>% 
  dplyr::select(gene) %>% distinct() %>% dplyr::rename(Gene = "gene")

# Retrieve the mouse Ensembl IDs for the given gene symbols
mouse_ensembl_ids <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = all_genes$Gene,
  mart = mouse_mart
)

# Merge the retrieved Ensembl IDs with your all_genes table
all_genes_with_ensembl <- all_genes %>%
  left_join(mouse_ensembl_ids, by = c("Gene" = "external_gene_name"))

# 1C. Convert to mouse and manually curate PRC enzyme symbols
PRC_mouse <- loh_PRCenzymes_ensemblID %>%
  dplyr::filter(!is.na(Updated_gene)) %>% 
  rowwise() %>%
  dplyr::mutate(mouse_gene = ifelse(mmusculus_homolog_ensembl_gene %in% all_genes_with_ensembl$ensembl_gene_id,
                                    all_genes_with_ensembl$Gene[match(mmusculus_homolog_ensembl_gene, all_genes_with_ensembl$ensembl_gene_id)],
                                    mouse_gene)) %>%
  dplyr::mutate(mouse_gene = case_when(Updated_gene == "CBX3" ~ "Cbx3", 
                                       Updated_gene == "PHF1" ~ "Phf1", 
                                       Updated_gene == "EPOP" ~ "Epop", 
                                       Updated_gene == "PCGF2" ~ "Pcgf2", 
                                       Updated_gene == "RYBP" ~ "Rybp", 
                                       Updated_gene == "SKP1" ~ "Skp1a",
                                       Updated_gene == "PALI1" ~ "Pali1", 
                                       Updated_gene == "PALI2" ~ "Pali2", 
                                       Updated_gene == "SCMH2" ~ "Scmh2", TRUE ~ mouse_gene)) %>% ungroup() %>% 
  dplyr::mutate(Complex, Substructure, mouse_gene) %>% distinct()

#Save table. 
#write_tsv(PRC_mouse, "240620.PRC_mouse.tsv")


# ===========================================
# 2) Add RPKM / DEG information
# ===========================================
PRC_mouse <- read_tsv("240620.PRC_mouse.tsv") %>% 
  dplyr::filter(complete.cases(mouse_gene))

RNA_and_ATAC <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/1_ATAC_masterTable_continued/240415_TFcleaning/240621_RNA_and_ATAC_toPlot.rds")[, c(1,3,4)] 

RNA_allgenes_RPKM <- RNA_and_ATAC %>% 
  group_by(condition) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = condition, values_from = rna_value) %>%
  select(-row)

PRC_mouse_with_RPKM <- left_join(RNA_allgenes_RPKM, PRC_mouse, by =c("gene_symbol" = "mouse_gene"))

# Add differential info:
memory_DEGs_complete <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/240606_DEG_kmeans_info_coords_GAM_final_table.tsv")
PRC_mouse_with_RPKM_DEGinfo <- left_join(PRC_mouse_with_RPKM, memory_DEGs_complete, by =c("gene_symbol" = "gene")) # This table is completed. 

#write_tsv(PRC_mouse_with_RPKM_DEGinfo, "240703.PRC_mouse_with_RPKM_DEGinfo.tsv")

# ===========================================
# 3) Plot PRC enzyme expression
# ===========================================
# PRC_mouse_with_RPKM_DEGinfo <- PRC_mouse_with_RPKM_DEGinfo %>%
#   dplyr::select(gene_symbol, saline, h1_cocaine, h4_cocaine, h8_cocaine, h24_cocaine, d14_cocaine, 
#                 Complex, Substructure,
#                 memory_class, direction) %>% 
#   distinct()

PRC_complete_df.log2RPKMs.log2FC <- PRC_mouse_with_RPKM_DEGinfo %>% 
  mutate(log2_saline = log2(saline + 1),
         log2_h1_cocaine = log2(h1_cocaine + 1),
         log2_h4_cocaine = log2(h4_cocaine + 1),
         log2_h8_cocaine = log2(h8_cocaine + 1),
         log2_h24_cocaine = log2(h24_cocaine + 1),
         log2_d14_cocaine = log2(d14_cocaine + 1)) %>% 
  mutate(log2FC_h1_cocaine = log2_h1_cocaine - log2_saline,
         log2FC_h4_cocaine = log2_h4_cocaine - log2_saline,
         log2FC_h8_cocaine = log2_h8_cocaine - log2_saline,
         log2FC_h24_cocaine = log2_h24_cocaine - log2_saline,
         log2FC_d14_cocaine = log2_d14_cocaine - log2_saline) %>% 
  mutate(max_log2 = pmax(log2_saline, log2_h1_cocaine, log2_h4_cocaine, log2_h8_cocaine, log2_h24_cocaine, log2_d14_cocaine), 
         max_logFC_abs = pmax(abs(log2FC_h1_cocaine), abs(log2FC_h4_cocaine), abs(log2FC_h8_cocaine), abs(log2FC_h24_cocaine), abs(log2FC_d14_cocaine)), 
         max_logFC = pmax(abs(log2FC_h1_cocaine), abs(log2FC_h4_cocaine), abs(log2FC_h8_cocaine), abs(log2FC_h24_cocaine), abs(log2FC_d14_cocaine)) *
           ifelse(pmax(abs(log2FC_h1_cocaine), abs(log2FC_h4_cocaine), abs(log2FC_h8_cocaine), abs(log2FC_h24_cocaine), abs(log2FC_d14_cocaine)) == abs(log2FC_h1_cocaine), sign(log2FC_h1_cocaine),
                  ifelse(pmax(abs(log2FC_h1_cocaine), abs(log2FC_h4_cocaine), abs(log2FC_h8_cocaine), abs(log2FC_h24_cocaine), abs(log2FC_d14_cocaine)) == abs(log2FC_h4_cocaine), sign(log2FC_h4_cocaine),
                         ifelse(pmax(abs(log2FC_h1_cocaine), abs(log2FC_h4_cocaine), abs(log2FC_h8_cocaine), abs(log2FC_h24_cocaine), abs(log2FC_d14_cocaine)) == abs(log2FC_h8_cocaine), sign(log2FC_h8_cocaine),
                                ifelse(pmax(abs(log2FC_h1_cocaine), abs(log2FC_h4_cocaine), abs(log2FC_h8_cocaine), abs(log2FC_h24_cocaine), abs(log2FC_d14_cocaine)) == abs(log2FC_h24_cocaine), sign(log2FC_h24_cocaine),
                                       ifelse(pmax(abs(log2FC_h1_cocaine), abs(log2FC_h4_cocaine), abs(log2FC_h8_cocaine), abs(log2FC_h24_cocaine), abs(log2FC_d14_cocaine)) == abs(log2FC_d14_cocaine), sign(log2FC_d14_cocaine), NA))))))

# 3A. PRC2 plots
#Plot 1: substructure based
PRC2_toplot_df <- PRC_complete_df.log2RPKMs.log2FC %>% 
  dplyr::mutate(Label = ifelse(Complex == "PRC2", Substructure, "no")) %>% 
  dplyr::mutate(Label = ifelse(is.na(Label), "no", Label)) %>%  
  dplyr::mutate(memory_class = ifelse(is.na(memory_class), "no", memory_class)) %>%  
  dplyr::select(gene_symbol, Label, max_log2, max_logFC, memory_class) %>% distinct()

# Define colors for different labels
label_colors <- c("no" = "gray90",  "PRC2.1" = "#30B06A", "Core" = "#C2AC1E", "PRC2.2" = "#EE6055")

# Separate data frames for "no" genes and other genes
no_genes <- PRC2_toplot_df %>% filter(Label == "no")
other_genes <- PRC2_toplot_df %>% filter(Label %in% c("PRC2.1", "Core", "PRC2.2"))

# Plot using ggplot
prc2_plot1 <- ggplot() +
  geom_point(data = no_genes, aes(x = max_log2, y = max_logFC, color = Label),
             size = 2.5, alpha = 0.8) +
  geom_point(data = other_genes, aes(x = max_log2, y = max_logFC, color = Label),
             size = 2.5, alpha = 0.8) +
  geom_text_repel(data = other_genes, aes(x = max_log2, y = max_logFC, label = gene_symbol),
                  box.padding = 0.5, point.padding = 0.3,
                  size = 3, color = "black",
                  segment.color = "grey", segment.size = 0.2,
                  force = 3, min.segment.length = 0.1,
                  arrow = arrow(length = unit(0.02, "npc"))) +
  scale_color_manual(values = label_colors, na.value = "gray90") +
  labs(x = "Max log2", y = "Max logFC", color = "Label", 
       title = "PRC2: All components coloured by substructure")  +
  theme_minimal()

#Plot 2: memory class based
memory_colors <- c("transient"="#D9481C", "recovered"="#F3A712", "memory"="#008AB8", "delayed" = "#954B77") 

# Separate data frames for "no" genes and other genes
no_genes <- PRC2_toplot_df %>% filter(Label == "no")
other_genes <- PRC2_toplot_df %>% filter(Label %in% c("PRC2.1", "Core", "PRC2.2"), 
                                         memory_class %in% c("transient", "recovered", "memory", "delayed"))

# Plot using ggplot
prc2_plot2 <- ggplot() +
  geom_point(data = no_genes, aes(x = max_log2, y = max_logFC, color = Label),
             size = 2.5, alpha = 0.8) +
  geom_point(data = other_genes, aes(x = max_log2, y = max_logFC, color = memory_class),
             size = 2.5, alpha = 0.8) +
  geom_text_repel(data = other_genes, aes(x = max_log2, y = max_logFC, label = gene_symbol),
                  box.padding = 0.5, point.padding = 0.3,
                  size = 3, color = "black",
                  segment.color = "grey", segment.size = 0.2,
                  force = 3, min.segment.length = 0.1,
                  arrow = arrow(length = unit(0.02, "npc"))) +
  scale_color_manual(values = memory_colors, na.value = "gray90") +
  labs(x = "Max log2", y = "Max logFC", color = "memory_class", 
       title = "PRC2: Only DEGs coloured by memory_class") +
  theme_minimal()

prc2_plot1 + prc2_plot2

# 3B. PRC1-canonical plots
#Plot 1: substructure based
cPRC1_toplot_df <- PRC_complete_df.log2RPKMs.log2FC %>% 
  dplyr::mutate(Label = ifelse(Complex == "cPRC1", Substructure, "no")) %>% 
  dplyr::mutate(Label = ifelse(is.na(Label), "no", Label)) %>%  
  dplyr::mutate(memory_class = ifelse(is.na(memory_class), "no", memory_class)) %>%  
  dplyr::select(gene_symbol, Label, max_log2, max_logFC, memory_class) %>% distinct()

# Define colors for different labels
label_colors <- c("no" = "gray90",  "C_specific" = "#30B06A", "Core" = "#C2AC1E")

# Separate data frames for "no" genes and other genes
no_genes <- cPRC1_toplot_df %>% filter(Label == "no")
other_genes <- cPRC1_toplot_df %>% filter(Label %in% c("C_specific", "Core"))

# Plot using ggplot
cPRC1_plot1 <- ggplot() +
  geom_point(data = no_genes, aes(x = max_log2, y = max_logFC, color = Label),
             size = 2.5, alpha = 0.8) +
  geom_point(data = other_genes, aes(x = max_log2, y = max_logFC, color = Label),
             size = 2.5, alpha = 0.8) +
  geom_text_repel(data = other_genes, aes(x = max_log2, y = max_logFC, label = gene_symbol),
                  box.padding = 0.5, point.padding = 0.3,
                  size = 3, color = "black",
                  segment.color = "grey", segment.size = 0.2,
                  force = 3, min.segment.length = 0.1,
                  arrow = arrow(length = unit(0.02, "npc"))) +
  scale_color_manual(values = label_colors, na.value = "gray90") +
  labs(x = "Max log2", y = "Max logFC", color = "Label", 
       title = "cPRC1: All components coloured by substructure")  +
  theme_minimal()

#Plot 2: memory class based
memory_colors <- c("transient"="#D9481C", "recovered"="#F3A712", "memory"="#008AB8", "delayed" = "#954B77") 

# Separate data frames for "no" genes and other genes
no_genes <- cPRC1_toplot_df %>% filter(Label == "no")
other_genes <- cPRC1_toplot_df %>% filter(Label %in% c("C_specific", "Core"), 
                                          memory_class %in% c("transient", "recovered", "memory", "delayed"))

# Plot using ggplot
cPRC1_plot2 <- ggplot() +
  geom_point(data = no_genes, aes(x = max_log2, y = max_logFC, color = Label),
             size = 2.5, alpha = 0.8) +
  geom_point(data = other_genes, aes(x = max_log2, y = max_logFC, color = memory_class),
             size = 2.5, alpha = 0.8) +
  geom_text_repel(data = other_genes, aes(x = max_log2, y = max_logFC, label = gene_symbol),
                  box.padding = 0.5, point.padding = 0.3,
                  size = 3, color = "black",
                  segment.color = "grey", segment.size = 0.2,
                  force = 3, min.segment.length = 0.1,
                  arrow = arrow(length = unit(0.02, "npc"))) +
  scale_color_manual(values = memory_colors, na.value = "gray90") +
  labs(x = "Max log2", y = "Max logFC", color = "memory_class", 
       title = "cPRC1: Only DEGs coloured by memory_class") +
  theme_minimal()

cPRC1_plot1 + cPRC1_plot2


# 3C. PRC1-noncanonical plots
#Plot 1: substructure based
ncPRC1_toplot_df <- PRC_complete_df.log2RPKMs.log2FC %>% 
  dplyr::mutate(Label = ifelse(Complex == "ncPRC1", Substructure, "no")) %>% 
  dplyr::mutate(Label = ifelse(is.na(Label), "no", Label)) %>%  
  dplyr::mutate(memory_class = ifelse(is.na(memory_class), "no", memory_class)) %>%  
  dplyr::select(gene_symbol, Label, max_log2, max_logFC, memory_class) %>% distinct()

# Define colors for different labels
label_colors <- c("no" = "gray90",  "NC_specific" = "#30B06A", "Core" = "#C2AC1E")

# Separate data frames for "no" genes and other genes
no_genes <- ncPRC1_toplot_df %>% filter(Label == "no")
other_genes <- ncPRC1_toplot_df %>% filter(Label %in% c("NC_specific", "Core"))

# Plot using ggplot
ncPRC1_plot1 <- ggplot() +
  geom_point(data = no_genes, aes(x = max_log2, y = max_logFC, color = Label),
             size = 2.5, alpha = 0.8) +
  geom_point(data = other_genes, aes(x = max_log2, y = max_logFC, color = Label),
             size = 2.5, alpha = 0.8) +
  geom_text_repel(data = other_genes, aes(x = max_log2, y = max_logFC, label = gene_symbol),
                  box.padding = 0.5, point.padding = 0.3,
                  size = 3, color = "black",
                  segment.color = "grey", segment.size = 0.2,
                  force = 3, min.segment.length = 0.1,
                  arrow = arrow(length = unit(0.02, "npc"))) +
  scale_color_manual(values = label_colors, na.value = "gray90") +
  labs(x = "Max log2", y = "Max logFC", color = "Label", 
       title = "ncPRC1: All components coloured by substructure")  +
  theme_minimal()

#Plot 2: memory class based
memory_colors <- c("transient"="#D9481C", "recovered"="#F3A712", "memory"="#008AB8", "delayed" = "#954B77") 

# Separate data frames for "no" genes and other genes
no_genes <- ncPRC1_toplot_df %>% filter(Label == "no")
other_genes <- ncPRC1_toplot_df %>% filter(Label %in% c("NC_specific", "Core"), 
                                           memory_class %in% c("transient", "recovered", "memory", "delayed"))

# Plot using ggplot
ncPRC1_plot2 <- ggplot() +
  geom_point(data = no_genes, aes(x = max_log2, y = max_logFC, color = Label),
             size = 2.5, alpha = 0.8) +
  geom_point(data = other_genes, aes(x = max_log2, y = max_logFC, color = memory_class),
             size = 2.5, alpha = 0.8) +
  geom_text_repel(data = other_genes, aes(x = max_log2, y = max_logFC, label = gene_symbol),
                  box.padding = 0.5, point.padding = 0.3,
                  size = 3, color = "black",
                  segment.color = "grey", segment.size = 0.2,
                  force = 3, min.segment.length = 0.1,
                  arrow = arrow(length = unit(0.02, "npc"))) +
  scale_color_manual(values = memory_colors, na.value = "gray90") +
  labs(x = "Max log2", y = "Max logFC", color = "memory_class", 
       title = "ncPRC1: Only DEGs coloured by memory_class") +
  theme_minimal()

ncPRC1_plot1 + ncPRC1_plot2




# ===========================================
# 4) MA plots for PRC genes across cocaine contrasts
# ===========================================

# 250330 - For figures only
#Plot only when they are not DEGs
DEG_complete_results <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002_muscat_VTAvsSN_all_subtypes/Analysis2/deg_results/231002_DEG_complete_results_Analysis2.tsv") %>% 
  dplyr::filter(cluster_id == "VTA")

loh_PRCenzymes <- read_tsv("240620.PRC_mouse.tsv") %>% 
  dplyr::select(Complex, mouse_gene) %>% 
  dplyr::mutate(Complex = ifelse(Complex %in% c("cPRC1", "ncPRC1"), "PRC1", Complex)) %>% 
  dplyr::distinct()

#Establish all possible contrasts:
all_contrasts <- c("h1_cocaine-saline",
                   "h4_cocaine-saline",
                   "h8_cocaine-saline",
                   "h24_cocaine-saline",
                   "d14_cocaine-saline")


# Create a list to store all plots
MAplots <- list()

for (i in seq_along(all_contrasts)) {
  comparison <- all_contrasts[i]
  
  # Merge PRC info
  plot.data <- DEG_complete_results %>%
    filter(contrast == comparison) %>%
    left_join(loh_PRCenzymes, by = c("gene" = "mouse_gene")) %>%
    mutate(
      gene_class = case_when(
        Complex == "PRC2" ~ "PRC2",
        Complex == "PRC1" ~ "PRC1",
        TRUE ~ "Other"
      ),
      to_label = ifelse(gene_class != "Other" & significant != "No significant", "yes", "no")
    )
  
  # Define manual colors for overlay
  gene_colors <- c("PRC2" = "#FF6663", "PRC1" = "#2364AA")
  
  # Plot
  MA_plot <- ggplot(plot.data, aes(x = logCPM, y = logFC)) +
    
    # Plot all genes as gray dots first
    geom_point(data = plot.data, color = "gray90", size = 1.5) +
    
    # Overlay colored dots for PRC1 / PRC2 genes
    geom_point(
      data = subset(plot.data, gene_class != "Other" & significant != "No significant"),
      aes(color = gene_class),
      size = 1.5
    ) +
    geom_hline(yintercept = c(-0.5, 0.5), colour = "goldenrod", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "goldenrod") +
    
    # Label those same genes with colored text and white background
    geom_label_repel(
      data = subset(plot.data, to_label == "yes"),
      aes(label = gene, color = gene_class),
      fill = "white",  # neutral box background
      box.padding = 1, size = 5, show.legend = FALSE
    ) +
    
    # Apply color scale for points and label text
    scale_color_manual(values = gene_colors, name = "Gene Class") +
    
    # Remove the legend
    guides(color = "none") +
    
    # Theme and layout
    theme_classic(base_size = 12) +
    ylim(-6, 9) +
    theme(
      text = element_text(size = 15),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )
  
  # Store the plot
  MAplots[[paste("MAplot", i, sep = "_")]] <- MA_plot
}

# Combine the first five MA plots into one row
MA_row1 <- (MAplots$MAplot_1 + MAplots$MAplot_2 + MAplots$MAplot_3 + MAplots$MAplot_4 + MAplots$MAplot_5) + 
  plot_layout(nrow = 1)

# Display the row
MA_row1


# Define TIFF output using ggsave
# ggsave(
#   filename = "/fast/AG_Pombo/luna/2025_pdf_files/chapter4_PRCexpression_MAplots_DEGs.tif",
#   plot = MA_row1,
#   width = 14.65359,
#   height = 4.62745,
#   dpi = 300,
#   units = "in",
#   device = "tiff"
# )

# tiff("/fast/AG_Pombo/luna/2025_pdf_files/chapter4_PRCexpression_MAplots_DEGs.tif", width = 14.65359, height = 4.62745, res = 300)
# MA_row1
# dev.off()

