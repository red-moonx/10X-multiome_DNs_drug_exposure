# ===========================================
# Script Title: 8 GEX Memory Classes
# Author: [Your Name]
# Date: 2025-12-15
# Description:
#   This script classify VTA-DN DEGs in eigh unique dynamic classes: late and early response.
#   It assigns memory-related categories (memory, recovered, delayed) and directionality (up/down) to DEG clusters.
#   It also extracts transient DEGs and saves a final unified classification table.
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
library(dplyr)
library(tidyr)
library(tibble)
library(glue)
library(ggplot2)
library(readr)
library(stringr)
library(SeuratDisk)

# ========== Set Working Directory ==========
dir <- "/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/8_240510_GEX_memory_only_v2/1_240529_GEX_kmeans"
setwd(dir)

# ========== 1. Load and Subset Seurat Object ==========
# Load full RNA object and subset to VTA cells from selected late timepoints
load("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002.RNAonly.excluding_putative_SN_cells.rds")
DNs.RNA.seu #2411 cells; "region" separates VTA from SN cells


#a) Select time points and treatments of interest: we are excluding h1_saline
Idents(DNs.RNA.seu) <- "orig.ident"
samples_totest <- c("h4_saline_R1", "h8_saline_R1", "h24_saline_R1", "d14_saline_R1", #salines
                    "h1_cocaine_R1", "h1_cocaine_R2", "h4_cocaine_R1", "h4_cocaine_R2", "h8_cocaine_R1", "h8_cocaine_R2", #ETP
                    "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3") #LTP

sample_totest_colors <- sample_colors[names(sample_colors) %in% samples_totest]
sample_totest_colors <- sample_totest_colors[samples_totest]
DNs.RNA.seu.2reps <- subset(x = DNs.RNA.seu, idents = samples_totest)


#b) Extract VTA-DN data
DNs.RNA.seu.VTA <- subset(DNs.RNA.seu.2reps, subset = region == "VTA")

#c) NEW: extract only LATE DEGs (24/05/29):
DNs.RNA.seu.VTA.late <- subset(DNs.RNA.seu.VTA, subset = simpleIdent %in% c("saline", "h24_cocaine", "d14_cocaine") & orig.ident != "h1_saline_R1")


# ========== 2. Metadata cleanup and export  ==========
# Fix replicate metadata
# Remove scale.data before saving H5Seurat/h5ad formats
# 24/02/10
selected_metadata <- DNs.RNA.seu.VTA.late@meta.data[, c("timepoint", "treatment", "replicate")]
selected_metadata_mod <- selected_metadata %>% 
  dplyr::mutate(replicate = ifelse(timepoint == "h4" & treatment == "saline", "R1", 
                                   ifelse(timepoint == "h8" & treatment == "saline", "R2",
                                          ifelse(timepoint == "h24" & treatment == "saline", "R3",
                                                 ifelse(timepoint == "d14" & treatment == "saline", "R4", replicate)))))


# Update the Seurat object with the new metadata
DNs.RNA.seu.VTA.late@meta.data <- selected_metadata_mod
DNs.RNA.seu.VTA.late[["RNA"]]$scale.data <- NULL 

#SeuratDisk is saving the scale.data. We have to delete it, so Seuratdisk uses "data" 
DNs.RNA.seu.VTA.late.data <- DietSeurat(DNs.RNA.seu.VTA.late, counts = TRUE, data = TRUE, scale.data = FALSE)
SaveH5Seurat(DNs.RNA.seu.VTA.data, filename = "240529.DNs.VTA.RNA.seu.sub.data.893.h5Seurat")
Convert("240529.DNs.VTA.RNA.seu.sub.data.893.h5Seurat", dest = "h5ad")


# ========== 3. Extract Late DEGs ==========
# Load Muscat DEG table, filter for late contrasts (h24/d14), save gene list
vta.late.degs <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/231024_discrete_updated/231026_muscat/Analysis2/deg_results/231002_DEG_complete_results_Analysis2.tsv") %>% 
  dplyr::filter(cluster_id == "VTA", significant != "No significant", 
                contrast %in% c("h24_cocaine-saline", "d14_cocaine-saline", "d14_cocaine-h24_cocaine")) %>%
  dplyr::select(gene) %>% distinct()
write_tsv(vta.late.degs, "240529_vta.late.DEGs.txt")




# ========== $%$%$%$% ==========

# The work from this script will be input for the jupyter notebook (2.5)
# Visual inspection of clusters in 2.5, to classify in different dynamic classes (this step could be improved and automatized in future versions, if needed)

# ========== $%$%$%$% ==========


# ========== 4. Define K-means Memory Clusters ==========
# Define groupings (memory, recovered, delayed) and assign direction (up/down)

kgex_recovered <- c(1,6,7,11,14,17,20,29,30,31,32,33,36,40,49)
kgex_memory <- c(3,8,13,16,18,22,26,27,28,34,38,41,46,48)
kgex_delayed <- c(2,4,5,9,10,12,15,19,21,23,24,25,35,37,39,42,43,44,45,47,50)

kgex_up <- c(1,2,5,6,8,10,14,15,16,22,23,25,27,28,30,31,32,33,35,38,39,43,45,46,50)
kgex_down <-  c(3,4,7,9,11,12,13,17,18,19,20,21,24,26,29,34,36,37,40,41,42,44,47,48,49)

#Sanity checks
combined_vector <- c(kgex_up, kgex_down)
all(1:50 %in% combined_vector)
no_duplicates <- length(unique(combined_vector)) == 50

# Add info to master table 

memory_DEGs_kmeans <- read_csv("240530_dynamic_cluster_df_version2.csv") %>% 
  dplyr::rename(memory_GEX_kmeans = "level_0", gene = "0") %>% 
  dplyr::mutate(memory_GEX_kmeans = str_replace(memory_GEX_kmeans, "cluster", ""), 
                memory_class = case_when(
                  memory_GEX_kmeans %in% kgex_recovered ~ "recovered",
                  memory_GEX_kmeans %in% kgex_memory ~ "memory",
                  memory_GEX_kmeans %in% kgex_delayed ~ "delayed"
                ),
                direction = case_when(
                  memory_GEX_kmeans %in% kgex_up ~ "up",
                  memory_GEX_kmeans %in% kgex_down ~ "down"
                )
  )

# ========== 5. Merge DEG Clusters with Transient Genes ==========
# Extract remaining DEGs not in memory clusters ('transient')
# Combine and export full DEG classification table

all_DEGs <- read_tsv("/fast/AG_Pombo/luna/2023_vta_multiome/3_conos/1_230418_conos_DN_eval/2_temporalRNA/230907.discrete.analysis.final/230926.excluding_putative_SN_cells/230929.max.resolution.RNAonly/231002_muscat_VTAvsSN_all_subtypes/Analysis2/deg_results/231002_DEG_complete_results_Analysis2.tsv") %>% 
  dplyr::filter(cluster_id == "VTA", significant !=  "No significant") %>% 
  dplyr::select(gene, contrast, logFC, significant)

# Reorder the contrast column as an ordered factor with the desired order
all_DEGs$contrast <- factor(all_DEGs$contrast, levels = all_contrasts, ordered = TRUE)

# Reorder the rows of the DataFrame based on the ordered contrast column
all_DEGs <- all_DEGs[order(all_DEGs$contrast), ]

transient_DEGs <- all_DEGs %>%
  dplyr::filter(gene %notin% memory_DEGs_kmeans$gene) %>% #Connect with previous code (240510_GEX_memory)
  dplyr::select(gene, contrast, logFC) %>% 
  dplyr::mutate(direction = ifelse(logFC > 0, "up", "down"))

transient_DEGs_first_occurrences <- transient_DEGs %>%
  group_by(gene) %>%
  slice_head(n=1) %>%
  ungroup() %>% 
  dplyr::mutate(memory_class = "transient", 
                memory_GEX_kmeans = NA)  %>%
  dplyr::select(memory_GEX_kmeans, gene, memory_class, direction)


# Save DEgs table with kmeans info alone.
DEG_complete_results_kmeans <- rbind(transient_DEGs_first_occurrences, memory_DEGs_kmeans, by = "gene")
write_tsv(DEG_complete_results_kmeans, "240601_DEG_complete_results_kmeans.tsv")











