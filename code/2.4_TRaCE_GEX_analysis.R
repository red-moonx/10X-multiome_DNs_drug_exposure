# ===========================================
# Script Title: 2.4. TRaCE analysis
# Author: Luna Zea Redondo
# Date: 2026-05-15
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

source("/fast/AG_Pombo/luna/2026_rebuttal/config.R", local = FALSE)

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
dir <- "/fast/AG_Pombo/luna/2026_rebuttal/6_TRACE"
setwd(dir)

# ========== 1. Load and Subset Seurat Object ==========
# Load full RNA object and subset to VTA cells from selected late timepoints
seu.VTA_DNs <- readRDS("/fast/AG_Pombo/luna/2026_rebuttal/2_DN-evaluation/250507_seu.VTA_DNs.rds")

seu.VTA_DNs@meta.data <- seu.VTA_DNs@meta.data %>%
  tidyr::separate(
    orig.ident,
    into = c("timepoint", "treatment", "replicate"),
    sep = "_",
    remove = FALSE
  )

#a) Select time points and treatments of interest: we are excluding h1_saline
Idents(seu.VTA_DNs) <- "orig.ident"
samples_totest <- c("h4_saline_R1", "h8_saline_R1", "h24_saline_R1", "d14_saline_R1", #salines
                    "h1_cocaine_R1", "h1_cocaine_R2", "h4_cocaine_R1", "h4_cocaine_R2", "h8_cocaine_R1", "h8_cocaine_R2", #ETP
                    "h24_cocaine_R1", "h24_cocaine_R2", "h24_cocaine_R3", "d14_cocaine_R1", "d14_cocaine_R2", "d14_cocaine_R3") #LTP

sample_totest_colors <- sample_colors[names(sample_colors) %in% samples_totest]
sample_totest_colors <- sample_totest_colors[samples_totest]
seu.VTA_DNs.2reps <- subset(x = seu.VTA_DNs, idents = samples_totest)


#c) NEW: extract only LATE DEGs (24/05/29):
seu.VTA_DNs.late <- subset(seu.VTA_DNs.2reps, subset = simpleIdent %in% c("saline", "h24_cocaine", "d14_cocaine") & orig.ident != "h1_saline_R1")


# ========== 2. Metadata cleanup and export  ==========
# Fix replicate metadata
# Remove scale.data before saving H5Seurat/h5ad formats
# 24/02/10
selected_metadata <- seu.VTA_DNs.late@meta.data[, c("timepoint", "treatment", "replicate")]
selected_metadata_mod <- selected_metadata %>% 
  dplyr::mutate(replicate = ifelse(timepoint == "h4" & treatment == "saline", "R1", 
                                   ifelse(timepoint == "h8" & treatment == "saline", "R2",
                                          ifelse(timepoint == "h24" & treatment == "saline", "R3",
                                                 ifelse(timepoint == "d14" & treatment == "saline", "R4", replicate)))))


# Update the Seurat object with the new metadata
seu.VTA_DNs.late@meta.data <- selected_metadata_mod
seu.VTA_DNs.late[["RNA"]]$scale.data <- NULL 

obj <- seu.VTA_DNs.late

# convert Assay5 -> Assay
obj[["RNA"]] <- as(obj[["RNA"]], Class = "Assay")

# keep only counts + data
obj <- DietSeurat(
  obj,
  assays = "RNA",
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE
)


SaveH5Seurat(obj, filename = "250515.DNs.VTA.RNA.seu.sub.data.879.h5Seurat")
Convert("250515.DNs.VTA.RNA.seu.sub.data.879.h5Seurat", dest = "h5ad")


# ========== 3. Extract Late DEGs ==========
# Load Muscat DEG table, filter for late contrasts (h24/d14), save gene list
DEG_complete_results <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/3_muscat/deg_results/260513_DEG_complete_results_Analysis1.tsv")
vta.late.degs <- DEG_complete_results %>%
  dplyr::filter(significant != "No significant", 
                contrast %in% c("h24_cocaine-saline", "d14_cocaine-saline", "d14_cocaine-h24_cocaine")) %>%
  dplyr::select(gene) %>% distinct()
#write_tsv(vta.late.degs, "250615_vta.late.DEGs.txt")




# ========== $%$%$%$% ==========

# The work from this script will be input for the jupyter notebook (2.5)
# Visual inspection of clusters in 2.5, to classify in different dynamic classes (this step could be improved and automatized in future versions, if needed)

# ========== $%$%$%$% ==========


# ========== 4. Define K-means Memory Clusters ==========
# Define groupings (memory, recovered, delayed) and assign direction (up/down)



# kgex_recovered <- c(2, 6, 8, 13, 18, 25, 29, 32, 33, 38, 40)
# kgex_memory <- c(1, 5, 9, 11, 14, 19, 20, 21, 26, 27, 28, 31, 37)
# kgex_delayed <- c(3, 4, 7, 10, 12, 15, 16, 17, 22, 23, 24, 30, 34, 35, 36, 39)
# 
# kgex_down <- c(1, 4, 6, 9, 10, 12, 13, 16, 18, 21, 23, 24, 26, 27, 30, 32, 35, 37, 40)
# kgex_up <- c(2, 3, 5, 7, 8, 11, 14, 15, 17, 19, 20, 22, 25, 28, 29, 31, 33, 34, 36, 38, 39)


kgex_recovered <- c(2, 5, 8, 9, 18, 20, 28, 29, 34, 39)
kgex_memory <- c(1, 7, 11, 12, 14, 16, 19, 26, 27, 30, 32, 33, 35, 40)
kgex_delayed <- c(3, 4, 6, 10, 13, 15, 17, 21, 22, 23, 24, 25, 31, 36, 37, 38)

kgex_up <- c(2, 3, 7, 8, 10, 11, 12, 13, 15, 18, 20, 21, 24, 25, 29, 30, 32, 33, 40)
kgex_down <- c(1, 4, 5, 6, 9, 14, 16, 17, 19, 22, 23, 26, 27, 28, 31, 34, 35, 36, 37, 38, 39)


#Sanity checks
combined_vector <- c(kgex_recovered, kgex_memory, kgex_delayed)
all(1:40 %in% combined_vector)
length(unique(combined_vector)) == 40

# Add info to master table 
memory_DEGs_kmeans <- read_csv("260515_dynamic_cluster_df_version2.csv") %>% 
  dplyr::rename(memory_GEX_kmeans = "level_0", gene = "0") %>% 
  dplyr::mutate(memory_GEX_kmeans = str_replace(memory_GEX_kmeans, "cluster", ""), 
                memory_class = case_when(
                  memory_GEX_kmeans %in% kgex_recovered ~ "recovered",
                  memory_GEX_kmeans %in% kgex_memory ~ "memory",
                  memory_GEX_kmeans %in% kgex_delayed ~ "delayed"
                ),
                direction = case_when(
                  memory_GEX_kmeans %in% kgex_up ~ "up",
                  memory_GEX_kmeans %in% kgex_down ~ "down"))

memory_DEGs_kmeans <- memory_DEGs_kmeans %>%
  dplyr::mutate(gene = gsub(":", "-", gene)) %>%  # fix naming inconsistency
  dplyr::filter(!grepl("^mt-", gene))              # remove mitochondrial genes

# ========== 5. Merge DEG Clusters with Transient Genes ==========
# Extract remaining DEGs not in memory clusters ('transient')
# Combine and export full DEG classification table

all_DEGs <- DEG_complete_results %>% 
  dplyr::filter(significant !=  "No significant") %>% 
  dplyr::select(gene, contrast, logFC, significant) %>%
  dplyr::filter(!grepl("^mt-", gene))

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


# Save DEGs table with kmeans info alone.
DEG_complete_results_kmeans <- rbind(transient_DEGs_first_occurrences, memory_DEGs_kmeans, by = "gene")
write_tsv(DEG_complete_results_kmeans, "260516_DEG_complete_results_kmeans.tsv")


# gene_track <- c("Cartpt", "Penk", "Npy", "Bdnf", "Drd2", "Th", "Foxa1")
# DEG_complete_results_kmeans[DEG_complete_results_kmeans$gene %in% gene_track, ]

test <- read_tsv("/fast/AG_Pombo/luna/2026_rebuttal/6_TRACE/260516_DEG_complete_results_kmeans.tsv")

