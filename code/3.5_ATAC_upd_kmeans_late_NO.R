## ATAC - Memory
## 241114 - kmeans- updated work
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() 
#.rs.restartR()
`%notin%` <- Negate(`%in%`)
.libPaths(c("~/profiles/r_multiome230913/site-library"))
httr::set_config(httr::config(ssl_verifypeer = 0L))

library(muscat)
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(ggplot2)
library(glue)
library(tidyr)
library(harmony)
library(PupillometryR)
library(ggrepel)
library(ggExtra)
library(tidyverse)
library(Cairo)
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(parallel)
library(glue)
library(ggpubr)
library(hexbin)
library(scales)
library(presto)
library(enrichR)
library(DOSE)
library(chromVARmotifs)
library(RColorBrewer)
library(ComplexHeatmap)
library(Pando)
library(patchwork)
library(scater)
library(limma)
library(UpSetR)
library(genomation)
library(dendsort)
library(memes)
library(SeuratDisk)

source("/fast/AG_Pombo/luna/2023_vta_multiome/2024/1_GEX/7_240510_GEX_memory_only/240510_functions_upset.sh")
setwd("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/3_241113_kmeans")

####################################################
#1. Modified from 240314_kmeans_DARs_preparation.R
####################################################

#1.1. Get only late cells (saline, 24h and 14d)
################################################
wnn.seu.VTA <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/6_231107_updated_analyses/2_ATAC_alternative_processing/3_231205_ArchR_DARsandMotifs/231219_rpkm_based_DARs/240119_thresholds_setC/240122.wnn.seu.VTA.withMotifs.rds")
wnn.seu.VTA.1291 <- subset(wnn.seu.VTA, subset = orig.ident %notin% c("m30_cocaine_R1", "h1_saline_R1"))
wnn.seu.VTA.late <- subset(wnn.seu.VTA.1291, subset = simpleIdent %in% c("saline", "h24_cocaine", "d14_cocaine") & orig.ident != "h1_saline_R1")

#1.2 Generate adata file for version 2
######################################
wnn.seu.peaks <- wnn.seu.VTA.late
wnn.seu.peaks@assays[["RNA"]] <- NULL
wnn.seu.peaks@assays[["SCT"]] <- NULL
wnn.seu.peaks@assays[["ATAC"]] <- NULL
wnn.seu.peaks@assays$peaks@motifs <- NULL# Addied on 241115: because https://github.com/mojaveazure/seurat-disk/issues/15
DefaultAssay(wnn.seu.peaks) <- "peaks"

#Adding this little bit from GEX kmeans (on 24/06/01: allow for replicate identification in saline samples)
selected_metadata <- wnn.seu.peaks@meta.data[, c("timepoint", "treatment", "replicate")]
selected_metadata_mod <- selected_metadata %>% 
  dplyr::mutate(replicate = ifelse(timepoint == "h4" & treatment == "saline", "R1", 
                                   ifelse(timepoint == "h8" & treatment == "saline", "R2",
                                          ifelse(timepoint == "h24" & treatment == "saline", "R3",
                                                 ifelse(timepoint == "d14" & treatment == "saline", "R4", replicate)))))


# Update the Seurat object with the new metadata
wnn.seu.peaks@meta.data <- selected_metadata_mod
wnn.seu.peaks.data <- DietSeurat(wnn.seu.peaks, counts = TRUE, data = TRUE, scale.data = FALSE)
SaveH5Seurat(wnn.seu.peaks.data, filename = "241113.wnn.seu.peaks.data.805.h5Seurat")
Convert("241113.wnn.seu.peaks.data.805.h5Seurat", dest = "h5ad")

#1.3. Save DARs: Late-DARs: I used the ones obtained in 240515 - are the same-
###########################
#load("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/1_ATAC_masterTable_continued/240314.ongoing.ATACwork.rds")
# VTA_RPKM_formatted_metadata_withDARs <- VTA_RPKM_formatted_metadata_withDARs %>% 
#   dplyr::mutate(diff = ifelse((p_val < pval & logFC >=logfc), "Upregulated (pval < 0.05 and logFC >= 0.25)", 
#                               ifelse((p_val < pval & logFC <=logfc), "Downregulated (pval < 0.05 and logFC <= -0.25)", "No significant")))
# 
# vta.late.DARs <- VTA_RPKM_formatted_metadata_withDARs %>% 
#   dplyr::filter(contrast %in% c("h24_cocaine-saline", "d14_cocaine-h24_cocaine", "d14_cocaine-saline"),
#                 diff != "No significant") %>%
#   dplyr::select(peakID) %>% distinct()
# 
# write_tsv(vta.late.DARs, "240515_vta.late.DARs.txt")

#Now they are not the same, new DARs (from 241113: )
DARs_corrected <- readRDS("/fast/AG_Pombo/luna/2023_vta_multiome/2024/2_ATAC/6_241113_ATACmemory_woh1sal/241113.DARs_corrected.rds")
DARs_corrected_late_toKmeans <- DARs_corrected %>% 
  dplyr::filter(diff == "Yes", 
                contrast %in% c("h24_cocaine-saline","d14_cocaine-saline", "d14_cocaine-h24_cocaine")) %>% 
  dplyr::select(peakID) %>% distinct()

#write_tsv(DARs_corrected_late_toKmeans, "241115_vta.late.DARs.txt")  


#1.4. Get all samples
################################################
wnn.seu.VTA.all <- wnn.seu.VTA.1291

#1.5 Generate adata file for version 2: to plot all timepoints
################################################################
wnn.seu.peaks <- wnn.seu.VTA.all
wnn.seu.peaks@assays[["RNA"]] <- NULL
wnn.seu.peaks@assays[["SCT"]] <- NULL
wnn.seu.peaks@assays[["ATAC"]] <- NULL
wnn.seu.peaks@assays$peaks@motifs <- NULL# Addied on 241115: because https://github.com/mojaveazure/seurat-disk/issues/15
DefaultAssay(wnn.seu.peaks) <- "peaks"

#Adding this little bit from GEX kmeans (on 24/06/01: allow for replicate identification in saline samples)
selected_metadata <- wnn.seu.peaks@meta.data[, c("timepoint", "treatment", "replicate")]
selected_metadata_mod <- selected_metadata %>% 
  dplyr::mutate(replicate = ifelse(timepoint == "h4" & treatment == "saline", "R1", 
                                   ifelse(timepoint == "h8" & treatment == "saline", "R2",
                                          ifelse(timepoint == "h24" & treatment == "saline", "R3",
                                                 ifelse(timepoint == "d14" & treatment == "saline", "R4", replicate)))))


# Update the Seurat object with the new metadata
wnn.seu.peaks@meta.data <- selected_metadata_mod
wnn.seu.peaks.data <- DietSeurat(wnn.seu.peaks, counts = TRUE, data = TRUE, scale.data = FALSE)
# SaveH5Seurat(wnn.seu.peaks.data, filename = "241115.wnn.seu.peaks.data.1291.h5Seurat")
# Convert("241115.wnn.seu.peaks.data.1291.h5Seurat", dest = "h5ad")
# load("241115.kmeans_ATAC_first_part.rds")

####################################################
#2. Create a table with "memory info" for each class.
####################################################
#kmeans selected: late DARs, late Samples, 50k
#After running kmeans and manually classified them. 

kdars_recovered <- c(1,2,10,11,12,19,20,23,25,31,32,33,35,37,38,42,47)
kdars_memory <- c(5,9,13,14,22,24,27,28,30,43,45,50)
kdars_delayed <- c(3,4,6,7,8,15,16,17,18,21,26,29,34,36,39,40,41,44,46,48,49)

kdars_up <-  c(2,4,7,8,9,11,12,14,16,18,19,21,23,24,25,27,29,33,35,37,38,40,41,44,46,49,50)
kdars_down <-  c(1,3,5,6,10,13,15,17,20,22,26,28,30,31,32,34,36,39,42,43,45,47,48)


#Sanity checks
# combined_vector <- c(kdars_up, kdars_down)
# all(1:50 %in% combined_vector)
# length(unique(combined_vector)) == 50

#2.1. Memory classes
######################
memory_DARs_kmeans <- read_csv("241115_DARs_memory_kmeans_version2_smooth.csv") %>% 
  dplyr::rename(memory_DARs_kmeans = "level_0", peakID = "0") %>% 
  dplyr::mutate(memory_DARs_kmeans = str_replace(memory_DARs_kmeans, "cluster", ""), 
                memory_class = case_when(
                  memory_DARs_kmeans %in% kdars_recovered ~ "recovered",
                  memory_DARs_kmeans %in% kdars_memory ~ "memory",
                  memory_DARs_kmeans %in% kdars_delayed ~ "delayed"
                ),
                direction = case_when(
                  memory_DARs_kmeans %in% kdars_up ~ "up",
                  memory_DARs_kmeans %in% kdars_down ~ "down"
                )
  )


#2.2. Transient DARs:
######################
all_DARs <- VTA_RPKM_formatted_metadata_withDARs %>% 
  dplyr::filter(diff !=  "No") %>% 
  dplyr::select(peakID, contrast, logFC, diff)

# Reorder the contrast column as an ordered factor with the desired order
all_DARs$contrast <- factor(all_DARs$contrast, levels = all_contrasts, ordered = TRUE)

# Reorder the rows of the DataFrame based on the ordered contrast column
all_DARs <- all_DARs[order(all_DARs$contrast), ]

transient_DARs <- all_DARs %>%
  dplyr::filter(peakID %notin% memory_DARs_kmeans$peakID) %>% #Connect with previous code (240510_GEX_memory)
  dplyr::select(peakID, contrast, logFC) %>% 
  dplyr::mutate(direction = ifelse(logFC > 0, "up", "down"))

transient_DARs_first_occurrences <- transient_DARs %>%
  group_by(peakID) %>%
  slice_head(n=1) %>%
  ungroup() %>% 
  dplyr::mutate(memory_class = "transient", 
                memory_DARs_kmeans = NA)  %>%
  dplyr::select(memory_DARs_kmeans, peakID, memory_class, direction)


###########################################
#3. Save DARs table with kmeans info alone.
###########################################
DAR_complete_results_kmeans <- rbind(transient_DARs_first_occurrences, memory_DARs_kmeans)
write_tsv(DAR_complete_results_kmeans, "241114_DAR_complete_results_kmeans.tsv")
#save.image("241114.DARs.kmeans.final.rds")
