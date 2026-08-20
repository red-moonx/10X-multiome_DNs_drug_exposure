# ===========================================
# Script Title: 1.3 Additional markers -dotplot
# Author: Luna Zea Redondo
# Date: 2026-05-02
# Description:
#   This script explores a larger set of marker genes for different neuronal populations in the VTA
# ===========================================


library(patchwork)
library(ggplot2)

da_markers <- c(
  # Classic
  "Th", "Ddc", "Slc6a3", "Slc18a2",
  # Reliable in snRNA-seq
  "Rit2", "Sv2c", "Cacna2d2",
  # TFs
  "Nr4a2", "Pitx3", "Foxa2", "Lmx1b", "En1", "Otx2",
  # Subtype
  "Aldh1a1", "Etv1", "Calb1", "Calb2", "Grp", "Ntsr1"
)

gaba_markers <- c(
  # Core
  "Gad1", "Gad2", "Slc32a1",
  # Supportive
  "Slc6a1", "Nxph1", "Crhbp",
  # TFs
  "Dlx1", "Dlx2", "Nr2f2", "Meis2",
  # Subtypes
  "Pvalb", "Sst", "Vip", "Reln", "Calb1"
)

gluta_markers <- c(
  # Core
  "Slc17a6", "Slc17a7", "Slc17a8",
  # Supportive
  "Gria1", "Grin1", "Camk2a", "Neurod6",
  # TFs
  "Tbr1", "Nfia", "Nfib",
  # Subtypes
  "Satb2", "Foxp2"
)

p1 <- DotPlot(seu.VTA_DNs,
              features = da_markers,
              group.by = "seurat_clusters",
              scale = FALSE) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  ggtitle("Dopaminergic Markers") +
  theme(axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- DotPlot(seu.VTA_DNs,
              features = gaba_markers,
              group.by = "seurat_clusters",
              scale = FALSE) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  ggtitle("GABAergic Markers") +
  theme(axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

p3 <- DotPlot(seu.VTA_DNs,
              features = gluta_markers,
              group.by = "seurat_clusters",
              scale = FALSE) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "darkgreen") +
  ggtitle("Glutamatergic Markers") +
  theme(axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

p1 / p2 / p3 +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    title = "Neurotransmitter Identity Across VTA Clusters",
    theme = theme(plot.title = element_text(hjust = 0.5, 
                                            face = "bold", 
                                            size = 14))
  )

# Save
ggsave("VTA_neurotransmitter_dotplots.pdf", 
       width = 14, height = 16, 
       device = "pdf")