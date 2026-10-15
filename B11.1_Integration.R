library(qs2)
library(tidyverse)
library(Seurat)
library(glmGamPoi)

# Input and Output
# Input is filtered object
output_base       <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Plots_dir         <- file.path(output_base, "Plots")

# Load the merged, QC-filtered, annotated library-level object
FANS <- qs_read(file.path(SeuratObjects_dir, "FANS_Merged_QCFiltered.qs2"))

# Sanity check 
unique(FANS$ID)[1:5]
stopifnot(!any(grepl("_donor", unique(FANS$ID))))

# -------------------------------------------------------------------------
# Split Layers, Normalize, and Integrate
FANS <- SCTransform(FANS, vst.flavor = "v2", method = "glmGamPoi")
FANS <- RunPCA(FANS)

FANS <- IntegrateLayers(
  object = FANS,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = TRUE
)

# -------------------------------------------------------------------------
# Cluster cells
FANS[["RNA"]] <- JoinLayers(FANS[["RNA"]])

FANS <- FindNeighbors(FANS, reduction = "harmony", dims = 1:30)
FANS <- FindClusters(FANS, cluster.name = "Seurat_Clusters")
FANS <- RunUMAP(FANS, reduction = "harmony", dims = 1:30, reduction.name = "Umap")

# -------------------------------------------------------------------------
# Plot - combined
p_umap <- DimPlot(
  FANS,
  reduction = "Umap",
  group.by = c("tdp_status", "neun_status", "ID", "Seurat_Clusters"),
  combine = TRUE, label.size = 2
)
ggsave(file.path(Plots_dir, "FANS_Integrated_UMAP.svg"), plot = p_umap, device = "svg", width = 14, height = 10, dpi = 300)

# -------------------------------------------------------------------------
# Plot UMAPs Individually
plot_groups <- c("tdp_status", "neun_status", "ID", "Seurat_Clusters")

for (group in plot_groups) {
  show_label <- ifelse(group == "Seurat_Clusters", TRUE, FALSE)
  
  p_umap <- DimPlot(
    FANS,
    reduction = "Umap",
    group.by = group,
    label = show_label,
    label.size = 4
  ) +
    theme_classic() +
    labs(title = paste("UMAP by", group))
  
  if (group %in% c("ID", "Seurat_Clusters")) {
    p_umap <- p_umap + guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
  }
  
  file_name <- paste0("11_", group, ".svg")
  
  ggsave(
    filename = file.path(Plots_dir, file_name),
    plot = p_umap,
    device = "svg",
    width = 9,
    height = 6,
    dpi = 300
  )
}

# -------------------------------------------------------------------------
# Save data
qs_save(FANS, file.path(SeuratObjects_dir, "FANS_Integrated.qs2"))
