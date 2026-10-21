library(qs2)
library(tidyverse)
library(Seurat)
library(glmGamPoi)

# Input and Output 
# Input is annotated object
output_base       <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Plots_dir         <- file.path(output_base, "Plots")
Annotations_dir   <- file.path(output_base, "Annotations")
Donors_path       <- file.path(output_base, "Samples", "Singularity_Donors.qs2")


FANS_Donors <- qs_read(file.path(SeuratObjects_dir, "FANS_Donor_Annotated.qs2"))

# -------------------------------------------------------------------------
# Filter Donors (< 50 cells) to Prevent SCTransform Errors
donor_counts    <- table(FANS_Donors$ID)
low_cell_donors <- names(donor_counts[donor_counts < 50])

if (length(low_cell_donors) > 0) {
  message("Excluding donors with < 50 cells: ", paste(low_cell_donors, collapse = ", "))
  valid_donors <- names(donor_counts[donor_counts >= 50])
  FANS_Donors  <- subset(FANS_Donors, subset = ID %in% valid_donors)
}

# -------------------------------------------------------------------------
# Normalize, Run PCA, and Integrate Layers
FANS_Donors <- SCTransform(FANS_Donors)
FANS_Donors <- RunPCA(FANS_Donors)

FANS_Donors <- IntegrateLayers(
  object = FANS_Donors,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = TRUE
)

# -------------------------------------------------------------------------
# Clustering & Dimensionality Reduction
FANS_Donors <- FindNeighbors(FANS_Donors, reduction = "harmony", dims = 1:30)
FANS_Donors <- FindClusters(FANS_Donors, cluster.name = "Seurat_Clusters")
FANS_Donors <- RunUMAP(FANS_Donors, reduction = "harmony", dims = 1:30, reduction.name = "Umap")

# -------------------------------------------------------------------------
# Plot
p_umap <- DimPlot(
  FANS_Donors,
  reduction = "Umap",
  group.by = c("tdp_status", "neun_status", "ID", "Library", "Seurat_Clusters"),
  combine = TRUE, label.size = 2
)
ggsave(file.path(Plots_dir, "FANS_Donor_Integrated_UMAP.svg"), plot = p_umap, device = "svg", width = 14, height = 10, dpi = 300)


# -------------------------------------------------------------------------
# Plot UMAPs Individually
plot_groups <- c("tdp_status", "neun_status", "ID", "Library", "Seurat_Clusters")

for (group in plot_groups) {
  # Add labels directly onto the UMAP for Seurat_Clusters
  show_label <- ifelse(group == "Seurat_Clusters", TRUE, FALSE)
  
  p_umap <- DimPlot(
    FANS_Donors,
    reduction = "Umap",
    group.by = group,
    label = show_label,
    label.size = 4
  ) + 
    theme_classic() +
    labs(title = paste("UMAP by", group))
  
  # Format legends 
  if (group %in% c("ID", "Seurat_Clusters", "Library")) {
    p_umap <- p_umap + guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
  }
  
  # Construct the file name 
  file_name <- paste0("11_donor_", group, ".svg")
  
  # Save the individual plot
  ggsave(
    filename = file.path(Plots_dir, file_name), 
    plot = p_umap, 
    device = "svg", 
    width = 9,
    height = 6, 
    dpi = 300
  )
}

# Save
qs_save(FANS_Donors, file.path(SeuratObjects_dir, "FANS_Donor_Integrated.qs2"))
