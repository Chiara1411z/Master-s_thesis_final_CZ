library(qs2)
library(tidyverse)
library(Seurat)

# Input and Output
# Input is integrated object
output_base  <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Donors_path <- file.path(output_base, "Samples", "Singularity_Donors.qs2")
Plots_dir <- file.path(output_base, "Plots")

M0_RNA <- qs_read(file.path(SeuratObjects_dir, "M0_RNA.qs2"))
FANS_Donors <- qs_read(file.path(SeuratObjects_dir, "FANS_Donor_Integrated.qs2"))
Donors <- qs_read(Donors_path)

FANS_Donors <- NormalizeData(FANS_Donors, assay = "RNA")
FANS_Donors <- JoinLayers(FANS_Donors, assay = "RNA")

# -------------------------------------------------------------------------
# Find anchors and transfer labels
FANS_Donors.anchors <- FindTransferAnchors(
  reference = M0_RNA,
  reference.assay = "RNA",
  query.assay = "RNA",
  query = FANS_Donors,
  dims = 1:30,
  reference.reduction = NULL
)

FANS_Donors_WNN_L1_Predictions  <- TransferData(anchorset = FANS_Donors.anchors, refdata = M0_RNA$WNN_L1,  query.assay = "RNA", dims = 1:30)
FANS_Donors_WNN_L15_Predictions <- TransferData(anchorset = FANS_Donors.anchors, refdata = M0_RNA$WNN_L15, query.assay = "RNA", dims = 1:30)
FANS_Donors_WNN_L2_Predictions  <- TransferData(anchorset = FANS_Donors.anchors, refdata = M0_RNA$WNN_L2,  query.assay = "RNA", dims = 1:30)
FANS_Donors_WNN_L25_Predictions <- TransferData(anchorset = FANS_Donors.anchors, refdata = M0_RNA$WNN_L25, query.assay = "RNA", dims = 1:30)
FANS_Donors_WNN_L3_Predictions  <- TransferData(anchorset = FANS_Donors.anchors, refdata = M0_RNA$WNN_L3,  query.assay = "RNA", dims = 1:30)
FANS_Donors_RNA_L4_Predictions  <- TransferData(anchorset = FANS_Donors.anchors, refdata = M0_RNA$RNA_L4,  query.assay = "RNA", dims = 1:30)

FANS_Donors <- AddMetaData(FANS_Donors, metadata = FANS_Donors_WNN_L1_Predictions$predicted.id,  col.name = "ID_WNN_L1_Predicted")
FANS_Donors <- AddMetaData(FANS_Donors, metadata = FANS_Donors_WNN_L15_Predictions$predicted.id, col.name = "ID_WNN_L15_Predicted")
FANS_Donors <- AddMetaData(FANS_Donors, metadata = FANS_Donors_WNN_L2_Predictions$predicted.id,  col.name = "ID_WNN_L2_Predicted")
FANS_Donors <- AddMetaData(FANS_Donors, metadata = FANS_Donors_WNN_L25_Predictions$predicted.id, col.name = "ID_WNN_L25_Predicted")
FANS_Donors <- AddMetaData(FANS_Donors, metadata = FANS_Donors_WNN_L3_Predictions$predicted.id,  col.name = "ID_WNN_L3_Predicted")
FANS_Donors <- AddMetaData(FANS_Donors, metadata = FANS_Donors_RNA_L4_Predictions$predicted.id,  col.name = "ID_RNA_L4_Predicted")

# -------------------------------------------------------------------------
# Plot
FANS_Donors$tdp_status <- dplyr::recode(FANS_Donors$tdp_status, "TDP_negative" = "TDP-43-", 
                                 "TDP_positive" = "TDP-43+")

# Ensure it's a factor so the legend order is consistent
FANS_Donors$tdp_status <- factor(FANS_Donors$tdp_status, levels = c("TDP-43-", "TDP-43+"))

# Target colors
tdp_colors <- c("TDP-43-" = "#800000", "TDP-43+" = "#DAA520")

# Plot and save each prediction level
prediction_cols <- c(
  "ID_WNN_L1_Predicted",
  "ID_WNN_L15_Predicted",
  "ID_WNN_L2_Predicted",
  "ID_WNN_L25_Predicted",
  "ID_WNN_L3_Predicted",
  "ID_RNA_L4_Predicted"
)

for (col in prediction_cols) {
  # Create the left UMAP
  p1 <- DimPlot(FANS_Donors, reduction = "Umap", group.by = "tdp_status") + 
    scale_color_manual(values = tdp_colors) +
    ggtitle("tdp_status")
  
  # Create the right UMAP
  p2 <- DimPlot(FANS_Donors, reduction = "Umap", group.by = col) +
    ggtitle(col)
  
  # Combine the plots
  p_combined <- p1 | p2
  
  ggsave(
    filename = file.path(Plots_dir, paste0("12_Donor_LabelTransfer_", col, ".svg")),
    plot = p_combined,
    device = "svg",
    width = 18,
    height = 6,
    dpi = 300
  )
}

# Save
qs_save(FANS_Donors, file.path(SeuratObjects_dir, "FANS_Donor_Integrated_LabelsTransferred.qs2"))
