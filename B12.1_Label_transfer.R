library(qs2)
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# Input and Output
# Input is integrated object
output_base  <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Samples_path <- file.path(output_base, "Samples", "Samples.qs2")
Plots_dir <- file.path(output_base, "Plots")

M0_RNA <- qs_read(file.path(SeuratObjects_dir, "M0_RNA.qs2"))
FANS <- qs_read(file.path(SeuratObjects_dir, "FANS_Integrated.qs2"))
Samples <- qs_read(Samples_path)

FANS <- NormalizeData(FANS, assay = "RNA")
FANS <- JoinLayers(FANS, assay = "RNA")


# -------------------------------------------------------------------------
# Find anchors and transfer labels
FANS.anchors <- FindTransferAnchors(
  reference = M0_RNA,
  reference.assay = "RNA",
  query.assay = "RNA",
  query = FANS,
  dims = 1:30,
  reference.reduction = NULL
)

length(FANS.anchors@anchor.features)
head(FANS.anchors@anchor.features, 100)

FANS_WNN_L1_Predictions <- TransferData(anchorset = FANS.anchors, refdata = M0_RNA$WNN_L1, query.assay = "RNA", dims = 1:30)
FANS_WNN_L15_Predictions <- TransferData(anchorset = FANS.anchors, refdata = M0_RNA$WNN_L15, query.assay = "RNA", dims = 1:30)
FANS_WNN_L2_Predictions <- TransferData(anchorset = FANS.anchors, refdata = M0_RNA$WNN_L2, query.assay = "RNA", dims = 1:30)
FANS_WNN_L25_Predictions <- TransferData(anchorset = FANS.anchors, refdata = M0_RNA$WNN_L25, query.assay = "RNA", dims = 1:30)
FANS_WNN_L3_Predictions <- TransferData(anchorset = FANS.anchors, refdata = M0_RNA$WNN_L3, query.assay = "RNA", dims = 1:30)
FANS_RNA_L4_Predictions <- TransferData(anchorset = FANS.anchors, refdata = M0_RNA$RNA_L4, query.assay = "RNA", dims = 1:30)

FANS <- AddMetaData(FANS, metadata = FANS_WNN_L1_Predictions$predicted.id, col.name = "ID_WNN_L1_Predicted")
FANS <- AddMetaData(FANS, metadata = FANS_WNN_L15_Predictions$predicted.id, col.name = "ID_WNN_L15_Predicted")
FANS <- AddMetaData(FANS, metadata = FANS_WNN_L2_Predictions$predicted.id, col.name = "ID_WNN_L2_Predicted")
FANS <- AddMetaData(FANS, metadata = FANS_WNN_L25_Predictions$predicted.id, col.name = "ID_WNN_L25_Predicted")
FANS <- AddMetaData(FANS, metadata = FANS_WNN_L3_Predictions$predicted.id, col.name = "ID_WNN_L3_Predicted")
FANS <- AddMetaData(FANS, metadata = FANS_RNA_L4_Predictions$predicted.id, col.name = "ID_RNA_L4_Predicted")


# -------------------------------------------------------------------------
# Rename the labels in the metadata
FANS$tdp_status <- dplyr::recode(
  FANS$tdp_status, 
  "TDP_negative" = "TDP-43_low", 
  "TDP_positive" = "TDP-43_high",
  "TDP-43-"      = "TDP-43_low",
  "TDP-43+"      = "TDP-43_high"
)

# Ensure it's a factor so the legend order is consistent
FANS$tdp_status <- factor(FANS$tdp_status, levels = c("TDP-43_low", "TDP-43_high"))

# Target colors
tdp_colors <- c("TDP-43_low" = "#800000", "TDP-43_high" = "#DAA520")

# Expression vector for mathematical superscript in ggplot legend
legend_labels_tdp <- expression(
  "TDP-43"^"Low", 
  "TDP-43"^"High"
)
names(legend_labels_tdp) <- c("TDP-43_low", "TDP-43_high")

# Plot and save each prediction level
prediction_cols <- c(
  "ID_WNN_L1_Predicted",
  "ID_WNN_L15_Predicted",
  "ID_WNN_L2_Predicted",
  "ID_WNN_L25_Predicted",
  "ID_WNN_L3_Predicted",
  "ID_RNA_L4_Predicted"
)

titles <- c(
  "ID_WNN_L1_Predicted"  = "Level 1 resolution",
  "ID_WNN_L15_Predicted" = "Level 1.5 resolution",
  "ID_WNN_L2_Predicted"  = "Level 2 resolution",
  "ID_WNN_L25_Predicted" = "Level 2.5 resolution",
  "ID_WNN_L3_Predicted"  = "Level 3 resolution",
  "ID_RNA_L4_Predicted"  = "Level 4 resolution"
)

for (col in prediction_cols) {
  # Create the left UMAP 
  p1 <- DimPlot(FANS, reduction = "Umap", group.by = "tdp_status") + 
    scale_color_manual(values = tdp_colors, labels = legend_labels_tdp) +
    ggtitle("TDP-43 status")
  
  # Create the right UMAP
  p2 <- DimPlot(FANS, reduction = "Umap", group.by = col) +
    ggtitle(titles[col])
  
  # Combine the plots 
  p_combined <- p1 | p2
  
  ggsave(
    filename = file.path(Plots_dir, paste0("12_LabelTransfer_", col, ".svg")),
    plot     = p_combined,
    device   = "svg",
    width    = 18,
    height   = 6,
    dpi      = 300
  )
}

# Save
qs_save(FANS, file.path(SeuratObjects_dir, "FANS_Integrated_LabelsTransferred.qs2"))
