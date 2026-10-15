library(qs2)
library(tidyverse)
library(Seurat)

# Input and Output
# Input is merged object
output_base  <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Samples_path <- file.path(output_base, "Samples", "Samples.qs2")

Samples <- qs_read(Samples_path)
FANS <- qs_read(file.path(SeuratObjects_dir, "FANS_Merged.qs2"))

# -------------------------------------------------------------------------
# Define AUX functions
filter_by_MADs <- function(Var, nMads){
  MAD <- mad(Var, center = median(Var, na.rm = TRUE), constant = 1, na.rm = TRUE)
  Var < median(Var, na.rm = TRUE) - nMads*MAD | Var > median(Var, na.rm = TRUE) + nMads*MAD
}

# -------------------------------------------------------------------------
# QC Filter Cells
FANS$Filter <- !(
  filter_by_MADs(FANS$nCount_RNA_log1p, nMads = 3) |
    filter_by_MADs(FANS$Features_by_Counts, nMads = 3) |
    filter_by_MADs(FANS$Pct_Counts_In_Top_20_Genes, nMads = 5) |
    FANS$percent.mt > 1 |
    FANS$percent.ribo > 1
)

table(FANS$Filter, useNA = "always")

FANS <- subset(FANS, subset = Filter)
FANS

FANS@meta.data <- FANS@meta.data %>%
  dplyr::select(-c(MAD_Filter_Flag, scDblFinder.class, Souporcell_Status, Filter))

# -------------------------------------------------------------------------
# Assay Uncorrected filter
FANS@assays$RNA_Uncorrected <- NULL
FANS

# Save
qs_save(FANS, file.path(SeuratObjects_dir, "FANS_Merged_QCFiltered.qs2"))


ggsave(file.path(Plots_dir, "9_8_FeaturesByCounts_vs_nCount_scatter.svg"), p8, device = "svg", width = 9, height = 6, dpi = 300)
