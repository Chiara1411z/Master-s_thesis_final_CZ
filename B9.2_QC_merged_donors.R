library(qs2)
library(tidyverse)
library(Seurat)

# Input and Output
# Input is merged object
output_base  <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Donors_path <- file.path(output_base, "Samples", "Singularity_Donors.qs2")

Donors <- qs_read(Donors_path)
FANS_Donors <- qs_read(file.path(SeuratObjects_dir, "FANS_Donor_Merged.qs2"))

# -------------------------------------------------------------------------
# Define AUX functions
filter_by_MADs <- function(Var, nMads){
  MAD <- mad(Var, center = median(Var, na.rm = TRUE), constant = 1, na.rm = TRUE)
  Var < median(Var, na.rm = TRUE) - nMads*MAD | Var > median(Var, na.rm = TRUE) + nMads*MAD
}

# -------------------------------------------------------------------------
# QC Filter Cells
FANS_Donors$Filter <- !(
  filter_by_MADs(FANS_Donors$nCount_RNA_log1p, nMads = 3) |
    filter_by_MADs(FANS_Donors$Features_by_Counts, nMads = 3) |
    filter_by_MADs(FANS_Donors$Pct_Counts_In_Top_20_Genes, nMads = 5) |
    FANS_Donors$percent.mt > 1 |
    FANS_Donors$percent.ribo > 1
)

table(FANS_Donors$Filter, useNA = "always")

FANS_Donors <- subset(FANS_Donors, subset = Filter)
FANS_Donors

FANS_Donors@meta.data <- FANS_Donors@meta.data %>%
  dplyr::select(-c(MAD_Filter_Flag, scDblFinder.class, Filter))

# Save
qs_save(FANS_Donors, file.path(SeuratObjects_dir, "FANS_Donor_Merged_QCFiltered.qs2"))
