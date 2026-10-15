library(data.table)
library(qs2)
library(tidyverse)
library(Seurat)

# Input and Output
# Input is doublet filtered object
output_base  <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Donors_path <- file.path(output_base, "Samples", "Singularity_Donors.qs2")

Donors <- qs_read(Donors_path)
FANS_Seurats_Donors <- qs_read(file.path(SeuratObjects_dir, "FANS_Donor_SeuratObjects_Filtered_Singlets.qs2"))

# -------------------------------------------------------------------------
# Merge Seurat objects 
FANS_Donors <- merge(
  FANS_Seurats_Donors[[1]],
  y = FANS_Seurats_Donors[c(2:length(FANS_Seurats_Donors))],
  add.cell.ids = names(FANS_Seurats_Donors)
)

FANS_Donors$ID <- FANS_Donors$orig.ident
FANS_Donors$tdp_status  <- Donors$tdp_status[match(FANS_Donors$ID, Donors$ID)]
FANS_Donors$neun_status <- Donors$neun_status[match(FANS_Donors$ID, Donors$ID)]
FANS_Donors$Library     <- Donors$Library[match(FANS_Donors$ID, Donors$ID)]
FANS_Donors$scDblFinder_Skipped <- Donors$scDblFinder_Skipped[match(FANS_Donors$ID, Donors$ID)]

table(FANS_Donors$tdp_status)
table(FANS_Donors$neun_status)

# Save
qs_save(FANS_Donors, file.path(SeuratObjects_dir, "FANS_Donor_Merged.qs2"))
