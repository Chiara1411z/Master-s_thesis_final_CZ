library(data.table)
library(qs2)
library(tidyverse)
library(Seurat)

# Input and Output
# Input is doublet filtered object
output_base  <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Samples_path <- file.path(output_base, "Samples", "Samples.qs2")

Samples <- qs_read(Samples_path)
FANS_Seurats <- qs_read(file.path(SeuratObjects_dir, "FANS_Single_SeuratObjects_Filtered_Singlets.qs2"))

# -------------------------------------------------------------------------
# Merge Seurat objects
FANS <- merge(
  FANS_Seurats[[1]],
  y = FANS_Seurats[c(2:length(FANS_Seurats))],
  add.cell.ids = names(FANS_Seurats)
)

FANS$ID  <- FANS$orig.ident
FANS$tdp_status  <- Samples$tdp_status[match(FANS$ID, Samples$ID)]
FANS$neun_status <- Samples$neun_status[match(FANS$ID, Samples$ID)]
FANS$scDblFinder_Skipped <- Samples$scDblFinder_Skipped[match(FANS$ID, Samples$ID)]

table(FANS$tdp_status)
table(FANS$neun_status)

# Save
qs_save(FANS, file.path(SeuratObjects_dir, "FANS_Merged.qs2"))
