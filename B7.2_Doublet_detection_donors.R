library(qs2)
library(Seurat)
library(tidyverse)
library(data.table)

# Input and Output
# Input is filtered and enriched object
output_base  <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Donors_path <- file.path(output_base, "Samples", "Singularity_Donors.qs2")

Donors <- qs_read(Donors_path)
FANS_Seurats_Donors <- qs_read(file.path(SeuratObjects_dir, "FANS_Donor_SeuratObjects_Filtered_scDblFinder_Marked.qs2"))

# -------------------------------------------------------------------------
# Add data to Donors data
missing_scDbl_donors <- names(FANS_Seurats_Donors)[!sapply(FANS_Seurats_Donors, function(x) "scDblFinder.class" %in% colnames(x@meta.data))]
if (length(missing_scDbl_donors) > 0) {
  message("No scDblFinder call for: ", paste(missing_scDbl_donors, collapse = ", "),
          " (too few cells -- scDblFinder filter skipped for these donors)")
}

nDoublets_scDblFinder <- sapply(FANS_Seurats_Donors, function(x){
  if (!"scDblFinder.class" %in% colnames(x@meta.data)) return(NA)
  sum(x$scDblFinder.class == "doublet")
})
Donors$nDoublets_scDblFinder <- nDoublets_scDblFinder[match(Donors$ID, names(nDoublets_scDblFinder))]
Donors$scDblFinder_Skipped   <- Donors$ID %in% missing_scDbl_donors

# -------------------------------------------------------------------------
# Filter doublets
FANS_Seurats_Donors_Singlets <- lapply(FANS_Seurats_Donors, function(Seurat){
  if (!"scDblFinder.class" %in% colnames(Seurat@meta.data)) return(Seurat)
  subset(Seurat, subset = scDblFinder.class == "singlet")
})

sapply(FANS_Seurats_Donors, function(x) dim(x)[2]) |> sum()
sapply(FANS_Seurats_Donors_Singlets, function(x) dim(x)[2]) |> sum()

lapply(FANS_Seurats_Donors_Singlets, dim) %>% data.frame() %>% t() %>% data.frame() %>%
  rownames_to_column(var = "ID") -> tmp

Donors$nCells_Filtered <- tmp$X2[match(Donors$ID, tmp$ID)]
rm(tmp)

# Save
qs_save(FANS_Seurats_Donors_Singlets, file.path(SeuratObjects_dir, "FANS_Donor_SeuratObjects_Filtered_Singlets.qs2"))
qs_save(Donors, Donors_path)
