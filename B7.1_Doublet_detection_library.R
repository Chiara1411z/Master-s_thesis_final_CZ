library(qs2)
library(Seurat)
library(tidyverse)
library(data.table)

# Input and Output
# Input is filtered and enriched object
output_base  <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Samples_path <- file.path(output_base, "Samples", "Samples.qs2")

Samples <- qs_read(Samples_path)
FANS <- qs_read(file.path(SeuratObjects_dir, "FANS_Single_SeuratObjects_Filtered_scDblFinder_Marked_Genotypes_Assigned.qs2"))

# -------------------------------------------------------------------------
# Flag libraries with no scDblFinder call (too few cells for scDblFinder to run)
missing_scDbl <- names(FANS)[!sapply(FANS, function(x) "scDblFinder.class" %in% colnames(x@meta.data))]
if (length(missing_scDbl) > 0) {
  message("No scDblFinder call for: ", paste(missing_scDbl, collapse = ", "),
          " (too few cells -- scDblFinder filter skipped for these; Souporcell filter still applied)")
}

# -------------------------------------------------------------------------
# Add data to Samples data
nDoublets_scDblFinder <- sapply(FANS, function(x){
  if (!"scDblFinder.class" %in% colnames(x@meta.data)) return(NA)
  sum(x$scDblFinder.class == "doublet")
})
nDoublets_Souporcell  <- sapply(FANS, function(x) sum(x$Souporcell_Status == "doublet", na.rm = TRUE))
PctDoublets_Souporcell <- sapply(FANS, function(x) sum(x$Souporcell_Status == "doublet", na.rm = TRUE)*100/nrow(x@meta.data))

Samples$nDoublets_scDblFinder  <- nDoublets_scDblFinder[match(Samples$ID, names(nDoublets_scDblFinder))]
Samples$nDoublets_Souporcell   <- nDoublets_Souporcell[match(Samples$ID, names(nDoublets_Souporcell))]
Samples$PctDoublets_Souporcell <- PctDoublets_Souporcell[match(Samples$ID, names(PctDoublets_Souporcell))]
Samples$scDblFinder_Skipped    <- Samples$ID %in% missing_scDbl

# -------------------------------------------------------------------------
# Filter doublets
# Apply BOTH filters in sequence: scDblFinder singlets, then Souporcell singlets.
# Libraries lacking scDblFinder.class are passed through unfiltered at this step.

FANS_Seurats_Singlets <- lapply(FANS, function(Seurat){
  if (!"scDblFinder.class" %in% colnames(Seurat@meta.data)) return(Seurat)
  subset(Seurat, subset = scDblFinder.class == "singlet")
})

sapply(FANS, function(x) dim(x)[2]) |> sum()
sapply(FANS_Seurats_Singlets, function(x) dim(x)[2]) |> sum()

FANS_Seurats_Singlets <- lapply(FANS_Seurats_Singlets, function(Seurat){
  Seurat$Souporcell_Status[is.na(Seurat$Souporcell_Status)] <- "NA"
  subset(Seurat, subset = Souporcell_Status %in% c("NA", "singlet"))
})

sapply(FANS, function(x) dim(x)[2]) |> sum()
sapply(FANS_Seurats_Singlets, function(x) dim(x)[2]) |> sum()

lapply(FANS_Seurats_Singlets, dim) %>% data.frame() %>% t() %>% data.frame() %>%
  rownames_to_column(var = "Sample") -> tmp

Samples$nCells_Filtered <- tmp$X2[match(Samples$ID, tmp$Sample)]
rm(tmp)

# Save
qs_save(FANS_Seurats_Singlets, file.path(SeuratObjects_dir, "FANS_Single_SeuratObjects_Filtered_Singlets.qs2"))
qs_save(Samples, Samples_path)
