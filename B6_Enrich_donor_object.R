library(qs2)
library(Seurat)
library(tidyverse)
library(data.table)

# Input and Output
# Input is filtered seurat object
output_base  <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Samples_path <- file.path(output_base, "Samples", "Samples.qs2")

Souporcell_Output_Dir <- "Z.../"

Samples <- qs_read(Samples_path)
FANS <- qs_read(file.path(SeuratObjects_dir, "FANS_Single_SeuratObjects_Filtered_scDblFinder_Marked.qs2"))

Souporcell_Output_Dirs <- list.dirs(Souporcell_Output_Dir, recursive = FALSE, full.names = FALSE)

setdiff(names(FANS), Souporcell_Output_Dirs)

Souporcell_Clusters <- lapply(
  setNames(Souporcell_Output_Dirs, nm = Souporcell_Output_Dirs),
  function(x){
    fread(file.path(Souporcell_Output_Dir, x, "clusters.tsv"), data.table = FALSE)
  }
)

Souporcell_Clusters <- Souporcell_Clusters[
  match(names(FANS), names(Souporcell_Clusters))[!is.na(match(names(FANS), names(Souporcell_Clusters)))]
]

Souporcell_AmbientRNA <- lapply(
  setNames(Souporcell_Output_Dirs, nm = Souporcell_Output_Dirs),
  function(x){
    fread(file.path(Souporcell_Output_Dir, x, "ambient_rna.txt"), data.table = FALSE)
  }
)
Souporcell_AmbientRNA <- sapply(
  Souporcell_AmbientRNA,
  function(x) str_replace_all(colnames(x)[5], "%", "") %>% as.numeric()
)

# -------------------------------------------------------------------------
# Add Souporcell data
sapply(
  setNames(names(Souporcell_Clusters), nm = names(Souporcell_Clusters)),
  function(x){
    message(x)
    tmp <- Souporcell_Clusters[[x]]
    all_found <- all(rownames(FANS[[x]]@meta.data) %in% tmp$barcode)
    FANS[[x]]$Souporcell_Cluster <<- tmp$assignment[match(rownames(FANS[[x]]@meta.data), tmp$barcode)]
    FANS[[x]]$Souporcell_Status  <<- tmp$status[match(rownames(FANS[[x]]@meta.data), tmp$barcode)]
    FANS[[x]]$DoubleCheck.tmp    <<- tmp$barcode[match(rownames(FANS[[x]]@meta.data), tmp$barcode)]
    return(all_found)
  }
) |> all()

# Check correct assignment
sapply(
  FANS,
  function(x){
    if (any(colnames(x@meta.data) == "DoubleCheck.tmp")) {
      return(all(rownames(x@meta.data) == x$DoubleCheck.tmp))
    } else { return(TRUE) }
  }
) |> all()

FANS <- lapply(
  FANS,
  function(x){
    if ("DoubleCheck.tmp" %in% colnames(x@meta.data)) {
      x@meta.data <- x@meta.data %>% dplyr::select(-DoubleCheck.tmp)
    } else {
      x$Souporcell_Cluster <- NA
      x$Souporcell_Status  <- NA
    }
    return(x)
  }
)

Samples$Souporcell_PctAmbientRNA <- Souporcell_AmbientRNA[match(Samples$ID, names(Souporcell_AmbientRNA))]

# -------------------------------------------------------------------------
# Save
qs_save(
  FANS,
  file.path(SeuratObjects_dir, "FANS_Single_SeuratObjects_Filtered_scDblFinder_Marked_Genotypes_Assigned.qs2")
)

qs_save(Samples, Samples_path)
