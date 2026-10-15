library(Seurat)
library(data.table)
library(qs2)
library(tidyverse)

# Input and Output
# Input: filtered_feature_bc_matrix/, Sample table
output_base  <- "/.../"
Samples_path <- file.path(output_base, "Samples", "Samples.qs2")

if (file.exists(Samples_path)) {
  Samples <- qs_read(Samples_path)
  
} else {
  Samples <- data.frame(
    ID          = c(...),
    tdp_status  = c(...),    
    neun_status = c(...)
  )
  qs_save(Samples, Samples_path)
}

SoupX_Objects <- qs_read(
  file.path(output_base, "SoupX", "SoupX_Objects.qs2")
)

Cellranger_Output_Dir <- "/.../"

SeuratObjects_dir <- file.path(output_base, "SeuratObjects")


# -------------------------------------------------------------------------
# Generate Seurat objects
FANS_Seurats <- lapply(
  setNames(Samples$ID, nm = Samples$ID),
  function(Sample) {
    
    message(
      paste0(
        "Sample ",
        Sample,
        ": ..."
      )
    )
    
    filtered_data <- Read10X_h5(
      file.path(Cellranger_Output_Dir, Sample, "outs", "filtered_feature_bc_matrix.h5")
    )
    
    corrected <- SoupX_Objects[[Sample]]$AdjustedCounts
    
    if (is.null(corrected)) {
      stop(paste0(
        "No SoupX-corrected counts found for ", Sample,
      ))
    }
    
    corrected <- corrected[
      , which(colnames(corrected) %in% colnames(filtered_data))
    ]
    
    if (all(colnames(corrected) == colnames(filtered_data))) {
      message("Corrected and uncorrected data sorted concordantly!")
    } else {
      message("Corrected and uncorrected data NOT sorted concordantly! Sorting corrected data... ")
      corrected <- corrected[, match(colnames(filtered_data), colnames(corrected))]
    }
    
    Seurat <- CreateSeuratObject(
      counts = corrected,
      assay = "RNA",
      project = Sample,
      min.cells = 0,
      min.features = 0
    )
    
    Seurat[["RNA_Uncorrected"]] <- CreateAssayObject(
      counts = filtered_data,
      min.cells = 0,
      min.features = 0
    )
    
    if (all(rownames(Seurat@assays$RNA) == rownames(Seurat@assays$RNA_Uncorrected)) &
        all(colnames(Seurat@assays$RNA) == colnames(Seurat@assays$RNA_Uncorrected))) {
      message("'RNA' and 'RNA_Uncorrected' Assays concordant! ")
    } else {
      warning("Warning! 'RNA' and 'RNA_Uncorrected' Assays not concordant! ")
    }
    
    message("")
    
    return(Seurat)
    
  }
)


# -------------------------------------------------------------------------
# Calculate SoupX Ambient RNA
SoupX_Pct_Ambient <- sapply(
  FANS_Seurats,
  function(x) {
    (1 - (sum(x@assays$RNA$counts) / sum(x@assays$RNA_Uncorrected$counts))) * 100
  }
)
all(
  names(SoupX_Pct_Ambient) ==
    Samples$ID
)

Samples$SoupX_Pct_Ambient <- SoupX_Pct_Ambient

# -------------------------------------------------------------------------
# Add number of called cells to Samples data

all(names(FANS_Seurats) == Samples$ID)
nCells <- sapply(
  FANS_Seurats,
  function(x) {
    return(
      dim(x)[2]
    )
  }
)
Samples$nCells_CR <- nCells[match(Samples$ID, names(nCells))]


# -------------------------------------------------------------------------
# Save
qs_save(
  FANS_Seurats,
  file.path(SeuratObjects_dir, "FANS_Single_SeuratObjects_Raw.qs2")
)

qs_save(
  Samples,
  Samples_path
)
