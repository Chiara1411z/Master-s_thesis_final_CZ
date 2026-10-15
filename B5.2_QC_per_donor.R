library(Seurat)
library(data.table)
library(qs2)
library(tidyverse)

# Input and Output 
# Input is single Seurat object per donor
output_base  <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Samples_path <- file.path(output_base, "Samples", "Samples.qs2")

FANS_Seurats_Raw <- qs_read(file.path(SeuratObjects_dir, "FANS_Single_SeuratObjects_Raw.qs2"))
Samples <- qs_read(Samples_path)

# -------------------------------------------------------------------------
# AUX functions
listCols <- function(m){
  lapply(seq_len(ncol(m)), function(j) {
    idx_start <- m@p[j] + 1
    idx_end <- m@p[j + 1]
    if (idx_start <= idx_end) m@x[idx_start:idx_end] else numeric(0)
  })
}

percentage_in_top_Pct <- function(Var, Pct=NULL, N=NULL){
  sum(sort(Var, decreasing = TRUE)[1:if(is.null(N)) ceiling(Pct*length(Var)) else N], na.rm = TRUE) /
    sum(Var, na.rm = TRUE)
}

filter_by_MADs <- function(Var, nMads){
  MAD <- mad(Var, center = median(Var, na.rm = TRUE), constant = 1, na.rm = TRUE)
  Var < median(Var, na.rm = TRUE) - nMads*MAD | Var > median(Var, na.rm = TRUE) + nMads*MAD
}

# -------------------------------------------------------------------------
# Calculate QC metrics
FANS_Seurats <- lapply(FANS_Seurats_Raw, function(x){
  x$percent.mt   <- PercentageFeatureSet(x, assay = "RNA", pattern = "^MT-")
  x$percent.ribo <- PercentageFeatureSet(x, assay = "RNA", pattern = "^RP[SL]")
  x$percent.soup <- 1 - (colSums(x@assays$RNA$counts, na.rm = TRUE) / colSums(x@assays$RNA_Uncorrected$counts, na.rm = TRUE))
  x$percent.mt_withSoup   <- PercentageFeatureSet(x, assay = "RNA_Uncorrected", pattern = "^MT-")
  x$percent.ribo_withSoup <- PercentageFeatureSet(x, assay = "RNA_Uncorrected", pattern = "^RP[SL]")
  x$Features_by_Counts <- x$nFeature_RNA / x$nCount_RNA
  x$nCount_RNA_log1p   <- log(x$nCount_RNA + 1)
  x$Pct_Counts_In_Top_20_Genes <- vapply(listCols(x@assays$RNA$counts), percentage_in_top_Pct, N = 20, FUN.VALUE = 0.0)
  
  x$MAD_Filter_Flag <- filter_by_MADs(x$nCount_RNA_log1p, nMads = 5) |
    x$percent.mt_withSoup >= 5 |
    x$percent.ribo >= 5
  return(x)
})

# -------------------------------------------------------------------------
# Filter cells by QC
FANS_Seurats <- lapply(FANS_Seurats, function(x){
  subset(x, subset = MAD_Filter_Flag, invert = TRUE)
})

lapply(FANS_Seurats_Raw, dim) %>% data.frame() %>% t() %>% data.frame() %>%
  rownames_to_column(var = "Sample") -> tmp
Samples$nCells_CR <- tmp$X2[match(Samples$ID, tmp$Sample)]
rm(tmp)

lapply(FANS_Seurats, dim) %>% data.frame() %>% t() %>% data.frame() %>%
  rownames_to_column(var = "Sample") -> tmp
Samples$nCells_QC_filtered <- tmp$X2[match(Samples$ID, tmp$Sample)]
rm(tmp)

# -------------------------------------------------------------------------
# Save
qs_save(FANS_Seurats, file.path(SeuratObjects_dir, "FANS_Single_SeuratObjects_Filtered.qs2"))
qs_save(Samples, Samples_path)
