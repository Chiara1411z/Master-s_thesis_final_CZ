library(SoupX)
library(Seurat)
library(data.table)
library(qs2)

# Input and Output
# Input is Output from CellRanger: raw_feature_bc_matrix/, filtered_feature_bc_matrix/
input_base  <- "/.../"
# Output: Corrected xount matrix
output_base <- "/.../"


qs_path  <- file.path(output_base, "SoupX_Objects.qs2")
csv_path <- file.path(output_base, "SoupX_Contamination_Fractions.csv")

# Relaxed thresholds for samples will few cells
relaxed_samples     <- c("FC51_GEX", "FC55_GEX")
manual_fallback_rho <- 0.015

# ---------------------------------------------------------------------------
# Detect sample folders
sample_ids <- list.files(input_base)
sample_ids <- sample_ids[grepl("^FC", sample_ids)]
cat("Found", length(sample_ids), "valid sample folders:\n")
print(sample_ids)

# ---------------------------------------------------------------------------
# SoupX RNA count matrix correction
SoupX_Objects  <- list()
failed_samples <- character()

for (Sample in sample_ids) {
  
  cat("\n---\nProcessing Sample:", Sample, "\n")
  
  cellranger_outs_raw      <- file.path(input_base, Sample, "outs", "raw_feature_bc_matrix")
  cellranger_outs_filtered <- file.path(input_base, Sample, "outs", "filtered_feature_bc_matrix")
  clusters_path            <- file.path(input_base, Sample, "outs", "analysis",
                                        "clustering", "gene_expression_graphclust", "clusters.csv")
  
  # Prechecks
  if (!dir.exists(cellranger_outs_raw)) {
    warning(paste("Raw matrix not found for", Sample, "- skipping."))
    failed_samples <- c(failed_samples, Sample)
    next
  }
  if (!dir.exists(cellranger_outs_filtered)) {
    warning(paste("Filtered matrix not found for", Sample, "- skipping."))
    failed_samples <- c(failed_samples, Sample)
    next
  }
  if (!file.exists(clusters_path)) {
    warning(paste("clusters.csv not found for", Sample,
                  "(check whether secondary analysis ran) - skipping."))
    failed_samples <- c(failed_samples, Sample)
    next
  }
  
  result <- tryCatch({
    
    raw_data      <- Read10X(cellranger_outs_raw)
    filtered_data <- Read10X(cellranger_outs_filtered)
    clusters      <- fread(clusters_path, data.table = FALSE)
    
    sc <- SoupChannel(raw_data, filtered_data)
    sc <- setClusters(sc, clusters = setNames(clusters$Cluster, nm = clusters$Barcode))
    
    if (Sample %in% relaxed_samples) {
      
      sc_fit <- tryCatch({
        autoEstCont(sc, forceAccept = TRUE, tfidfMin = 0.5, soupQuantile = 0.5)
      }, error = function(e) {
        message("Relaxed autoEstCont failed for ", Sample, ": ", conditionMessage(e))
        NULL
      })
      
      if (!is.null(sc_fit)) {
        sc <- sc_fit
        cat("  (relaxed thresholds used)\n")
      } else {
        sc <- setContaminationFraction(sc, rho = manual_fallback_rho)
        cat("  (manual fallback used) rho =", manual_fallback_rho,
            "- note as limitation in methods\n")
      }
      
    } else {
      sc <- autoEstCont(sc, forceAccept = TRUE)
    }
    
    out <- adjustCounts(sc, roundToInt = TRUE)
    
    # Rho formula
    list(
      SoupChannel    = sc,
      AdjustedCounts = out,
      RNA_SoupX_Rho  = 1 - (sum(out) / sum(filtered_data))
    )
    
  }, error = function(e) {
    warning(paste("SoupX failed for", Sample, ":", conditionMessage(e)))
    NULL
  })
  
  if (is.null(result)) {
    failed_samples <- c(failed_samples, Sample)
    next
  }
  
  SoupX_Objects[[Sample]] <- result
  cat("  rho =", round(result$RNA_SoupX_Rho, 4),
      "| nCells =", ncol(result$AdjustedCounts), "\n")
  
  rm(result)
  gc()
}

cat("\n\nDone.", length(SoupX_Objects), "of", length(sample_ids), "samples corrected.\n")
if (length(failed_samples) > 0) {
  cat("Failed / skipped samples:\n")
  print(failed_samples)
}

# ---------------------------------------------------------------------------
# Sample contamination fractions
soup_fractions_df <- data.frame(
  Sample        = names(SoupX_Objects),
  RNA_SoupX_Rho = sapply(SoupX_Objects, function(x) x$RNA_SoupX_Rho),
  nCells        = sapply(SoupX_Objects, function(x) ncol(x$AdjustedCounts))
)
print(soup_fractions_df)


# Save data 
qs_save(SoupX_Objects, qs_path)
fwrite(soup_fractions_df, csv_path)
