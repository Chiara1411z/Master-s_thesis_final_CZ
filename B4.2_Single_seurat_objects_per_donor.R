library(Seurat)
library(data.table)
library(qs2)
library(tidyverse)
library(ggplot2)
library(svglite)

# Input and Output
# Input is Demultiplexing output
output_base    <- "/.../"
Samples_path   <- file.path(output_base, "Samples", "Samples.qs2")

# -------------------------------------------------------------------------
# Load existing Samples table, merge in any new samples not yet present
new_samples <- data.frame(
  ID          = c(...),
  tdp_status  = c(...), 
  neun_status = c(...)
)

if (file.exists(Samples_path)) {
  Samples <- qs_read(Samples_path)
  missing_new <- new_samples %>% filter(!ID %in% Samples$ID)
  if (nrow(missing_new) > 0) {
    message("Adding new samples to Samples table: ", paste(missing_new$ID, collapse = ", "))
    Samples <- bind_rows(Samples, missing_new)
    qs_save(Samples, Samples_path)
  }
} else {
  Samples <- new_samples
  qs_save(Samples, Samples_path)
}

# Cellranger output and Souporcell output files as input
Cellranger_Output_Dir  <- "/.../"
Souporcell_Output_Dir  <- "/.../"
SeuratObjects_dir      <- file.path(output_base, "SeuratObjects")

# -------------------------------------------------------------------------
# Confirm Souporcell has actually been run for every sample before proceeding
missing_souporcell <- Samples$ID[!file.exists(file.path(Souporcell_Output_Dir, Samples$ID, "clusters.tsv"))]
if (length(missing_souporcell) > 0) {
  stop("Missing Souporcell output (clusters.tsv) for: ", paste(missing_souporcell, collapse = ", "),
       " -- run Souporcell for these samples before continuing.")
}

# -------------------------------------------------------------------------
# Generate one Seurat object per donor per library

FANS_Seurats_Donors <- list()

for (Sample in Samples$ID) {
  
  message(paste0("Sample ", Sample, ": ..."))
  
  filtered_data <- Read10X_h5(
    file.path(Cellranger_Output_Dir, Sample, "outs", "filtered_feature_bc_matrix.h5")
  )
  
  clusters <- fread(
    file.path(Souporcell_Output_Dir, Sample, "clusters.tsv"),
    data.table = FALSE
  )
  
  # Keep only confidently assigned singlets -> drop doublets/unassigned
  clusters_singlet <- clusters %>% filter(status == "singlet")
  
  if (!all(clusters_singlet$barcode %in% colnames(filtered_data))) {
    warning(paste0(
      "Some Souporcell barcodes for ", Sample,
      " were not found in the filtered matrix!"
    ))
  }
  
  donor_ids <- sort(unique(clusters_singlet$assignment))
  
  for (donor in donor_ids) {
    
    donor_label <- paste0(Sample, "_donor", donor)
    message("  -> ", donor_label)
    
    meta_sub <- clusters_singlet %>% filter(assignment == donor)
    donor_barcodes <- intersect(meta_sub$barcode, colnames(filtered_data))
    
    donor_counts <- filtered_data[, donor_barcodes]
    
    Seurat_obj <- CreateSeuratObject(
      counts       = donor_counts,
      assay        = "RNA",
      project      = donor_label,
      min.cells    = 0,
      min.features = 0
    )
    rownames(meta_sub) <- meta_sub$barcode
    Seurat_obj$Souporcell_status     <- meta_sub[colnames(Seurat_obj), "status"]
    Seurat_obj$Souporcell_assignment <- meta_sub[colnames(Seurat_obj), "assignment"]
    Seurat_obj$Library               <- Sample
    
    FANS_Seurats_Donors[[donor_label]] <- Seurat_obj
  }
  message("")
}


# -------------------------------------------------------------------------
# Build a donor-level metadata table

Donors <- data.frame(
  ID      = names(FANS_Seurats_Donors),
  Library = sapply(FANS_Seurats_Donors, function(x) unique(x$Library)),
  nCells_Souporcell = sapply(FANS_Seurats_Donors, ncol)
)

Donors <- Donors %>%
  left_join(Samples %>% dplyr::select(ID, tdp_status, neun_status), by = c("Library" = "ID"))


# -------------------------------------------------------------------------
# Save
qs_save(
  FANS_Seurats_Donors,
  file.path(SeuratObjects_dir, "FANS_Singulatity_Donor_SeuratObjects_Raw.qs2")
)

qs_save(
  Donors,
  file.path(output_base, "Samples", "Singularity_Donors.qs2")
)



# -------------------------------------------------------------------------
View(Donors) 
Donors %>%
  group_by(Library) %>%
  summarise(n_donors = n(), total_cells = sum(nCells_Souporcell), .groups = "drop")

# -------------------------------------------------------------------------
# Plot

Donors$Donor_Group <- as.factor(str_extract(Donors$ID, "[0-9]+$"))

p_donors <- ggplot(Donors, aes(x = Library, y = nCells_Souporcell, fill = Donor_Group)) +
  geom_col(position = "stack") +
  scale_fill_manual(
    name = "Donor",
    values = c(
      "#2ca02c",
      "#9467bd",
      "#e377c2",
      "#ff7f0e"
    )
  ) +
  scale_x_discrete(labels = function(x) str_remove(x, "_GEX$")) +
  theme_minimal() +
  labs(y = "Cells per donor", x = "Library") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# Save
ggsave(
  filename = file.path(output_base, "Plots", "Singularity_cells_per_donor.svg"), 
  plot = p_donors, 
  device = "svg", 
  width = 8, 
  height = 6, 
  dpi = 300
)


# -------------------------------------------------------------------------
# Check for doublet, singlet, unassigned
clusters <- fread(file.path(Souporcell_Output_Dir, "FC31", "clusters.tsv"), data.table = FALSE)

table(clusters$status)   # singlet / doublet / unassigned counts

ggplot(clusters, aes(x = log_prob_singleton, y = log_prob_doublet, color = status)) +
  geom_point(alpha = 0.4) +
  theme_minimal()

