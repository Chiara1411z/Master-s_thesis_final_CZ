library(qs2)
library(tidyverse)
library(Seurat)
library(glmGamPoi)
library(org.Hs.eg.db)
library(writexl)
library(rhdf5)

# Input and Output
# Input is filtered object
output_base       <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Plots_dir         <- file.path(output_base, "Plots")
Annotations_dir   <- file.path(output_base, "Annotations")
Donors_path       <- file.path(output_base, "Samples", "Singularity_Donors.qs2")

# Cellranger output as input
Cellranger_Output_Dir <- "/.../"

FANS_Donors <- qs_read(file.path(SeuratObjects_dir, "FANS_Donor_Merged_QCFiltered.qs2"))
Donors      <- qs_read(Donors_path)

# -------------------------------------------------------------------------
# Cell Counts per Donor
print("Cell counts per donor (Pre-filtering):")
print(table(FANS_Donors$ID))
write_xlsx(as.data.frame(table(FANS_Donors$ID)), "FANS_Donors_ID_Counts.xlsx")

p_donor_counts <- ggplot(FANS_Donors@meta.data, aes(x = ID, fill = ID)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 4) +
  theme_minimal() +
  labs(
    title = "Number of Cells per Donor (Pre-Filter)",
    x = "Donor ID",
    y = "Total Cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(file.path(Plots_dir, "Cells_per_donor.svg"), p_donor_counts, device = "svg", width = 8, height = 6, dpi = 300)

# -------------------------------------------------------------------------
# Filter Out Donors with Fewer Than 50 Cells
donor_counts    <- table(FANS_Donors$ID)
low_cell_donors <- names(donor_counts[donor_counts < 50])

if (length(low_cell_donors) > 0) {
  message("Excluding donors with < 50 cells: ", paste(low_cell_donors, collapse = ", "))
  valid_donors <- names(donor_counts[donor_counts >= 50])
  FANS_Donors  <- subset(FANS_Donors, subset = ID %in% valid_donors)
}

# -------------------------------------------------------------------------
# Feature Annotation for FANS_Donors - Route A: alias/symbol matching
FANS_RNA_Features <- data.frame(ID_10X = rownames(FANS_Donors))
all_OrgHs_Symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")

ids_to_check <- FANS_RNA_Features$ID_10X[!FANS_RNA_Features$ID_10X %in% all_OrgHs_Symbols]
ids_to_check <- ifelse(str_detect(ids_to_check, "^MT-"), str_replace(ids_to_check, "^MT-", "MT"), ids_to_check)

tmp <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = ids_to_check, columns = "SYMBOL", keytype = "ALIAS", multiVals = 'list')
)
tmp <- unique(tmp[!is.na(tmp$SYMBOL), ])

FANS_RNA_Features$SYMBOL <- str_replace_all(FANS_RNA_Features$ID_10X, "MT-", "MT")
FANS_RNA_Features$SYMBOL[FANS_RNA_Features$SYMBOL %in% tmp$ALIAS] <- sapply(
  tmp$SYMBOL[match(FANS_RNA_Features$SYMBOL[FANS_RNA_Features$SYMBOL %in% tmp$ALIAS], tmp$ALIAS)],
  function(x) x[1]
)
FANS_RNA_Features$SYMBOL[which(!FANS_RNA_Features$SYMBOL %in% keys(org.Hs.eg.db, keytype = "SYMBOL"))] <- NA

message("Alias-route resolved: ", sum(!is.na(FANS_RNA_Features$SYMBOL)), " / ", nrow(FANS_RNA_Features))

# -------------------------------------------------------------------------
# Feature Annotation - Route B: Ensembl-ID matching (supplements Route A)
h5_path <- file.path(Cellranger_Output_Dir, "FC31", "outs", "filtered_feature_bc_matrix.h5")
feature_ids     <- h5read(h5_path, "matrix/features/id")
feature_symbols <- h5read(h5_path, "matrix/features/name")


feature_symbols_unique <- make.unique(feature_symbols)

if (!all(feature_symbols_unique == FANS_RNA_Features$ID_10X)) {
  stop("Feature ID mismatch between h5 reference and FANS_Donors rownames -- ",
       "check whether FANS_Donors has a reduced feature set (e.g. from min.cells/",
       "min.features filtering) before proceeding.")
}

id_symbol_map <- data.frame(ID_10X = feature_symbols_unique, ENSEMBL_CellRanger = feature_ids)

CHRs_by_ensembl <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = id_symbol_map$ENSEMBL_CellRanger,
                        keytype = "ENSEMBL", columns = c("SYMBOL", "GENETYPE", "CHR"))
)

CHRs_by_ensembl_dedup <- CHRs_by_ensembl %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

FANS_RNA_Features_v2 <- id_symbol_map %>%
  dplyr::left_join(CHRs_by_ensembl_dedup, by = c("ENSEMBL_CellRanger" = "ENSEMBL"))

stopifnot(nrow(FANS_RNA_Features_v2) == nrow(FANS_RNA_Features))
message("Ensembl-route resolved: ", sum(!is.na(FANS_RNA_Features_v2$SYMBOL)), " / ", nrow(FANS_RNA_Features_v2))

# -------------------------------------------------------------------------
# Combine both routes -- prefer alias-route symbol, fall back to Ensembl-route
FANS_RNA_Features <- FANS_RNA_Features %>%
  dplyr::left_join(
    FANS_RNA_Features_v2 %>% dplyr::select(ID_10X, SYMBOL_ensembl = SYMBOL),
    by = "ID_10X"
  ) %>%
  dplyr::mutate(SYMBOL = coalesce(SYMBOL, SYMBOL_ensembl)) %>%
  dplyr::select(-SYMBOL_ensembl)

message("Combined resolved: ", sum(!is.na(FANS_RNA_Features$SYMBOL)), " / ", nrow(FANS_RNA_Features))

# -------------------------------------------------------------------------
# Mark GenBank Accessions
FANS_RNA_Features$Accession <- ifelse(grepl("\\.", FANS_RNA_Features$ID_10X) & is.na(FANS_RNA_Features$SYMBOL), "Accession", "NotAccession")

# Annotate Genetypes
GENETYPES <- unique(suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = na.omit(FANS_RNA_Features$SYMBOL), keytype = "SYMBOL", columns = "GENETYPE")
))
GENETYPES <- GENETYPES[!(GENETYPES$SYMBOL %in% GENETYPES$SYMBOL[duplicated(GENETYPES$SYMBOL)] & GENETYPES$GENETYPE == "unknown"), ]

FANS_RNA_Features$GENETYPE <- GENETYPES$GENETYPE[match(FANS_RNA_Features$SYMBOL, GENETYPES$SYMBOL)]
FANS_RNA_Features$GENETYPE[grep("^MT-", FANS_RNA_Features$ID_10X)] <- "mitochondrial-encoded"
FANS_RNA_Features$GENETYPE <- ifelse(FANS_RNA_Features$Accession == "Accession", "Accession",
                                     ifelse(is.na(FANS_RNA_Features$GENETYPE), "unknown", FANS_RNA_Features$GENETYPE))

# Annotate ENSEMBL, ENTREZ, and Chromosomes
ENSEMBLS <- unique(suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = na.omit(FANS_RNA_Features$SYMBOL), keytype = "SYMBOL", columns = "ENSEMBL")
))
FANS_RNA_Features_ENSEMBL_ENTREZ <- full_join(FANS_RNA_Features, ENSEMBLS, by = "SYMBOL", relationship = "many-to-many")

ENTREZIDS <- unique(suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = na.omit(FANS_RNA_Features_ENSEMBL_ENTREZ$ENSEMBL), keytype = "ENSEMBL", columns = "ENTREZID")
))
FANS_RNA_Features_ENSEMBL_ENTREZ <- full_join(FANS_RNA_Features_ENSEMBL_ENTREZ, ENTREZIDS, by = "ENSEMBL", relationship = "many-to-many")

CHRs <- unique(suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = na.omit(FANS_RNA_Features$SYMBOL), keytype = "SYMBOL", columns = "CHR")
))
FANS_RNA_Features$CHR <- CHRs$CHR[match(FANS_RNA_Features$SYMBOL, CHRs$SYMBOL)]

# Add to metadata
FANS_Donors[["RNA"]]@meta.data <- left_join(
  data.frame(ID_10X = rownames(FANS_Donors[["RNA"]])),
  FANS_RNA_Features,
  by = "ID_10X"
)
rownames(FANS_Donors[["RNA"]]@meta.data) <- FANS_Donors[["RNA"]]@meta.data$ID_10X

# Save the filtered & annotated Seurat object
qs_save(FANS_Donors, file.path(SeuratObjects_dir, "FANS_Donor_Annotated.qs2"))

# Save Annotation Tables
qs_save(FANS_RNA_Features, file.path(Annotations_dir, "FANS_Donor_RNA_Features.qs2"))
qs_save(FANS_RNA_Features_ENSEMBL_ENTREZ, file.path(Annotations_dir, "FANS_Donor_RNA_Features_ENSEMBL_ENTREZ.qs2"))

# -------------------------------------------------------------------------
# Save Feature Annotation Diagnostic Plots
p_genetype <- ggplot(FANS_RNA_Features, aes(x = fct_infreq(GENETYPE))) +
  geom_bar(fill = "pink") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.3, size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Number of features", title = "Gene type composition")
ggsave(file.path(Plots_dir, "10a_Gene_type_composition_Donors.svg"), plot = p_genetype, device = "svg", width = 8, height = 6, dpi = 300)

p_resolution <- FANS_RNA_Features %>%
  mutate(Resolved = case_when(
    Accession == "Accession" ~ "Accession-style (unresolved)",
    is.na(SYMBOL) ~ "No valid symbol",
    TRUE ~ "Resolved symbol"
  )) %>%
  dplyr::count(Resolved) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ggplot(aes(x = "", y = n, fill = Resolved)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(n, "\n(", pct, "%)")), position = position_stack(vjust = 0.5), size = 3) +
  theme_void() +
  labs(title = "Feature resolution status", fill = NULL)
ggsave(file.path(Plots_dir, "10b_Feature_resolution_status_Donors.svg"), plot = p_resolution, device = "svg", width = 8, height = 6, dpi = 300)

p_chromosomes <- FANS_RNA_Features %>%
  filter(!is.na(CHR)) %>%
  mutate(CHR = factor(CHR, levels = c(as.character(1:22), "X", "Y", "MT"))) %>%
  ggplot(aes(x = CHR)) +
  geom_bar(fill = "orange") +
  theme_classic() +
  labs(x = "Chromosome", y = "Number of features", title = "Feature distribution across chromosomes")
ggsave(file.path(Plots_dir, "10c_Chromosome_distribution_Donors.svg"), plot = p_chromosomes, device = "svg", width = 10, height = 6, dpi = 300)

