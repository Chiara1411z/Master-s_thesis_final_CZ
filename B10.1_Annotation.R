library(qs2)
library(tidyverse)
library(Seurat)
library(org.Hs.eg.db)
library(rhdf5)

# Input and Output
# Input is merged object
output_base  <- "/.../"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Annotations_dir   <- file.path(output_base, "Annotations")
Samples_path <- file.path(output_base, "Samples", "Samples.qs2")
# Cellranger Output as Input
Cellranger_Output_Dir <- "/.../"

Samples <- qs_read(Samples_path)
FANS <- qs_read(file.path(SeuratObjects_dir, "FANS_Merged_QCFiltered.qs2"))

# -------------------------------------------------------------------------
# Define AUX functions
tab_intersect <- function(a, b, plot = FALSE){
  res <- c(
    "a"     = as.numeric(table(a %in% b)["FALSE"]),
    "a&b"   = as.numeric(table(a %in% b)["TRUE"]),
    "b"     = as.numeric(table(b %in% a)["FALSE"])
  )
  print(res)
  return(invisible(res))
}

# -------------------------------------------------------------------------
# Feature annotation
any(duplicated(rownames(FANS)) | is.na(rownames(FANS)))

FANS_RNA_Features <- data.frame(ID_10X = rownames(FANS))

# Collect Gene Symbols
all_OrgHs_Symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")

tab_intersect(a = FANS_RNA_Features$ID_10X, b = all_OrgHs_Symbols)

# -------------------------------------------------------------------------
# a: Replace ALIASes in RNA features with gene SYMBOLS
grep("^MT-", FANS_RNA_Features$ID_10X, value = TRUE)

ids_to_check <- FANS_RNA_Features$ID_10X[!FANS_RNA_Features$ID_10X %in% all_OrgHs_Symbols]

ids_to_check <- ifelse(
  str_detect(ids_to_check, "^MT-"),
  str_replace(ids_to_check, "^MT-", "MT"),
  ids_to_check
)

tmp <- select(
  org.Hs.eg.db,
  keys = ids_to_check,
  columns = "SYMBOL",
  keytype = "ALIAS",
  multiVals = 'list'
)

tmp <- tmp[!is.na(tmp$SYMBOL), ]
any(is.na(tmp$ALIAS))
tmp <- unique(tmp)

table(duplicated(tmp$ALIAS))
table(duplicated(tmp$SYMBOL))

View(tmp[tmp$ALIAS %in% tmp$ALIAS[duplicated(tmp$ALIAS)], ])
tmp$ALIAS[duplicated(tmp$ALIAS)]

FANS_RNA_Features$SYMBOL <- str_replace_all(FANS_RNA_Features$ID_10X, "MT-", "MT")
FANS_RNA_Features$SYMBOL[FANS_RNA_Features$SYMBOL %in% tmp$ALIAS] <- sapply(
  tmp$SYMBOL[match(
    FANS_RNA_Features$SYMBOL[FANS_RNA_Features$SYMBOL %in% tmp$ALIAS], tmp$ALIAS
  )],
  function(x) x[1]
)

FANS_RNA_Features$SYMBOL[which(!FANS_RNA_Features$SYMBOL %in% keys(org.Hs.eg.db, keytype = "SYMBOL"))] <- NA
table(is.na(FANS_RNA_Features$SYMBOL))


# b: Supplement with Ensembl-ID-based resolution
h5_path <- file.path(Cellranger_Output_Dir, "...", "outs", "filtered_feature_bc_matrix.h5")
feature_ids     <- h5read(h5_path, "matrix/features/id")
feature_symbols <- h5read(h5_path, "matrix/features/name")

feature_symbols_unique <- make.unique(feature_symbols)
stopifnot(all(feature_symbols_unique == FANS_RNA_Features$ID_10X))

id_symbol_map <- data.frame(ID_10X = feature_symbols_unique, ENSEMBL_CellRanger = feature_ids)

CHRs_by_ensembl <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = id_symbol_map$ENSEMBL_CellRanger,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENETYPE", "CHR")
)

nrow(id_symbol_map); nrow(CHRs_by_ensembl); length(unique(CHRs_by_ensembl$ENSEMBL))

CHRs_by_ensembl_dedup <- CHRs_by_ensembl %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

FANS_RNA_Features_v2 <- id_symbol_map %>%
  left_join(CHRs_by_ensembl_dedup, by = c("ENSEMBL_CellRanger" = "ENSEMBL"))

nrow(FANS_RNA_Features_v2) 
table(is.na(FANS_RNA_Features_v2$SYMBOL)) 

# c: Combine -- prefer alias-route symbol, fall back to Ensembl-route where alias failed
FANS_RNA_Features <- FANS_RNA_Features %>%
  dplyr::left_join(
    FANS_RNA_Features_v2 %>% dplyr::select(ID_10X, SYMBOL_ensembl = SYMBOL),
    by = "ID_10X"
  ) %>%
  dplyr::mutate(SYMBOL = coalesce(SYMBOL, SYMBOL_ensembl)) %>%
  dplyr::select(-SYMBOL_ensembl)

table(is.na(FANS_RNA_Features$SYMBOL)) 

# -------------------------------------------------------------------------
# Mark GenBank accessions
FANS_RNA_Features$Accession <- "NotAccession"
FANS_RNA_Features$Accession[grep("\\.", FANS_RNA_Features$ID_10X)] <- "Accession"

View(FANS_RNA_Features[FANS_RNA_Features$Accession == "Accession" & !is.na(FANS_RNA_Features$SYMBOL), ])
FANS_RNA_Features$Accession[FANS_RNA_Features$Accession == "Accession" & !is.na(FANS_RNA_Features$SYMBOL)] <- "NotAccession"

table(is.na(FANS_RNA_Features$SYMBOL), FANS_RNA_Features$Accession)

View(FANS_RNA_Features[FANS_RNA_Features$Accession == "NotAccession" & is.na(FANS_RNA_Features$SYMBOL), ])

# -------------------------------------------------------------------------
# Annotate genetypes
GENETYPES <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = FANS_RNA_Features$SYMBOL,
  keytype = "SYMBOL",
  columns = "GENETYPE"
)

GENETYPES <- unique(GENETYPES)
table(GENETYPES$GENETYPE)

table(duplicated(GENETYPES$SYMBOL))
View(GENETYPES[GENETYPES$SYMBOL %in% GENETYPES$SYMBOL[duplicated(GENETYPES$SYMBOL)], ])
GENETYPES <- GENETYPES[-which(GENETYPES$SYMBOL %in% GENETYPES$SYMBOL[duplicated(GENETYPES$SYMBOL)] & GENETYPES$GENETYPE == "unknown"), ]
table(duplicated(GENETYPES$SYMBOL))
table(GENETYPES$GENETYPE)

FANS_RNA_Features$GENETYPE <- NA
FANS_RNA_Features$GENETYPE <- GENETYPES$GENETYPE[match(FANS_RNA_Features$SYMBOL, GENETYPES$SYMBOL)]

table(is.na(FANS_RNA_Features$GENETYPE), FANS_RNA_Features$Accession)
FANS_RNA_Features$GENETYPE[grep("MT-", FANS_RNA_Features$ID_10X)] <- "mitochondrial-encoded"
View(FANS_RNA_Features[FANS_RNA_Features$Accession != "Accession" & is.na(FANS_RNA_Features$GENETYPE), ])

FANS_RNA_Features$GENETYPE2 <- FANS_RNA_Features$GENETYPE
FANS_RNA_Features$GENETYPE2[FANS_RNA_Features$Accession == "Accession"] <- "Accession"
FANS_RNA_Features$GENETYPE2[FANS_RNA_Features$Accession != "Accession" & is.na(FANS_RNA_Features$GENETYPE2)] <- "unknown"

tmp2 <- FANS_RNA_Features$GENETYPE
FANS_RNA_Features$GENETYPE <- FANS_RNA_Features$GENETYPE2
FANS_RNA_Features$GENETYPE2 <- tmp2
rm(tmp2)

# -------------------------------------------------------------------------
# Annotate ENSEMBL IDs
ENSEMBLS <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = FANS_RNA_Features$SYMBOL,
  keytype = "SYMBOL",
  columns = "ENSEMBL"
)
ENSEMBLS <- unique(ENSEMBLS)

FANS_RNA_Features_ENSEMBL_ENTREZ <- full_join(
  FANS_RNA_Features, ENSEMBLS,
  by = c("SYMBOL" = "SYMBOL"),
  relationship = "many-to-many"
)

# Annotate ENTREZIDs
ENTREZIDS <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = FANS_RNA_Features_ENSEMBL_ENTREZ$ENSEMBL,
  keytype = "ENSEMBL",
  columns = "ENTREZID"
)
ENTREZIDS <- unique(ENTREZIDS)

FANS_RNA_Features_ENSEMBL_ENTREZ <- full_join(
  FANS_RNA_Features_ENSEMBL_ENTREZ, ENTREZIDS,
  by = c("ENSEMBL" = "ENSEMBL"),
  relationship = "many-to-many"
)

# Annotate Chromosomes
CHRs <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = FANS_RNA_Features$SYMBOL,
  keytype = "SYMBOL",
  columns = "CHR"
)

FANS_RNA_Features$CHR <- CHRs$CHR[match(FANS_RNA_Features$SYMBOL, CHRs$SYMBOL)]

# Save
qs_save(FANS_RNA_Features, file.path(Annotations_dir, "FANS_RNA_Features.qs2"))
qs_save(FANS_RNA_Features_ENSEMBL_ENTREZ, file.path(Annotations_dir, "FANS_RNA_Features_ENSEMBL_ENTREZ.qs2"))
