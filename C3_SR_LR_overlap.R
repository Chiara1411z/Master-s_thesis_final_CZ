library(tidyverse)
library(svglite)
library(writexl)
library(patchwork)

# ------------------------------------------------------------------
# Input and Output
# Input are SR and LR files for QC
sr_dir     <- "/.../"
data_base  <- "/.../"
output_dir <- "/.../"

libraries <- c("FC32", "FC42", paste0("FC", 45:55))

lr_knee_cutoff <- 50  

use_whitelist_files_if_present <- TRUE

# ------------------------------------------------------------------
# Helpers
revcomp <- function(seqs) {
  comp <- chartr("ACGTN", "TGCAN", seqs)
  vapply(strsplit(comp, ""), function(s) paste(rev(s), collapse = ""), character(1))
}

sr_summary <- read.csv(file.path(sr_dir, "All_Cell_Barcodes_Summary.csv")) %>%
  rename_with(tolower)

# ------------------------------------------------------------------
# Per-library analysis
overview <- list()

for (lib in libraries) {
  cat("\n Processing Library:", lib, "\n")
  
  # Library directory (Phase3 subfolder vs. library root)
  lib_dir_phase3 <- file.path(data_base, lib, "Phase3")
  lib_dir <- if (dir.exists(lib_dir_phase3)) lib_dir_phase3 else file.path(data_base, lib)
  
  # Long-read barcode counts
  lr_candidates <- c(file.path(lib_dir, paste0("lr_barcode_counts_", lib, ".txt")),
                     file.path(lib_dir, "lr_barcode.txt"))
  lr_file <- lr_candidates[file.exists(lr_candidates)][1]
  if (is.na(lr_file)) {
    warning("Skipping ", lib, ": no long-read count file in ", lib_dir)
    next
  }
  
  lr_counts <- read.table(lr_file, col.names = c("long_read_counts", "barcode"),
                          stringsAsFactors = FALSE) %>%
    mutate(barcode = revcomp(barcode))
  
  # Short-read rows for this library
  lib_pat <- paste0("^", lib, "(_GEX)?($|[^0-9])")
  lib_rows <- sr_summary %>%
    filter(if_any(any_of(c("sample_id", "cell_type", "cell_barcode", "sample")),
                  ~ grepl(lib_pat, .x, ignore.case = TRUE))) %>%
    mutate(barcode = gsub("-.*", "", cell_barcode))
  
  # SR whitelist
  whitelist_file <- file.path(lib_dir, paste0("whitelist_", lib, ".tsv"))
  if (use_whitelist_files_if_present && file.exists(whitelist_file)) {
    sr_barcodes <- unique(gsub("-.*", "", readLines(whitelist_file)))
    wl_source   <- "whitelist file"
  } else {
    sr_barcodes <- unique(lib_rows$barcode)
    wl_source   <- "summary CSV"
  }
  
  if (length(sr_barcodes) == 0) {
    warning("Skipping ", lib, ": no SR barcodes found (check labels, e.g. ", lib, "_GEX)")
    next
  }
  
  sr_counts <- lib_rows %>% count(barcode, name = "short_read_counts")
  
  # Q1: share of long reads sitting in SR-called cells
  total_lr_reads  <- sum(lr_counts$long_read_counts)
  pct_lr_in_cells <- 100 * sum(lr_counts$long_read_counts[lr_counts$barcode %in% sr_barcodes]) /
    total_lr_reads
  
  # Q2/Q3: SR-LR correlation
  df_merged <- inner_join(sr_counts, lr_counts, by = "barcode")
  rho <- NA_real_
  if (nrow(df_merged) >= 3) {
    rho <- unname(cor.test(log10(df_merged$short_read_counts + 1),
                           log10(df_merged$long_read_counts + 1),
                           method = "spearman", exact = FALSE)$estimate)
  }
  
  # Q4: SR cells with at least one long read
  sr_cells_with_lr <- sum(lr_counts$barcode %in% sr_barcodes & lr_counts$long_read_counts > 0)
  
  # Q5: independent LR cell calling and overlap
  lr_called_barcodes <- lr_counts %>%
    filter(long_read_counts >= lr_knee_cutoff) %>%
    pull(barcode) %>% unique()
  
  n_both    <- length(intersect(sr_barcodes, lr_called_barcodes))
  n_sr_only <- length(setdiff(sr_barcodes, lr_called_barcodes))
  n_lr_only <- length(setdiff(lr_called_barcodes, sr_barcodes))
  
  # Mean / median LR reads per SR-called cell (SR cells without LR reads count as 0)
  lr_sr_cells <- lr_counts$long_read_counts[lr_counts$barcode %in% sr_barcodes]
  n_zero      <- length(sr_barcodes) - length(lr_sr_cells)
  lr_sr_all   <- c(lr_sr_cells, rep(0, n_zero))
  
  mean_lr_sr_cells <- mean(lr_sr_all)
  med_lr_sr_cells  <- median(lr_sr_all)
  
  # Mean LR reads per LR-called cell (>= cutoff), regardless of SR calling
  mean_lr_lr_called <- mean(lr_counts$long_read_counts[lr_counts$long_read_counts >= lr_knee_cutoff])
  
  overview[[lib]] <- tibble(
    library                          = lib,
    sr_whitelist_source              = wl_source,
    total_lr_reads                   = total_lr_reads,
    pct_lr_reads_in_sr_cells         = pct_lr_in_cells,
    n_sr_cells                       = length(sr_barcodes),
    pct_sr_cells_with_any_lr         = 100 * sr_cells_with_lr / length(sr_barcodes),
    mean_lr_reads_per_sr_cell        = mean_lr_sr_cells,
    median_lr_reads_per_sr_cell      = med_lr_sr_cells,
    mean_lr_reads_per_lr_called_cell = mean_lr_lr_called,
    spearman_rho                     = rho,
    n_lr_called                      = length(lr_called_barcodes),
    n_overlap                        = n_both,
    n_sr_only                        = n_sr_only,
    n_lr_only                        = n_lr_only,
    jaccard                          = n_both / (n_both + n_sr_only + n_lr_only),
    pct_sr_cells_in_lr_called        = 100 * n_both / length(sr_barcodes)
  )
}

overview_df <- bind_rows(overview)

# ------------------------------------------------------------------
# Summary table
summary_table <- overview_df %>%
  transmute(
    Library                          = library,
    `SR cells`                       = n_sr_cells,
    `LR barcodes (>= cutoff)`        = n_lr_called,
    Overlap                          = n_overlap,
    `SR only`                        = n_sr_only,
    `LR only`                        = n_lr_only,
    Jaccard                          = round(jaccard, 3),
    `% SR cells in LR-called`        = round(pct_sr_cells_in_lr_called, 1),
    `% SR cells with any LR`         = round(pct_sr_cells_with_any_lr, 1),
    `% LR reads in SR cells`         = round(pct_lr_reads_in_sr_cells, 1),
    `Mean LR reads / SR cell`        = round(mean_lr_reads_per_sr_cell, 1),
    `Median LR reads / SR cell`      = median_lr_reads_per_sr_cell,
    `Mean LR reads / LR-called cell` = round(mean_lr_reads_per_lr_called_cell, 1),
    `Spearman rho (SR vs LR)`        = round(spearman_rho, 3),
    `Total LR reads`                  = total_lr_reads,
    `SR whitelist from`              = sr_whitelist_source
  )

write.csv(summary_table,
          file.path(output_dir, "All_FCs_LR_SR_summary_table.csv"),
          row.names = FALSE)

write_xlsx(summary_table,
           file.path(output_dir, "All_FCs_LR_SR_summary_table.xlsx"))

# ------------------------------------------------------------------
# Plot Venn diagrams
col_lr <- "#C5D3E8"  
col_sr <- "#4F7CB5"

r <- 1 
d <- 0.55 
th0 <- acos(d / r)

lib_levels <- overview_df$library

circle_pts <- function(cx) {
  th <- seq(0, 2 * pi, length.out = 400)
  tibble(x = cx + r * cos(th), y = r * sin(th))
}

circ_lr <- circle_pts(-d)
circ_sr <- circle_pts(d)

lens <- bind_rows(
  tibble(th = seq(-th0, th0, length.out = 100)) %>% mutate(x = -d + r * cos(th), y = r * sin(th)),
  tibble(th = seq(pi - th0, pi + th0, length.out = 100)) %>% mutate(x = d + r * cos(th), y = r * sin(th))
) %>% select(x, y)

# Region labels: count and % of the union
lab_df <- overview_df %>%
  transmute(library, lr_only = n_lr_only, both = n_overlap, sr_only = n_sr_only) %>%
  mutate(total = lr_only + both + sr_only) %>%
  pivot_longer(c(lr_only, both, sr_only), names_to = "region", values_to = "n") %>%
  mutate(
    x     = unname(c(lr_only = -r, both = 0, sr_only = r)[region]),
    label = paste0(n, "\n(", round(100 * n / total), "%)")
  )

# One Venn panel per library
make_venn <- function(lib) {
  labs_lib <- lab_df %>% filter(library == lib)
  
  ggplot() +
    geom_polygon(data = circ_lr, aes(x, y), fill = col_lr, colour = NA) +
    geom_polygon(data = circ_sr, aes(x, y), fill = col_sr, colour = NA) +
    geom_polygon(data = lens,    aes(x, y), fill = "white", colour = NA) +
    geom_path(data = circ_lr, aes(x, y), colour = "black", linewidth = 0.5) +
    geom_path(data = circ_sr, aes(x, y), colour = "black", linewidth = 0.5) +
    geom_text(data = labs_lib, aes(x, 0, label = label),
              family = "serif", size = 4, lineheight = 0.95) +
    annotate("text", x = c(-d, d), y = -1.3,
             label = c("Long-reads", "Short-reads"),
             family = "serif", size = 4) +
    labs(title = lib) +
    coord_fixed(xlim = c(-1.7, 1.7), ylim = c(-1.55, 1.1)) +
    theme_void(base_family = "serif") +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5,
                                    margin = margin(b = 4)))
}

plots <- lapply(lib_levels, make_venn)

p_venn_all <- wrap_plots(plots, ncol = 3) +
  plot_annotation(
    title      = "Barcode overlap between long-read and short-read cell calling",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5,
                                            family = "serif", margin = margin(b = 8)))
  ) &
  theme(
    plot.tag        = element_text(face = "bold", size = 16, family = "serif"),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(6, 6, 6, 6)
  )

n_rows <- ceiling(length(lib_levels) / 3)
fig_w  <- 13
fig_h  <- 3.3 * n_rows + 0.8

ggsave(file.path(output_dir, "All_FCs_venn_diagrams.svg"),
       plot = p_venn_all, width = fig_w, height = fig_h, device = "svg")
