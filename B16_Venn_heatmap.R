library(tidyverse)
library(ComplexHeatmap)
library(ggVennDiagram)
library(patchwork)
library(grid)

output_base <- "/.../"
Plots_dir   <- file.path(output_base, "Plots")


# -------------------------------------------------------------------------
# 1. Heatmap as a patchwork element
col_labels <- gsub("_", " ", colnames(heatmap_mat))

ht_fig <- Heatmap(
  heatmap_mat,
  name            = "avg_log2FC",
  col             = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
  na_col          = "grey80",
  cluster_rows    = TRUE,
  cluster_columns = FALSE,
  column_labels   = col_labels,
  column_title    = NULL     
)

pA <- wrap_elements(full = grid.grabExpr({
  # extra space on top for the title
  ComplexHeatmap::draw(ht_fig, padding = unit(c(2, 2, 12, 2), "mm"))
  grid.text("Candidate gene expression across neuronal subtypes",
            x = unit(12, "mm"), y = unit(1, "npc") - unit(3, "mm"),
            just = c("left", "top"),
            gp = gpar(fontface = "bold", fontsize = 14))
}))

# -------------------------------------------------------------------------
# 2. Venn diagrams
get_sig_genes <- function(df) {
  if (is.null(df)) return(character(0))
  df %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
    dplyr::pull(gene)
}

excitatory_pairs <- list(
  THEMIS    = list(ALS = deg_als_themis, FTD = deg_ftd_themis),
  RORB      = list(ALS = deg_als_rorb,   FTD = deg_ftd_rorb),
  LINC00507 = list(ALS = deg_als_linc,   FTD = deg_ftd_linc),
  FEZF2     = list(ALS = deg_als_fezf2,  FTD = deg_ftd_fezf2)
)

venn_list <- lapply(names(excitatory_pairs), function(cluster_name) {
  venn_input <- list(
    "ALS"     = get_sig_genes(excitatory_pairs[[cluster_name]]$ALS),
    "ALS-FTD" = get_sig_genes(excitatory_pairs[[cluster_name]]$FTD)
  )
  ggVennDiagram(venn_input, label_alpha = 0) +
    scale_fill_gradient(low = "white", high = "#4682B4") +
    coord_flip() +
    scale_x_continuous(expand = expansion(mult = 0.25)) +
    scale_y_continuous(expand = expansion(mult = 0.25)) +
    labs(title = cluster_name) +
    theme(
      legend.position = "none",
      plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
      aspect.ratio = 0.75
    )
})
names(venn_list) <- names(excitatory_pairs)

# -------------------------------------------------------------------------
# 3. Assemble
design <- "
AAAA
BBCC
DDEE
"

fig <- wrap_plots(
  A = pA,
  B = venn_list[["THEMIS"]],
  C = venn_list[["RORB"]],
  D = venn_list[["LINC00507"]],
  E = venn_list[["FEZF2"]],
  design = design
) +
  plot_layout(heights = c(1.1, 1, 1)) +
  plot_annotation(
    title      = "DEG comparison of ALS and ALS-FTD",
    tag_levels = "A",
    theme = theme(
      plot.tag   = element_text(size = 22, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

ggsave(file.path(Plots_dir, "14_Figure_DEG_comparison_heatmap_venn.svg"),
       fig, width = 13, height = 17, device = "svg")
