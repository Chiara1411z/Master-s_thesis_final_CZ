library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)

# Input and Output
output_base <- "/.../"
csv_dir     <- output_base
Plots_dir   <- file.path(output_base, "Plots")

# -------------------------------------------------------------------------
# Shared theme (matches the other thesis figures)
theme_thesis <- theme(
  plot.title           = element_text(hjust = 0.5, face = "bold", size = 11),
  axis.title.x         = element_text(face = "bold", size = 10),
  axis.title.y         = element_text(face = "bold", size = 10),
  axis.text            = element_text(color = "black", face = "bold", size = 9),
  axis.line            = element_line(linewidth = 1),
  strip.text           = element_text(face = "bold", size = 11),
  legend.title         = element_text(face = "bold", size = 10),
  legend.text          = element_text(size = 9),
  legend.justification = "top"
)

# -------------------------------------------------------------------------
# 1. Volcano plot function
plot_volcano_fig <- function(deg_df, title, y_quantile = 0.99) {
  
  plot_data <- deg_df %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(!is.na(p_val_adj)) %>%
    dplyr::mutate(
      sig = dplyr::case_when(
        p_val_adj < 0.05 & avg_log2FC >  0.5 ~ "Increased",
        p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Decreased",
        TRUE ~ "Not Significant"
      ),
      sig       = factor(sig, levels = c("Not Significant", "Decreased", "Increased")),
      underflow = p_val_adj == 0,      
      nlogp     = -log10(p_val_adj)     
    )
  
  n_inc <- sum(plot_data$sig == "Increased") 
  n_dec <- sum(plot_data$sig == "Decreased")
  
  # y-range from genes with a real p-value only
  real_vals <- plot_data$nlogp[!plot_data$underflow]
  y_cap <- quantile(real_vals, y_quantile, na.rm = TRUE)
  y_top <- y_cap * 1.25    
  
  # underflow genes are placed far above the visible range
  plot_data$nlogp[plot_data$underflow] <- y_top * 10
  
  # draw significant genes on top of the grey mass
  plot_data <- plot_data %>% dplyr::arrange(sig)
  
  lbl_dec <- paste0("atop(bolditalic('Decreased in TDP-43'^Low), bolditalic('n = ", n_dec, "'))")
  lbl_inc <- paste0("atop(bolditalic('Increased in TDP-43'^Low), bolditalic('n = ", n_inc, "'))")
  
  xmax <- max(abs(plot_data$avg_log2FC)) * 1.05
  
  ggplot(plot_data, aes(x = avg_log2FC, y = nlogp, color = sig)) +
    geom_point(alpha = 0.6, size = 0.7, stroke = 0) +
    scale_color_manual(values = c("Increased"       = "#D53E4F",
                                  "Decreased"       = "#3288BD",
                                  "Not Significant" = "grey35")) +
    geom_vline(xintercept = c(-0.5, 0.5), color = "grey50", linewidth = 0.4, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05),  color = "grey50", linewidth = 0.4, linetype = "dashed") +
    coord_cartesian(xlim = c(-xmax, xmax),
                    ylim = c(-y_cap * 0.02, y_top),
                    expand = FALSE) +
    annotate("text", x = -Inf, y = Inf, label = lbl_dec, parse = TRUE,
             color = "#3288BD", hjust = -0.05, vjust = 1.1, size = 3.5) +
    annotate("text", x =  Inf, y = Inf, label = lbl_inc, parse = TRUE,
             color = "#D53E4F", hjust = 1.05, vjust = 1.1, size = 3.5) +
    theme_classic() +
    theme_thesis +
    labs(title = title, x = "Log2 FC [mean]", y = "FDR [q-value, -log10]") +
    theme(legend.position = "none")
}

# -------------------------------------------------------------------------
# 2. GO ORA function that RETURNS the plot
make_go_plot <- function(deg_files, title_text, top_n = 5) {
  gene_lists <- list()
  
  for (subtype_name in names(deg_files)) {
    full_path <- file.path(csv_dir, deg_files[[subtype_name]])
    if (file.exists(full_path)) {
      df <- read.csv(full_path)
      gene_col <- intersect(c("Gene", "gene", "X"), colnames(df))[1]
      
      up_genes <- df %>%
        dplyr::filter(!is.na(p_val_adj) & p_val_adj < 0.05 & avg_log2FC > 0.5) %>%
        dplyr::pull(!!rlang::sym(gene_col)) %>% as.character()
      down_genes <- df %>%
        dplyr::filter(!is.na(p_val_adj) & p_val_adj < 0.05 & avg_log2FC < -0.5) %>%
        dplyr::pull(!!rlang::sym(gene_col)) %>% as.character()
      
      gene_lists[[paste0(subtype_name, "||Upregulated")]]   <- up_genes
      gene_lists[[paste0(subtype_name, "||Downregulated")]] <- down_genes
    }
  }
  
  ora_res <- compareCluster(
    geneClusters  = gene_lists,
    fun           = "enrichGO",
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  if (is.null(ora_res) || nrow(as.data.frame(ora_res)) == 0) stop("No significant enrichment found.")
  
  df_plot <- as.data.frame(ora_res) %>%
    dplyr::mutate(GeneRatio_num = sapply(strsplit(as.character(GeneRatio), "/"),
                                         function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::slice_min(order_by = p.adjust, n = top_n, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    tidyr::separate(Cluster, into = c("Subtype", "Direction"), sep = "\\|\\|", remove = FALSE)
  
  # facet labels without the "ALS Exc " / "ALS-FTD Exc " prefix
  short_labels <- sub("^ALS(-FTD)? Exc ", "", names(deg_files))
  
  df_plot$Subtype   <- factor(df_plot$Subtype, levels = names(deg_files), labels = short_labels)
  df_plot$Direction <- factor(df_plot$Direction, levels = c("Upregulated", "Downregulated"))
  
  x_labels <- c(
    "Upregulated"   = expression("Increased in TDP-43"^"Low"),
    "Downregulated" = expression("Decreased in TDP-43"^"Low")
  )
  
  ggplot(df_plot, aes(x = Direction, y = Description)) +
    geom_point(aes(size = GeneRatio_num, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
    scale_x_discrete(labels = x_labels) +
    scale_size_continuous(name = "GeneRatio") +
    facet_grid(~ Subtype, switch = "x", scales = "free_x") +
    theme_bw() +
    theme_thesis +
    theme(
      axis.line        = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      axis.title       = element_blank(),
      strip.placement  = "outside",
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = 11, margin = margin(t = 10)),
      panel.border     = element_rect(color = "black", fill = NA, linetype = "dotted"),
      panel.grid.major = element_line(color = "grey90"),
      panel.spacing    = unit(0.5, "lines")
    ) +
    labs(title = title_text)
}

# -------------------------------------------------------------------------
# 3. Figure builder: A-D = volcanoes (own block), E = GO ORA
tag_theme <- theme(plot.tag = element_text(size = 15, face = "bold"))
add_tag   <- function(p, letter) p + labs(tag = letter) + tag_theme

build_figure <- function(volcano_plots, go_plot, fig_title) {
  
  volcano_block <- wrap_plots(
    add_tag(volcano_plots[[1]], "A"),
    add_tag(volcano_plots[[2]], "B"),
    add_tag(volcano_plots[[3]], "C"),
    add_tag(volcano_plots[[4]], "D"),
    ncol = 2
  )
  
  (wrap_elements(full = volcano_block) / add_tag(go_plot, "E")) +
    plot_layout(heights = c(2, 1.15)) +
    plot_annotation(
      title = fig_title,
      theme = theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b = 10))
      )
    )
}

# -------------------------------------------------------------------------
# 4. ALS
volc_als <- list(
  plot_volcano_fig(deg_als_themis, "THEMIS"),
  plot_volcano_fig(deg_als_rorb,   "RORB"),
  plot_volcano_fig(deg_als_linc,   "LINC00507"),
  plot_volcano_fig(deg_als_fezf2,  "FEZF2")
)

als_exc_files <- list(
  "ALS Exc THEMIS"    = "DEG_MAST_ALS_Exc_THEMIS.csv",
  "ALS Exc RORB"      = "DEG_MAST_ALS_Exc_RORB.csv",
  "ALS Exc LINC00507" = "DEG_MAST_ALS_Exc_LINC00507.csv",
  "ALS Exc FEZF2"     = "DEG_MAST_ALS_Exc_FEZF2.csv"
)
go_als <- make_go_plot(als_exc_files, "GO BP ORA: excitatory neurons in ALS")

fig_als <- build_figure(volc_als, go_als, "DEG and GO analysis of excitatory neurons in ALS")

ggsave(file.path(Plots_dir, "Fig13_Figure_Exc_ALS_volcano_GO.svg"),
       fig_als, width = 16, height = 20, device = "svg")

# -------------------------------------------------------------------------
# 5. ALS-FTD
volc_ftd <- list(
  plot_volcano_fig(deg_ftd_themis, "THEMIS"),
  plot_volcano_fig(deg_ftd_rorb,   "RORB"),
  plot_volcano_fig(deg_ftd_linc,   "LINC00507"),
  plot_volcano_fig(deg_ftd_fezf2,  "FEZF2")
)

alsftd_exc_files <- list(
  "ALS-FTD Exc THEMIS"    = "DEG_MAST_ALSFTD_Exc_THEMIS.csv",
  "ALS-FTD Exc RORB"      = "DEG_MAST_ALSFTD_Exc_RORB.csv",
  "ALS-FTD Exc LINC00507" = "DEG_MAST_ALSFTD_Exc_LINC00507.csv",
  "ALS-FTD Exc FEZF2"     = "DEG_MAST_ALSFTD_Exc_FEZF2.csv"
)
go_ftd <- make_go_plot(alsftd_exc_files, "GO BP ORA: excitatory neurons in ALS-FTD")

fig_ftd <- build_figure(volc_ftd, go_ftd, "DEG and GO analysis of excitatory neurons in ALS-FTD")

ggsave(file.path(Plots_dir, "Fig14_Exc_ALSFTD_volcano_GO.svg"),
       fig_ftd, width = 16, height = 20, device = "svg")
