library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(scales)
library(patchwork)
library(ggtext)
library(qs2)

# -------------------------------------------------------------------------
# Input and Output
output_base       <- "/home/chiara/NAS/Labmembers/Zuehlke_Chiara/Masters_thesis/R/Code/Cellranger/Seurat_analysis"
SeuratObjects_dir <- file.path(output_base, "SeuratObjects")
Plots_dir         <- file.path(output_base, "Plots")

FANS <- qs_read(file.path(SeuratObjects_dir, "FANS_Integrated_LabelsTransferred.qs2"))

# Set TRUE to add the MASC forest plot as a third row
include_forest <- TRUE

# -------------------------------------------------------------------------
# Shared theme
my_theme <- theme_classic() +
  theme(
    plot.title             = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.title.position    = "panel",  
    axis.title.x           = element_text(face = "bold", size = 10),
    axis.title.y           = element_text(face = "bold", size = 10),
    axis.text              = element_text(color = "black", face = "bold", size = 9),
    axis.text.x            = element_text(face = "bold"),
    axis.text.y            = element_text(face = "bold"),
    axis.line              = element_line(linewidth = 1),
    legend.title           = element_text(face = "bold", size = 10),
    legend.text            = element_text(size = 11),
    legend.justification   = "top"
  )

# Forest Plot laden und Labels direkt im Objekt überschreiben (umgeht Scale-Konflikte)
p_forest <- qs_read(file.path(SeuratObjects_dir, "p_tdp_forest_clean.qs2"))
p_forest$labels$x <- "**Odds Ratio (Log Scale, 95% CI vs Healthy TDP-43<sup>High</sup>)**"
p_forest$labels$y <- "Subcluster"

p_forest <- p_forest +
  labs(
    title = "Selective vulnerability: TDP-43 pathology enrichment per neuronal type", 
    tag = "G"
  ) +
  my_theme +
  theme(
    plot.subtitle = element_text(size = 9, face = "italic", color = "gray30", hjust = 0.5),
    axis.title.x  = ggtext::element_markdown(face = "bold", size = 10, hjust = 0.5), # Markdown-fett erzwungen
    axis.title.y  = element_text(face = "bold", size = 10)
  )

# -------------------------------------------------------------------------
# Clean and extract metadata
meta_df <- FANS@meta.data %>%
  mutate(
    sample_num = as.numeric(str_extract(coalesce(ID, orig.ident), "\\d+")),
    disease = case_when(
      sample_num >= 31 & sample_num <= 47 ~ "ALS",
      sample_num >= 48 & sample_num <= 55 ~ "ALS-FTD",
      TRUE ~ NA_character_
    ),
    neun_status = dplyr::recode(neun_status,
                                "NeuN_negative" = "NeuN-",
                                "NeuN_positive" = "NeuN+")
  ) %>%
  filter(!is.na(disease))

# Global colors for L2.5 cell types
all_cell_types <- sort(unique(meta_df$ID_WNN_L25_Predicted))
global_colors  <- scales::hue_pal()(length(all_cell_types))
names(global_colors) <- all_cell_types

# Subsets (shared by both rows)
is_neuron <- function(x) grepl("^Exc_|^Inh_", x)

df_a <- meta_df %>% filter(disease == "ALS",     neun_status == "NeuN+", is_neuron(ID_WNN_L25_Predicted))
df_b <- meta_df %>% filter(disease == "ALS-FTD", neun_status == "NeuN+", is_neuron(ID_WNN_L25_Predicted))
df_c <- meta_df %>% filter(disease == "ALS-FTD", neun_status == "NeuN-", !is_neuron(ID_WNN_L25_Predicted))

# -------------------------------------------------------------------------
# ROW 1: Cell type composition (mit Tags A, B, C)
plot_comp <- function(df, title, x_label, tag_letter) {
  df_plot <- df %>%
    group_by(ID_WNN_L25_Predicted) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    mutate(proportion = n_cells / sum(n_cells), x_label = x_label)
  
  ggplot(df_plot, aes(x = x_label, y = proportion, fill = ID_WNN_L25_Predicted)) +
    geom_col(position = "fill", color = "black", linewidth = 0.2) +
    geom_text(
      data = df_plot %>% filter(proportion > 0.03),
      aes(label = paste0(ID_WNN_L25_Predicted, "\n(", round(proportion * 100, 1), "%)")),
      position = position_stack(vjust = 0.5),
      size = 2.5
    ) +
    scale_fill_manual(values = global_colors, drop = TRUE, name = "L2.5 cell type") +
    labs(title = title, x = "Subcluster", y = "Proportion", tag = tag_letter) +
    my_theme
}

p1_a <- plot_comp(df_a, "ALS",     "NeuN+", "A")
p1_b <- plot_comp(df_b, "ALS-FTD", "NeuN+", "B")
p1_c <- plot_comp(df_c, "ALS-FTD", "NeuN-", "C")

# -------------------------------------------------------------------------
# ROW 2: TDP-43 status per cell type (mit Tags D, E, F)
prep_data <- function(df) {
  df %>%
    group_by(ID_WNN_L25_Predicted, tdp_status) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(ID_WNN_L25_Predicted) %>%
    mutate(total = sum(n_cells), proportion = n_cells / total)
}

tdp_colors <- c("TDP-43_low" = "#800000", "TDP-43_high" = "#DAA520")
legend_labels_tdp <- expression("TDP-43"^"Low", "TDP-43"^"High")
names(legend_labels_tdp) <- c("TDP-43_low", "TDP-43_high")

create_plot <- function(df, title, tag_letter) {
  ggplot(df, aes(x = ID_WNN_L25_Predicted, y = proportion, fill = tdp_status)) +
    geom_col(position = "fill") +
    geom_text(
      data = df %>% filter(proportion > 0.04),
      aes(label = n_cells),
      position = position_stack(vjust = 0.5),
      size = 3, color = "white"
    ) +
    scale_fill_manual(values = tdp_colors, labels = legend_labels_tdp) +
    coord_flip() +
    labs(title = title, x = "Subcluster", y = "Proportion", fill = "TDP-43 status", tag = tag_letter) +
    my_theme
}

p2_a <- create_plot(prep_data(df_a), "ALS: NeuN+", "D")
p2_b <- create_plot(prep_data(df_b), "ALS-FTD: NeuN+", "E")
p2_c <- create_plot(prep_data(df_c), "ALS-FTD: NeuN-", "F")

# -------------------------------------------------------------------------
# LEGEND HANDLING
get_leg <- function(p) {
  g   <- ggplotGrob(p + theme(legend.position = "right",
                              legend.justification = "top"))
  idx <- which(g$layout$name == "guide-box-right")
  g$grobs[[idx]]
}

leg_l25    <- get_leg(p1_a + labs(fill = "L2.5 cell type (NeuN+)"))
leg_glia   <- get_leg(p1_c + labs(fill = "L2.5 cell type (NeuN-)"))
leg_tdp    <- get_leg(p2_a)
leg_forest <- get_leg(p_forest)

place_legend <- function(g, y = 0.95) {
  wrap_elements(
    full = grid::grobTree(
      g,
      vp = grid::viewport(
        x = 0, y = grid::unit(y, "npc"), just = c("left", "top"),
        width  = grid::grobWidth(g),
        height = grid::grobHeight(g)
      )
    )
  )
}

place_two_legends <- function(g1, g2, y = 0.95, gap = 0.5) {
  gap_u <- grid::unit(gap, "cm")
  wrap_elements(
    full = grid::grobTree(
      grid::grobTree(
        g1,
        vp = grid::viewport(
          x = 0, y = grid::unit(y, "npc"), just = c("left", "top"),
          width = grid::grobWidth(g1), height = grid::grobHeight(g1))
      ),
      grid::grobTree(
        g2,
        vp = grid::viewport(
          x = 0, y = grid::unit(y, "npc") - grid::grobHeight(g1) - gap_u,
          just = c("left", "top"),
          width = grid::grobWidth(g2), height = grid::grobHeight(g2))
      )
    )
  )
}

# -------------------------------------------------------------------------
# ASSEMBLY
no_leg <- theme(legend.position = "none")

plots <- list(
  p1_a + no_leg, p1_b + no_leg, p1_c + no_leg, 
  place_two_legends(leg_l25, leg_glia, y = 0.95)
  p2_a + no_leg, p2_b + no_leg, p2_c + no_leg, 
  place_legend(leg_tdp,    y = 0.93),       
  p_forest + no_leg,                 
  place_legend(leg_forest, y = 0.90)       
)

layout <- "
ABCD
EFGH
IIIJ
"

final_plot <- wrap_plots(plots, design = layout, widths = c(1, 1, 1, 0.5)) +
  plot_annotation(
    title = "Analysis of nuclei on resolution level 2.5",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10))
    )
  ) & 
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(
  filename = file.path(Plots_dir, "Fig12_Combined_Analysis_Nuclei_L25.svg"),
  plot     = final_plot,
  device   = "svg",
  width    = 22,
  height   = 22
)
