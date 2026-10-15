library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(ggpubr)

# ---------------------------------------------------------
# Input and Output
# Input is qPCR output
inp_path  <- "/.../"
outp_path <- "/.../"

df_genestudy <- read_csv(inp_path)

target_colors <- c(
  "HC"               = "#DAA520",
  "Carrier\nC9orf72" = "#A6D633",
  "fALS\nC9orf72"    = "#FF1493",
  "ALS"              = "#800000"
)

# Formatted Y-axis label with proper superscript
y_label_expr <- expression(bold("Relative Expression ("*2^{-ΔΔCq}*")"))

# ---------------------------------------------------------
# Data Cleaning & Sample Mapping
df_cleaned <- df_genestudy %>%
  mutate(
    Cq = suppressWarnings(as.numeric(`Mean Cq`)),
    gex_num = case_when(
      Sample %in% c("GEX-ALS", "ALS") ~ 5,
      Sample %in% c("GEX-HC", "HC") ~ 6,
      TRUE ~ as.numeric(str_extract(Sample, "\\d+"))
    ),
    tissue = case_when(
      gex_num >= 5 & gex_num <= 19 ~ "Brain M4\nMotor Cortex",
      gex_num >= 50 & gex_num <= 72 ~ "PBMCs",
      TRUE ~ NA_character_
    ),
    condition = case_when(
      gex_num %in% c(...) ~ "ALS",
      gex_num %in% c(...) ~ "HC",
      gex_num %in% c(...) ~ "Carrier\nC9orf72",
      gex_num %in% c(...) ~ "fALS\nC9orf72",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!gex_num %in% c(66, 67, 68)) %>%
  filter(!is.na(condition), !is.na(tissue), !is.na(Cq)) %>%
  group_by(gex_num, Target, tissue, condition) %>%
  summarise(Cq = mean(Cq, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    condition = factor(condition, levels = c("HC", "Carrier\nC9orf72", "fALS\nC9orf72", "ALS", "ALS\nFTD"))
  )

# ---------------------------------------------------------
# Delta-Delta Cq Normalization (CLTC & HPRT1)
df_ref <- df_cleaned %>%
  filter(toupper(Target) %in% c("CLTC", "HPRT1")) %>%
  group_by(gex_num, tissue) %>%
  summarise(Cq_ref = mean(Cq, na.rm = TRUE), .groups = "drop")

df_delta_cq <- df_cleaned %>%
  filter(!toupper(Target) %in% c("CLTC", "HPRT1")) %>%
  left_join(df_ref, by = c("gex_num", "tissue")) %>%
  mutate(dCq = Cq - Cq_ref)

df_hc_baseline <- df_delta_cq %>%
  filter(condition == "HC") %>%
  group_by(tissue, Target) %>%
  summarise(mean_dCq_HC = mean(dCq, na.rm = TRUE), .groups = "drop")

df_normalized <- df_delta_cq %>%
  left_join(df_hc_baseline, by = c("tissue", "Target")) %>%
  mutate(
    ddCq = dCq - mean_dCq_HC,
    Relative_Expression = 2^(-ddCq)
  )

# Ratio Calculations
df_ratio <- df_normalized %>%
  filter(Target %in% c("hTANK", "hTANK_CE")) %>%
  dplyr::select(gex_num, tissue, condition, Target, Relative_Expression) %>%
  pivot_wider(names_from = Target, values_from = Relative_Expression) %>%
  mutate(hTANK_Ratio = hTANK_CE / hTANK) %>%
  filter(!is.infinite(hTANK_Ratio), !is.na(hTANK_Ratio))

# ---------------------------------------------------------
# Plot
create_qpcr_plot <- function(data, y_var, title, tissue_label, tag, y_label = y_label_expr, y_lim = NULL) {
  
  data <- droplevels(data)   # drop unused groups so only the present conditions are shown
  is_mc <- any(grepl("Motor Cortex", data$tissue))
  
  # Position for the significance labels
  y_max <- max(data[[y_var]], na.rm = TRUE)
  label_pos_y <- y_max * 1.15
  
  p <- ggplot(data, aes(x = condition, y = .data[[y_var]], fill = condition)) +
    geom_boxplot(outlier.shape = NA, color = "black", linewidth = 0.5, width = 0.45) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.5) +
    scale_fill_manual(values = target_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    theme_classic() +
    labs(title = title, subtitle = tissue_label, x = NULL, y = y_label, tag = tag) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 12), 
      plot.tag = element_text(face = "bold", size = 15),
      axis.line = element_line(linewidth = 1),
      axis.text = element_text(color = "black", face = "bold", size = 9),
      axis.title.y = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 0, hjust = 1),
      aspect.ratio = 0.65,
      plot.margin = margin(t = 10, r = 25, b = 10, l = 25)
    )
  
  if (!is.null(y_lim)) {
    p <- p + coord_cartesian(ylim = y_lim)
  }
  
  # Significance statistics
  if (is_mc) {
    p <- p + stat_compare_means(
      method = "t.test",
      label = "p.signif",
      hide.ns = TRUE,
      label.x = 2,
      label.y = label_pos_y,
      size = 4.5
    )
  } else {
    p <- p + stat_compare_means(
      method = "t.test",
      ref.group = "HC",
      p.adjust.method = "BH",
      label = "p.signif",
      hide.ns = TRUE,
      label.y = label_pos_y,
      size = 4.5
    )
  }
  
  return(p)
}

# ---------------------------------------------------------
# Generate subplots
df_hTANK_M4  <- df_normalized %>% filter(Target == "hTANK", grepl("Motor Cortex", tissue))
df_hTANK_PB  <- df_normalized %>% filter(Target == "hTANK", grepl("PBMCs", tissue))

df_CE_M4     <- df_normalized %>% filter(Target == "hTANK_CE", grepl("Motor Cortex", tissue))
df_CE_PB     <- df_normalized %>% filter(Target == "hTANK_CE", grepl("PBMCs", tissue))

df_ratio_M4  <- df_ratio %>% filter(grepl("Motor Cortex", tissue))
df_ratio_PB  <- df_ratio %>% filter(grepl("PBMCs", tissue))

t_total <- "Total hTANK1 expression"
t_ce    <- "hTANK_CE expression"
t_ratio <- "hTANK_CE / total hTANK1 ratio"

p_hTANK_M4 <- create_qpcr_plot(df_hTANK_M4, "Relative_Expression", t_total, "MC",    tag = "A")
p_CE_M4    <- create_qpcr_plot(df_CE_M4,    "Relative_Expression", t_ce,    "MC",    tag = "B")
p_ratio_M4 <- create_qpcr_plot(df_ratio_M4, "hTANK_Ratio",         t_ratio, "MC",    tag = "C", y_label = "Ratio")

p_hTANK_PB <- create_qpcr_plot(df_hTANK_PB, "Relative_Expression", t_total, "PBMCs", tag = "D")
p_CE_PB    <- create_qpcr_plot(df_CE_PB,    "Relative_Expression", t_ce,    "PBMCs", tag = "E")
p_ratio_PB <- create_qpcr_plot(df_ratio_PB, "hTANK_Ratio",         t_ratio, "PBMCs", tag = "F", y_label = "Ratio")

# ---------------------------------------------------------
# Assemble figure
plots_list <- list(
  p_hTANK_M4, p_hTANK_PB,
  p_CE_M4,    p_CE_PB,
  p_ratio_M4, p_ratio_PB
)

design_string <- "
AB
CD
EF
"

final_master_figure <- wrap_plots(plots_list, design = design_string) +
  plot_annotation(
    title = "Expression of hTANK1 and hTANK_CE primer in MC and PBMCs",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(t = 10, b = 15))
    )
  )

# Save
ggsave(
  filename = "GeneStudy_hTANK_hTANK_CE.svg",
  plot = final_master_figure,
  path = outp_path,
  width = 28,
  height = 26,
  units = "cm"
)
