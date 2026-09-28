library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(svglite)

# Clear active workspace memory
rm(list = ls())

# ---------------------------------------------------------
# 1. File Paths & Setup
inp_path  <- "/home/chiara/NAS/Labmembers/Zuehlke_Chiara/Masters_thesis/qPCR/Gene_Studies/2026_04_30_GeneStudy_brainUL2508/genestudy_alldata.xlsx"
outp_path <- "/home/chiara/NAS/Labmembers/Zuehlke_Chiara/Masters_thesis/R/Code/GeneStudy"

# Formatted Y-axis label with proper superscript (2^-ΔΔCq)
y_label_expr <- expression(bold("Relative Expression ("*2^{-ΔΔCq}*")"))

# Read dataset
df_genestudy <- read_xlsx(inp_path)

# ---------------------------------------------------------
# 2. Clean Strings & Metadata Formatting
df_filtered <- df_genestudy %>%
  mutate(Cq = as.numeric(Cq)) %>%
  filter(Content == "Unkn") %>%
  mutate(
    `NeuN Status` = case_when(
      str_detect(Sample, "NeuN\\+") ~ "NeuN Positive",
      str_detect(Sample, "NeuN-")   ~ "NeuN Negative",
      TRUE ~ "Control"
    ),
    `TDP Status` = case_when(
      str_detect(Sample, "TDP43\\+") ~ "TDP-43+",
      str_detect(Sample, "TDP43-")   ~ "TDP-43-",
      TRUE ~ "Control"
    ),
    Concentration = if_else(str_detect(Sample, "10c"), "10)", "1")
  ) %>%
  filter(!Sample %in% c("NeuN-TDP43-", "NeuN-TDP43-10")) %>%
  dplyr::select(-Well, -any_of(c("Set Point", "...10", "Cq Mean", "Cq Std. Dev")))

# ---------------------------------------------------------
# 3. Delta-Delta Cq Normalization (2^-ddCq relative to TDP-43-)
df_averaged <- df_filtered %>%
  group_by(Sample, Target, `NeuN Status`, `TDP Status`, Concentration) %>%
  summarise(Mean_Cq = mean(Cq, na.rm = TRUE), .groups = "drop")

# Housekeeping reference genes (CLTC, HPRT1, huRNU6)
df_housekeeping <- df_averaged %>%
  filter(Target %in% c("CLTC", "HPRT1", "huRNU6")) %>%
  group_by(Sample) %>%
  summarise(HG_Ref = mean(Mean_Cq, na.rm = TRUE), .groups = "drop")

# Calculate Delta Cq
df_normalized <- df_averaged %>%
  left_join(df_housekeeping, by = "Sample") %>%
  mutate(`Delta Cq` = Mean_Cq - HG_Ref)

df_markers <- df_normalized %>% 
  filter(!Target %in% c("CLTC", "HPRT1", "huRNU6"))

# Calculate TDP-43- baseline mean per Target and NeuN Status
df_baseline <- df_markers %>%
  filter(`TDP Status` == "TDP-43-") %>%
  group_by(Target, `NeuN Status`) %>%
  summarise(mean_dCq_control = mean(`Delta Cq`, na.rm = TRUE), .groups = "drop")

# Calculate Delta-Delta Cq and Relative Expression (2^-ddCq)
df_markers <- df_markers %>%
  left_join(df_baseline, by = c("Target", "NeuN Status")) %>%
  mutate(
    ddCq = `Delta Cq` - mean_dCq_control,
    Relative_Expression = 2^(-ddCq)
  )

# Get ordered target list matching original figure
unique_markers <- c("hCX3CR1", "hMBP", "hRBFOX3", "hS100B", "hTANK_CE", "hTARDBP1", "hTARDBP2")
unique_markers <- intersect(unique_markers, unique(df_markers$Target))

# ---------------------------------------------------------
# 4. Generate Marker Boxplots (No Significance Labels)
marker_plots <- lapply(unique_markers, function(marker) {
  
  df_single_marker <- df_markers %>% filter(Target == marker)
  
  p <- ggboxplot(
    df_single_marker, 
    x = "TDP Status", 
    y = "Relative_Expression", 
    add = "jitter", 
    color = "TDP Status", 
    palette = c("TDP-43-" = "#DAA520", "TDP-43+" = "#800000"),         
    title = paste(marker, "Expression"),
    xlab = "TDP-43 Status"
  ) +
    labs(y = y_label_expr) +  # Applied here as standard ggplot2 layer to prevent ggpubr expression bug
    theme_classic() +
    theme(
      legend.position = "none", 
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 10),
      axis.text = element_text(color = "black", face = "bold", size = 9),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
      strip.text = element_text(face = "bold", size = 10)
    ) + 
    facet_wrap(~ `NeuN Status`) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20)))
  
  return(p)
})

names(marker_plots) <- unique_markers

# ---------------------------------------------------------
# 5. Qualitative Binary Cryptic Exon Profiles (hTANK_CE Barplot)
CQ_CUTOFF <- 37.5

df_exon <- df_filtered %>%
  filter(Target == "hTANK_CE") %>%
  filter(Concentration == "1") %>%
  mutate(Cq_clean = if_else(Cq >= CQ_CUTOFF | is.na(Cq), NA_real_, Cq)) %>%
  group_by(Sample, `NeuN Status`, `TDP Status`) %>%
  summarise(
    n_replicates = n(),
    n_detected   = sum(!is.na(Cq_clean)),
    Mean_Cq      = mean(Cq_clean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Exon_Presence = factor(
      if_else(n_detected >= ceiling(n_replicates / 2), "Present", "Absent"),
      levels = c("Absent", "Present")
    )
  )

plot_exon <- ggplot(df_exon, aes(x = `TDP Status`, fill = Exon_Presence)) +
  geom_bar(position = "fill", width = 0.4, color = "black") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_fill_manual(
    values = c("Absent" = "#8B1A1A", "Present" = "#00FF00"),
    drop = FALSE
  ) +
  labs(
    title = "hTANK_CE Cryptic Splicing Detection Rate",
    y = "Percentage of Samples",
    x = "TDP-43 Pathology Status",
    fill = "Cryptic Exon"
  ) +
  theme_pubr() + 
  theme(
    legend.position = "right", 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black", face = "bold", size = 9),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    strip.text = element_text(face = "bold", size = 10)
  ) +
  facet_wrap(~ `NeuN Status`)

# ---------------------------------------------------------
# 6. Combine and Save Figure
all_plots <- c(marker_plots, list(plot_exon))

cat("Total panels assembled:", length(all_plots), "\n")

final_figure <- wrap_plots(all_plots, ncol = 3) + 
  plot_annotation(tag_levels = "A")

setwd(outp_path)
ggsave("qPCR_Master_Thesis_Figure_alldata_final_new.svg", plot = final_figure, width = 14, height = 12)
