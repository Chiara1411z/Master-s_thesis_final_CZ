library(qs2)
library(tidyverse)
library(Seurat)
library(patchwork)
library(scales)

# Input and Output
# Input are single and merged, filtered seurat objects
output_base        <- "/.../"
SeuratObjects_dir  <- file.path(output_base, "SeuratObjects")
Annotations_dir    <- file.path(output_base, "Annotations")
Samples_path       <- file.path(output_base, "Samples", "Samples.qs2")
Donors_path        <- file.path(output_base, "Samples", "Singularity_Donors.qs2")
Plots_dir          <- file.path(output_base, "Plots")

Samples <- qs_read(Samples_path)
Donors  <- qs_read(Donors_path)

# Shared x-axis group labels
group_labels <- c(
  "NeuN_positive\nTDP_positive" = expression(atop("NeuN+", "TDP-43"^"High")),
  "NeuN_positive\nTDP_negative" = expression(atop("NeuN+", "TDP-43"^"Low")),
  "NeuN_negative\nTDP_positive" = expression(atop("NeuN-", "TDP-43"^"High")),
  "NeuN_negative\nTDP_negative" = expression(atop("NeuN-", "TDP-43"^"Low"))
)


# Shared text styles & bold axes theme
panel_title      <- element_text(hjust = 0.5, size = 15, face = "bold")
main_title_theme <- theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
tag_theme        <- theme(plot.tag   = element_text(size = 18, face = "bold"))

# Helper function for QC filtering
filter_by_MADs <- function(Var, nMads) {
  Var <- as.vector(Var)
  MAD <- mad(Var, center = median(Var, na.rm = TRUE), constant = 1, na.rm = TRUE)
  as.vector(Var < median(Var, na.rm = TRUE) - nMads * MAD | Var > median(Var, na.rm = TRUE) + nMads * MAD)
}

# Helper function for pie charts
plot_pie <- function(df, title) {
  df <- df %>%
    mutate(
      Percentage = round(Count / sum(Count) * 100),
      Label = paste0(format(Count, big.mark = ","), "\n(", Percentage, "%)")
    )
  
  ggplot(df, aes(x = "", y = Count, fill = Status)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 3.5) +
    scale_fill_manual(values = c("Kept" = "#7FFFD4", "Removed" = "#FF1493")) +
    theme_void() +
    labs(title = title, fill = "Status") +
    theme(
      plot.title = panel_title,
      legend.position = "right"
    )
}

# Helper function for funnel plots
make_funnel_plot <- function(df, stage_cols, stage_labels, title, fill_color = "#BF3EFF") {
  df %>%
    dplyr::select(dplyr::all_of(stage_cols)) %>%
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Stage", values_to = "N") %>%
    dplyr::group_by(Stage) %>%
    dplyr::summarise(Total_N = sum(N, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Stage = factor(stage_labels[Stage], levels = stage_labels)) %>%
    ggplot(aes(x = Stage, y = Total_N)) +
    geom_col(fill = fill_color) +
    geom_text(aes(label = Total_N), vjust = -0.3, size = 2.8) +
    theme_classic() +
    theme(
      axis.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
      plot.title = panel_title
    ) +
    labs(x = "Stage", y = "Nuclei [N]", title = title)
}


# -----------------------------------------------------------------------
# FIGURE 8 BLOCK
# 8A: Nuclei recovery per sample
fig8a_data <- Samples %>%
  mutate(ID_label = strip_gex(ID))

fig8a <- ggplot(fig8a_data, aes(x = ID_label, y = nCells_Filtered)) +
  geom_col(fill = "#A020F0") +
  geom_text(aes(label = nCells_Filtered), vjust = -0.3, size = 2.5) +
  theme_classic() +
  theme(
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
    plot.title = panel_title
  ) +
  labs(x = "Sample", y = "Nuclei [N]", title = "Nuclei recovery for each GEX sample")

# 8B: Nuclei before/after doublet removal per condition
fig8b_data <- Samples %>%
  group_by(tdp_status, neun_status) %>%
  summarise(Before = sum(nCells_QC_filtered, na.rm = TRUE),
            After  = sum(nCells_Filtered, na.rm = TRUE), .groups = "drop") %>%
  mutate(Group = factor(paste0(neun_status, "\n", tdp_status), levels = names(group_labels))) %>%
  pivot_longer(cols = c(Before, After), names_to = "Stage", values_to = "n") %>%
  mutate(Stage = factor(Stage, levels = c("Before", "After")))

fig8b <- ggplot(fig8b_data, aes(x = Group, y = n, fill = Stage)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Before" = "#FF1493", "After" = "#7FFFD4")) +
  scale_x_discrete(labels = group_labels) +
  theme_classic() + 
  theme(
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    plot.title = panel_title
  ) +
  labs(y = "Nuclei [N]", x = "Condition", title = "Nuclei distribution before and after doublet removal per condition")

# 8C: Doublet detection rate per sample
Doublet_Rates_path <- file.path(output_base, "Samples", "Doublet_Rates.qs2")
if (file.exists(Doublet_Rates_path)) {
  Doublet_Rates <- qs_read(Doublet_Rates_path)
  
  sample_levels <- if (is.factor(Doublet_Rates$Sample)) {
    unique(strip_gex(levels(Doublet_Rates$Sample)))
  } else {
    unique(strip_gex(Doublet_Rates$Sample))
  }
  Doublet_Rates <- Doublet_Rates %>%
    mutate(Sample = factor(strip_gex(Sample), levels = sample_levels))
  
  fig8c <- ggplot(Doublet_Rates, aes(x = n_ArtDoubl, y = Pct_Doublets * 100, color = Sample, fill = Sample)) +
    geom_vline(xintercept = 5000, col = "red", alpha = 0.8) +
    geom_line() +
    geom_point(size = 2, shape = 21, color = "black") +
    scale_x_continuous(trans = "log10") +
    theme_classic() +
    theme(
      axis.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = panel_title
    ) +
    labs(x = "Artificial doublets [N]", y = "Doublets [%]", title = "Doublet detection rate per GEX sample")
} else {
  fig8c <- ggplot() + theme_void() + labs(title = "(data not available)")
}

# 8D: Funnel plot
FANS_MergedQC <- qs_read(file.path(SeuratObjects_dir, "FANS_Merged_QCFiltered.qs2"))
merged_qc_counts_p1 <- as.data.frame(FANS_MergedQC@meta.data) %>%
  dplyr::count(ID, name = "nCells_MergedQC")

Samples_funnel <- Samples %>%
  left_join(merged_qc_counts_p1, by = "ID") %>%
  mutate(nCells_Final = nCells_MergedQC)

p1_stage_cols <- c("nCells_CR", "nCells_QC_filtered", "nCells_Filtered", "nCells_Final")
p1_stage_labels <- c(
  nCells_CR          = "Total",
  nCells_QC_filtered = "After initial QC",
  nCells_Filtered    = "After doublet removal",
  nCells_Final       = "Final"
)

fig8d <- make_funnel_plot(Samples_funnel, p1_stage_cols, p1_stage_labels,
                          "Retained nuclei at each processing stage", "#BF3EFF")

# 8E: Overall nuclei survival
p1_total   <- sum(Samples_funnel$nCells_CR, na.rm = TRUE)
p1_final   <- sum(Samples_funnel$nCells_Final, na.rm = TRUE)
fig8e_data <- data.frame(
  Status = c("Kept", "Removed"), 
  Count = c(p1_final, p1_total - p1_final)
)
fig8e <- plot_pie(fig8e_data, "Overall nuclei survival")

# 8F: Reasons for removed nuclei
FANS_PreA10 <- qs_read(file.path(SeuratObjects_dir, "FANS_Merged.qs2"))
fig8f_data <- FANS_PreA10@meta.data %>%
  transmute(
    features = filter_by_MADs(Features_by_Counts, 3),
    nCount   = filter_by_MADs(nCount_RNA_log1p, 3),
    mt       = as.vector(percent.mt > 1),
    top20    = filter_by_MADs(Pct_Counts_In_Top_20_Genes, 5),
    ribo     = as.vector(percent.ribo > 1)
  ) %>%
  pivot_longer(everything(), names_to = "Reason", values_to = "Flagged") %>%
  mutate(Flagged = tidyr::replace_na(Flagged, FALSE)) %>%
  filter(Flagged) %>%
  dplyr::count(Reason)

fig8f <- ggplot(fig8f_data, aes(x = fct_reorder(Reason, n, .desc = TRUE), y = n)) +
  geom_col(fill = "#BF3EFF") +
  geom_text(aes(label = n), vjust = -0.3, size = 2.5) +
  theme_classic() +
  theme(
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    plot.title = panel_title
  ) +
  labs(x = "Reasons", y = "Removed nuclei [N]", title = "Reasons for removed nuclei")

# Combined Figure 8
fig8 <- (fig8a | fig8b) /
  (fig8c | fig8d) /
  (fig8e | fig8f) +
  plot_annotation(
    title = "QC and doublet removal per sample", 
    tag_levels = "A",
    theme = main_title_theme
  ) &
  tag_theme

ggsave(file.path(Plots_dir, "Fig8_QC_Doublet_Path1_new.svg"), fig8,
       device = "svg", width = 16, height = 16, dpi = 300)


# -----------------------------------------------------------------------
# FIGURE 9 BLOCK
# 9A: Genotype distribution per GEX sample
fig9a_data <- Donors %>%
  mutate(
    Library = strip_gex(sub("_donor.*", "", ID)),
    Donor_raw = sub(".*_donor", "", ID),
    Donor = as.character(as.numeric(Donor_raw) + 1)
  )

fig9a <- ggplot(fig9a_data, aes(x = Library, y = nCells_Filtered, fill = Donor)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("1" = "#8DEEEE", "2" = "#00BFFF", "3" = "#1C86EE", "4" = "#104E8B")) +
  theme_classic() +
  theme(
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
    plot.title = panel_title
  ) +
  labs(x = "Sample", y = "Nuclei [N]", title = "Distribution of genotypes per GEX sample", fill = "Donor")

# 9B: Retained nuclei after souporcell and doublet removal
fig9b_data <- bind_rows(
  Samples %>%
    group_by(tdp_status, neun_status) %>%
    summarise(N = sum(nCells_QC_filtered, na.rm = TRUE), .groups = "drop") %>%
    mutate(Stage = "Before"),
  Donors %>%
    group_by(tdp_status, neun_status) %>%
    summarise(N = sum(nCells_QC_filtered, na.rm = TRUE), .groups = "drop") %>%
    mutate(Stage = "After Souporcell"),
  Donors %>%
    group_by(tdp_status, neun_status) %>%
    summarise(N = sum(nCells_Filtered, na.rm = TRUE), .groups = "drop") %>%
    mutate(Stage = "After scDblFinder")
) %>%
  mutate(
    Group = factor(paste0(neun_status, "\n", tdp_status), levels = names(group_labels)),
    Stage = factor(Stage, levels = c("Before", "After Souporcell", "After scDblFinder"))
  )

fig9b <- ggplot(fig9b_data, aes(x = Group, y = N, fill = Stage)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Before" = "#FF1493", "After Souporcell" = "#98F5FF", "After scDblFinder" = "#7FFFD4")) +
  scale_x_discrete(labels = group_labels) +
  theme_classic() + 
  theme(
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"), # Waagerecht & untereinander gestapelt
    plot.title = panel_title
  ) +
  labs(y = "Nuclei [N]", x = "Condition", title = "Retained nuclei after souporcell and doublet removal")

# 9C: Reasons for removed nuclei
FANS_Donor_PreQC_path <- file.path(SeuratObjects_dir, "FANS_Donor_Merged.qs2")
if (file.exists(FANS_Donor_PreQC_path)) {
  FANS_Donor_PreQC <- qs_read(FANS_Donor_PreQC_path)
  
  fig9c_data <- FANS_Donor_PreQC@meta.data %>%
    transmute(
      features = filter_by_MADs(Features_by_Counts, 3),
      nCount   = filter_by_MADs(nCount_RNA_log1p, 3),
      mt       = as.vector(percent.mt > 1),
      top20    = filter_by_MADs(Pct_Counts_In_Top_20_Genes, 5),
      ribo     = as.vector(percent.ribo > 1)
    ) %>%
    pivot_longer(everything(), names_to = "Reason", values_to = "Flagged") %>%
    mutate(Flagged = tidyr::replace_na(Flagged, FALSE)) %>%
    filter(Flagged) %>%
    dplyr::count(Reason)
  
  fig9c <- ggplot(fig9c_data, aes(x = fct_reorder(Reason, n, .desc = TRUE), y = n)) +
    geom_col(fill = "#97FFFF") +
    geom_text(aes(label = n), vjust = -0.3, size = 2.5) +
    theme_classic() +
    theme(
      axis.text = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      plot.title = panel_title
    ) +
    labs(x = "Reasons", y = "Nuclei [N]", title = "Reasons for removed nuclei")
  
} else {
  fig9c <- ggplot() + theme_void() + labs(title = "(data not available)")
}

# 9D: Overall nuclei survival
FANS_Donor_MergedQC <- qs_read(file.path(SeuratObjects_dir, "FANS_Donor_Merged_QCFiltered.qs2"))
merged_qc_counts_p2 <- as.data.frame(FANS_Donor_MergedQC@meta.data) %>%
  dplyr::count(ID, name = "nCells_MergedQC")

Donors_funnel <- Donors %>%
  left_join(merged_qc_counts_p2, by = "ID") %>%
  mutate(nCells_Final = nCells_MergedQC)

p2_total   <- sum(Donors_funnel$nCells_Souporcell, na.rm = TRUE)
p2_final   <- sum(Donors_funnel$nCells_Final, na.rm = TRUE)
fig9d_data <- data.frame(
  Status = c("Kept", "Removed"), 
  Count = c(p2_final, p2_total - p2_final)
)

fig9d <- plot_pie(fig9d_data, "Overall nuclei survival")

# Combined Figure 9
fig9 <- (fig9a | fig9b) /
  (fig9c | fig9d) +
  plot_annotation(
    title = "Retained nuclei per GEX sample with previous genotyping", 
    tag_levels = "A",
    theme = main_title_theme
  ) &
  tag_theme

ggsave(file.path(Plots_dir, "Fig9_QC_Doublet_Path2_new.svg"), fig9,
       device = "svg", width = 16, height = 10, dpi = 300)


# -----------------------------------------------------------------------
# FIGURE 10 BLOCK
# 10A: RNA contamination per sample
fig10a <- Samples %>%
  arrange(SoupX_Pct_Ambient) %>%
  mutate(ID_label = strip_gex(ID),
         ID_label = factor(ID_label, levels = unique(ID_label))) %>%
  ggplot(aes(x = ID_label, y = SoupX_Pct_Ambient)) +
  geom_col(fill = "#4682B4") +
  geom_text(aes(label = paste0(round(SoupX_Pct_Ambient, 1), "%")), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = expansion(mult = c(0, 0.15))) +
  theme_classic() +
  theme(
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = panel_title
  ) +
  labs(x = "Sample", y = "Estimated contamination [%]",
       title = "RNA contamination per sample")

# 10B: RNA contamination across conditions
fig10b <- Samples %>%
  mutate(
    PhenotypeGroup = case_when(
      neun_status == "NeuN_positive" & tdp_status == "TDP_positive" ~ "NeuN+, TDP-43 high",
      neun_status == "NeuN_positive" & tdp_status == "TDP_negative" ~ "NeuN+, TDP-43 low",
      neun_status == "NeuN_negative" & tdp_status == "TDP_positive" ~ "NeuN-, TDP-43 high",
      neun_status == "NeuN_negative" & tdp_status == "TDP_negative" ~ "NeuN-, TDP-43 low"
    ),
    PhenotypeGroup = factor(PhenotypeGroup, levels = c(
      "NeuN+, TDP-43 high",  # +/+
      "NeuN+, TDP-43 low",   # +/-
      "NeuN-, TDP-43 high",  # -/+
      "NeuN-, TDP-43 low"    # -/-
    ))
  ) %>%
  ggplot(aes(x = PhenotypeGroup, y = SoupX_Pct_Ambient, fill = PhenotypeGroup)) +
  scale_fill_manual(values = c(
    "NeuN+, TDP-43 high" = "#DAA520",
    "NeuN+, TDP-43 low"  = "#800000",
    "NeuN-, TDP-43 high" = "#F0E68C",
    "NeuN-, TDP-43 low"  = "#C00000"
  )) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, height = 0, size = 2, color = "black") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_x_discrete(
    labels = c(
      "NeuN+, TDP-43 high" = expression(atop("NeuN+", "TDP-43"^"High")),
      "NeuN+, TDP-43 low"  = expression(atop("NeuN+", "TDP-43"^"Low")),
      "NeuN-, TDP-43 high" = expression(atop("NeuN-", "TDP-43"^"High")),
      "NeuN-, TDP-43 low"  = expression(atop("NeuN-", "TDP-43"^"Low"))
    )
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    legend.position = "none",
    plot.title = panel_title
  ) +
  labs(x = "Condition", y = "Estimated contamination [%]",
       title = "RNA contamination across conditions")

# Combined Figure 10
fig10 <- (fig10a | fig10b) +
  plot_annotation(
    title = "Ambient RNA contamination assessment",
    tag_levels = "A",
    theme = main_title_theme
  ) &
  tag_theme

ggsave(file.path(Plots_dir, "Fig10_SoupX_Ambient_RNA_new.svg"), fig10,
       device = "svg", width = 14, height = 6, dpi = 300)
