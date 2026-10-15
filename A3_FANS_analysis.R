library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(svglite)

# ---------------------------------------------------------
# Input and Output Setup
# Input: FANS raw numbers
inp_path  <- "/.../"
outp_path <- "/.../"

SHOW_F <- FALSE

df_nuclei_all <- read_xlsx(inp_path) %>%
  filter(if_any(matches("date", ignore.case = TRUE), ~ !str_detect(as.character(.x), "2026-07-22|22/07/2026") | is.na(.x))) %>%
  mutate(Condition = ifelse(Condition == "ALSuALS-FTD", "ALS and ALS-FTD", Condition))

# Data Cleaning & Setup
condition_labels <- c(
  "HC"                = "HC\n(n = 5)",
  "ALS"               = "ALS\n(n = 14)",
  "ALS-FTD"           = "ALS-FTD\n(n = 4)",
  "ALS and ALS-FTD" = "ALS &\nALS-FTD\n(n = 2)"
)

# Apply labels to dataset and lock plotting order
df_nuclei <- df_nuclei_all %>%
  filter(!is.na(Condition)) %>%
  mutate(X_Label = condition_labels[Condition]) %>%
  mutate(
    Condition = factor(Condition, levels = names(condition_labels)),
    X_Label   = factor(X_Label, levels = condition_labels)
  )

# Color Palette
target_color_box <- c(
  "HC" = "#FFCC00",
  "ALS" = "#800000",
  "ALS-FTD" = "#FF1493",
  "ALS and ALS-FTD" = "#CD5555"
)

# ---------------------------------------------------------
# Shared theme
my_theme <- theme_classic() +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 8)),
    axis.title.x      = element_text(face = "bold", size = 10),
    axis.title.y      = element_text(face = "bold", size = 10),
    axis.text         = element_text(color = "black", face = "bold", size = 9),
    axis.text.x       = element_text(angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
    axis.line         = element_line(linewidth = 1),
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 9),
    plot.margin       = margin(t = 10, r = 10, b = 25, l = 10)
  )

# ---------------------------------------------------------
# Helper functions

# Global one-way ANOVA label
get_anova_label <- function(data, variable, show_F = FALSE) {
  f   <- as.formula(paste(variable, "~ Condition"))
  tab <- summary(aov(f, data = data))[[1]]
  p_val <- tab[["Pr(>F)"]][1]
  p_txt <- if (p_val < 0.001) "p < 0.001" else paste0("p = ", format(round(p_val, 2), nsmall = 2))
  if (show_F) {
    paste0("ANOVA, F(", tab$Df[1], ",", tab$Df[2], ") = ",
           format(round(tab[["F value"]][1], 2), nsmall = 2), ", ", p_txt)
  } else {
    paste0("ANOVA, ", p_txt)
  }
}

# Convert adjusted p-values to significance stars for the plot labels
format_p <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# Dunnett's test: every non-HC condition vs HC
run_dunnett <- function(data, variable) {
  fit <- aov(as.formula(paste(variable, "~ Condition")), data = data)
  multcomp::glht(fit, linfct = multcomp::mcp(Condition = "Dunnett"))
}

dunnett_labels <- function(data, variable, y_pos) {
  dun <- run_dunnett(data, variable)
  p   <- summary(dun)$test$pvalues
  tibble::tibble(x = 2:4, y = y_pos, label = format_p(p))
}

make_boxplot <- function(data, variable, title, show_F = SHOW_F) {
  max_y   <- max(data[[variable]], na.rm = TRUE)
  ann     <- dunnett_labels(data, variable, max_y * 1.15)
  aov_lab <- get_anova_label(data, variable, show_F)
  
  ggplot(data, aes(x = X_Label, y = .data[[variable]], fill = Condition)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5, color = "black") +
    geom_jitter(aes(color = Condition), width = 0.15, size = 2.5, alpha = 0.9) +
    scale_fill_manual(values = target_color_box) +
    scale_color_manual(values = target_color_box) +
    geom_text(data = ann, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 5) +
    annotate("text", x = -Inf, y = Inf, label = aov_lab,
             hjust = -0.08, vjust = 1.6, size = 3.6) +
    coord_cartesian(ylim = c(0, max_y * 1.45)) +
    labs(x = "Condition", y = "Events [%]", title = title) +
    my_theme + theme(legend.position = "none")
}

# ---------------------------------------------------------
# Plot A: Total DAPI Events Distributed Across Conditions
df_dapi_condition <- df_nuclei %>%
  group_by(X_Label, Condition) %>%
  summarise(Total_DAPI = sum(`Total DAPI Amount`, na.rm = TRUE), .groups = 'drop')

p_A <- ggplot(df_dapi_condition, aes(x = X_Label, y = Total_DAPI, fill = Condition)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  scale_fill_manual(name = "Condition", values = target_color_box) +
  scale_y_continuous(labels = scales::scientific, expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Condition", y = "Events [N]", title = "Total DAPI events") +
  my_theme

# ---------------------------------------------------------
# Data Prep for Boxplots (Plots B, C, D, E)
df_sample_level <- df_nuclei %>%
  mutate(
    NeuNpos_pct          = (`# NeuN+,TPD43+` + `# NeuN+,TPD43-`) / `Total DAPI Amount` * 100,
    NeuNneg_pct          = (`# NeuN-,TPD43+` + `# NeuN-,TPD43-`) / `Total DAPI Amount` * 100,
    TDP43low_NeuNpos_pct = `# NeuN+,TPD43-` / (`# NeuN+,TPD43+` + `# NeuN+,TPD43-`) * 100,
    TDP43low_NeuNneg_pct = `# NeuN-,TPD43-` / (`# NeuN-,TPD43+` + `# NeuN-,TPD43-`) * 100
  ) %>%
  filter(!(Condition == "ALS" & str_detect(`Sample Name`, regex("pool", ignore_case = TRUE)))) %>%
  group_by(`Sample Name`, Condition, X_Label) %>%
  summarise(
    NeuNpos_pct          = mean(NeuNpos_pct, na.rm = TRUE),
    NeuNneg_pct          = mean(NeuNneg_pct, na.rm = TRUE),
    TDP43low_NeuNpos_pct = mean(TDP43low_NeuNpos_pct, na.rm = TRUE),
    TDP43low_NeuNneg_pct = mean(TDP43low_NeuNneg_pct, na.rm = TRUE),
    .groups = 'drop'
  )


vars <- c(
  "NeuN+ event proportion"                 = "NeuNpos_pct",
  "NeuN- event proportion"                 = "NeuNneg_pct",
  "TDP-43 low proportion of NeuN+"         = "TDP43low_NeuNpos_pct",
  "TDP-43 low proportion of NeuN-"         = "TDP43low_NeuNneg_pct"
)

cat("\n=======================================================\n")
cat("  ONE-WAY ANOVA + DUNNETT'S TEST (COMPARISONS VS HC)   \n")
cat("=======================================================\n")
for (nm in names(vars)) {
  cat("\n---", nm, "---\n")
  fit <- aov(as.formula(paste(vars[[nm]], "~ Condition")), data = df_sample_level)
  print(summary(fit))
  print(summary(run_dunnett(df_sample_level, vars[[nm]])))
}

# ---------------------------------------------------------
# ALL NUMBERS FOR THE FIGURE (console + CSV)
panel_vars <- c(B = "NeuNpos_pct", C = "NeuNneg_pct",
                D = "TDP43low_NeuNpos_pct", E = "TDP43low_NeuNneg_pct")

stat_table <- bind_rows(lapply(names(panel_vars), function(pn) {
  v   <- panel_vars[[pn]]
  fit <- aov(as.formula(paste(v, "~ Condition")), data = df_sample_level)
  a   <- summary(fit)[[1]]
  dun <- summary(run_dunnett(df_sample_level, v))$test
  bind_rows(
    tibble::tibble(Panel = pn, Test = "One-way ANOVA", Comparison = "all groups",
                   Estimate = NA_real_, SE = NA_real_, Statistic = a[["F value"]][1],
                   df = paste0(a$Df[1], ",", a$Df[2]), p = a[["Pr(>F)"]][1]),
    tibble::tibble(Panel = pn, Test = "Dunnett vs HC", Comparison = names(dun$coefficients),
                   Estimate = unname(dun$coefficients), SE = unname(dun$sigma),
                   Statistic = unname(dun$tstat), df = as.character(a$Df[2]),
                   p = unname(dun$pvalues))
  )
}))

desc_table <- bind_rows(lapply(names(panel_vars), function(pn) {
  v <- panel_vars[[pn]]
  df_sample_level %>%
    group_by(Condition) %>%
    summarise(n      = n(),
              Mean   = mean(.data[[v]]),
              SD     = sd(.data[[v]]),
              Median = median(.data[[v]]),
              Q1     = quantile(.data[[v]], 0.25),
              Q3     = quantile(.data[[v]], 0.75),
              .groups = "drop") %>%
    mutate(Panel = pn, .before = 1)
}))

panelA_table <- df_dapi_condition %>%
  mutate(n_acquisitions = as.integer(table(df_nuclei$Condition)[as.character(Condition)])) %>%
  dplyr::select(Condition, Total_DAPI, n_acquisitions)

options(width = 200, pillar.sigfig = 4)
cat("\n== Panel A: total DAPI events ==\n");         print(panelA_table)
cat("Overall total DAPI events:", sum(panelA_table$Total_DAPI), "\n")
cat("\n== Panels B-E: ANOVA and Dunnett ==\n");       print(stat_table, n = Inf)
cat("\n== Panels B-E: descriptive statistics ==\n"); print(desc_table, n = Inf)

write.csv(stat_table,   file.path(outp_path, "Fig6_statistics.csv"),   row.names = FALSE)
write.csv(desc_table,   file.path(outp_path, "Fig6_descriptives.csv"), row.names = FALSE)
write.csv(panelA_table, file.path(outp_path, "Fig6_panelA.csv"),       row.names = FALSE)

# ---------------------------------------------------------
# Generate Boxplots
p_neunpos    <- make_boxplot(df_sample_level, "NeuNpos_pct", "NeuN+ event proportion")
p_neunneg    <- make_boxplot(df_sample_level, "NeuNneg_pct", "NeuN- event proportion")
p_TDPlow_pos <- make_boxplot(df_sample_level, "TDP43low_NeuNpos_pct",
                             expression(bold("TDP-43"^"Low" * " event proportion of NeuN+")))
p_TDPlow_neg <- make_boxplot(df_sample_level, "TDP43low_NeuNneg_pct",
                             expression(bold("TDP-43"^"Low" * " event proportion of NeuN-")))

# ---------------------------------------------------------
# Combine Figures and Save
layout_design <- "
  AL
  BC
  DE
"

LEGEND_SHIFT_LEFT <- -30
LEGEND_SHIFT_UP   <- -30

combine_figure <- p_A + p_neunpos + p_neunneg + p_TDPlow_pos + p_TDPlow_neg + guide_area() +
  plot_layout(design = layout_design, guides = "collect") +
  plot_annotation(
    title = "Nuclei distribution and TDP-43 alterations across conditions",
    tag_levels = "A",
    theme = theme(
      plot.title            = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 15)),
      legend.title          = element_text(face = "bold", size = 10),
      legend.text           = element_text(size = 9),
      legend.background     = element_rect(fill = "white", color = "grey90"),
      legend.justification  = c(0, 1),
      legend.box.margin     = margin(t = LEGEND_SHIFT_UP, l = LEGEND_SHIFT_LEFT)
    )
  ) & 
  theme(plot.tag = element_text(face = "bold", size = 15))

ggsave(
  filename = "Fig6_Nuclei_Distribution_and_TDP43_Alterations.svg",
  plot = combine_figure,
  path = outp_path,
  width = 30,
  height = 34,
  units = "cm"
)

