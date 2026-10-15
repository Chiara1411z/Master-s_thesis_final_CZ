library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork) 

# Input and Output
# Input is CellRanger outpur
inp_path <- "/.../"
outp_path <- "/.../"

df <- read_xlsx(inp_path)
colnames(df) <- c("GEX", "Num_Cells", "Mean_reads", "Median_genes", "Condition")

# Data Wrangling
df_clean <- df %>%
  separate(Condition, into = c("NeuN_status", "TDP43_status"), sep = "/") %>%
  mutate(
    NeuN_status = case_when(
      NeuN_status == "+" ~ "NeuN+",
      NeuN_status == "-" ~ "NeuN-",
      TRUE ~ paste0("NeuN ", NeuN_status)
    ),
    TDP43_status = case_when(
      TDP43_status == "+" ~ "TDP-43 High",
      TDP43_status == "-" ~ "TDP-43 Low",
      TRUE ~ paste0("TDP-43 ", TDP43_status)
    )
  )

# Define target colors and superscript labels
target_colors <- c(
  "TDP-43 Low"  = "#800000",  
  "TDP-43 High" = "#DAA520"
)  

tdp_expr_labels <- c(
  "TDP-43 Low"  = expression("TDP-43"^"Low"), 
  "TDP-43 High" = expression("TDP-43"^"High")
)

# Plot A: Mean reads per cell (Log Scaled)
p1 <- ggplot(df_clean, aes(x = NeuN_status, y = Mean_reads, fill = TDP43_status)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, color = "black") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), 
             size = 2, color = "black", alpha = 0.7) +
  scale_fill_manual(values = target_colors, labels = tdp_expr_labels) +
  scale_y_log10(labels = scales::comma) +
  theme_light(base_size = 14) + 
  labs(
    title = "Mean Reads per Cell",
    x = "Condition",
    y = "Mean reads (log10)",
    fill = "TDP-43 Event"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank() 
  )

# Plot B: Median genes per cell
p2 <- ggplot(df_clean, aes(x = NeuN_status, y = Median_genes, fill = TDP43_status)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, color = "black") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), 
             size = 2, color = "black", alpha = 0.7) +
  scale_fill_manual(values = target_colors, labels = tdp_expr_labels) +
  theme_light(base_size = 14) +
  labs(
    title = "Median Genes per Cell",
    x = "Condition",
    y = "Median genes",
    fill = "TDP-43 Event"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

# Combine plots 
combined_figure <- (p1 + p2) + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A') & 
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 18, face = "bold")
  )

# Save
ggsave(filename = file.path(outp_path, "mean_median_boxplot.svg"), plot = combined_figure, width = 12, height = 6)
