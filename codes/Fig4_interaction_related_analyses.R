# ==============================================================================
# Script Name: Figure4.R
# Description: Ligand-receptor interactions across tumor ecosystems
# Author: Zhang Lab
# ==============================================================================

# ==============================================================================
# Figure 4B: Chord plots illustrating cell-cell interactions across TME subtypes
# ==============================================================================

library(CellChat)

# 1. Path Definitions (Using Relative Paths for GitHub Reproducibility)
# Note: Ensure the files follow consistent case-sensitivity for Linux compatibility
file_list <- c(
  TLS      = "data/processed/rank_TLS_interaction_tier3.csv",
  Lymphoid = "data/processed/rank_Lymphoid_interaction_tier3.csv",
  Myeloid  = "data/processed/rank_Myeloid_interaction_tier3.csv",
  Steroid  = "data/processed/rank_Steroid_interaction_tier3.csv",
  Stromal  = "data/processed/rank_Stromal_interaction_tier3.csv", 
  Desert   = "data/processed/rank_Desert_interaction_tier3.csv"
)

output_pdf <- "results/figures/New_Rank_6TME_network_plot.pdf"
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

# 2. Canonical Cell Type Ordering
celltype_order <- c(
  "B", "CD4", "CD8", "NK", "DC", "Macrophage", 
  "Neutrophil", "Mast", "EC", "Fb", "Tumor"
)

# 3. Defined Color Palette for Cell Lineages
cell_colors <- c(
  "B"          = "#f08940",
  "CD4"        = "#64ba9f",
  "CD8"        = "#c9e3c1",
  "NK"         = "#f7e673",
  "DC"         = "#f6cbde",
  "Macrophage" = "#b37cb3",
  "Neutrophil" = "#d6251f",
  "Mast"       = "#f4af63",
  "EC"         = "#bcb7d6",
  "Fb"         = "#653b90",
  "Tumor"      = "#2572a9"
)

# 4. Multi-panel Visualization
pdf(output_pdf, width = 14, height = 3.5)

# Configure a 1-row, 6-column layout matrix
par(
  mfrow = c(1, 6),
  mar = c(1, 2.5, 4, 1) # Bottom, Left, Top, Right
)

for (nm in names(file_list)) {
  
  if (!file.exists(file_list[nm])) {
    stop(paste("Required input file not found for subtype:", nm, "-", file_list[nm]))
  }
  
  # Load interaction strength matrix
  interaction_raw <- read.csv(file_list[nm], row.names = 1, check.names = FALSE)
  interaction_matrix <- as.matrix(interaction_raw)
  
  # Robust filtering: Keep only cell types present in both the matrix and the canonical list
  valid_cells <- celltype_order[celltype_order %in% rownames(interaction_matrix) & 
                                  celltype_order %in% colnames(interaction_matrix)]
  
  if (length(valid_cells) == 0) {
    stop(paste("No matching canonical cell types found in dataset for:", nm))
  }
  
  # Align and subset matrix to the specified order safely
  sub_matrix <- interaction_matrix[valid_cells, valid_cells, drop = FALSE]
  color_use  <- cell_colors[rownames(sub_matrix)]
  
  # Render Circle Plot via CellChat
  netVisual_circle(
    sub_matrix,
    weight.scale      = TRUE,
    edge.width.max    = 4,
    arrow.size        = 0.35,
    label.edge        = FALSE,
    vertex.label.cex  = 1.3,
    color.use         = color_use,
    title.name        = paste0("TME-", nm)
  )
}

dev.off()







# ==============================================================================
# Figure 4F, S9F and S9G. Gini importance Bar plots
# ==============================================================================


library(ggplot2)
library(dplyr)
library(tidytext)

#Dataset Initialization and Logic Transformation
input_file  <- "data/processed/TME_importance_analysis_top20.csv"
output_pdf  <- "results/figures/TME_Subcluster_Contribution_Barplot.pdf"

if (!file.exists(input_file)) {
  stop(paste("Input data file missing:", input_file))
}

plot_data <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(Subtype = factor(Subtype, levels = c(
    "TLS", "Lymphoid", "Myeloid", "Steroid", "Stromal", "Desert"
  ))) %>%
  mutate(Direction = if_else(Contribution_Score > 0, "Up", "Down")) %>%
  mutate(Stars = case_when(
    Pvalue < 0.001 ~ "***",
    Pvalue < 0.01  ~ "**",
    Pvalue < 0.05  ~ "*",
    TRUE           ~ ""
  )) %>%
  mutate(CellType = reorder_within(CellType, Contribution_Score, Subtype))

#Multi-panel Visualization Configuration
p <- ggplot(plot_data, aes(x = Contribution_Score, y = CellType, fill = Direction)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.1) +
  
  # Fixed Alignment Bug: Embedded hjust expression within aes() to prevent facet re-indexing mismatch
  geom_text(
    aes(
      label = Stars, 
      hjust = if_else(Contribution_Score >= 0, -0.2, 1.2)
    ), 
    size = 3.5, 
    family = "sans"
  ) +
  
  facet_wrap(~Subtype, scales = "free_y", ncol = 6) +
  scale_y_reordered() + 
  scale_fill_manual(values = c("Up" = "#da6763", "Down" = "#4991c5")) +
  
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),     
    panel.grid.minor = element_blank(),     
    panel.border     = element_rect(linewidth = 0.8, color = "black", fill = NA),
    axis.text        = element_text(color = "black", size = 9),
    axis.title       = element_text(size = 11),
    strip.text       = element_text(size = 13, face = "bold"),
    strip.background = element_blank(),     
    legend.position  = "none"        
  ) +
  
  # Reserve symmetric margins for flanking significance labels
  scale_x_continuous(expand = expansion(mult = c(0.25, 0.25))) +
  labs(x = "Contribution Score", y = NULL)

#Save Figure with Publication-grade Parameters
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(output_pdf, p, width = 18, height = 3.5, device = "pdf")


# ==============================================================================
# Figure S9B-E and H. Half-circle dot plots
# ==============================================================================


















