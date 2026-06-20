rm(list = ls())
gc()

# Fix parallelization conflicts - MUST be set before loading bigsnpr
Sys.setenv(OMP_NUM_THREADS = 1)

# Install once (comment out later)
# install.packages("bigsnpr")
# install.packages("bigstatsr")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("plotly")
# install.packages("viridis")
# install.packages("RColorBrewer")
# install.packages("scales")
# install.packages("cowplot")
# install.packages("patchwork")
# install.packages("data.table")

library(bigsnpr)
library(bigstatsr)
library(ggplot2)
library(dplyr)
library(plotly)
library(viridis)
library(RColorBrewer)
library(scales)
library(cowplot)
library(patchwork)
library(data.table)

options(bigstatsr.symmetry.check = FALSE)  # speed-up

# Create output directories
if (!dir.exists("PCA_plots_bigsnpr")) dir.create("PCA_plots_bigsnpr")
if (!dir.exists("BigSNP_pngs")) dir.create("BigSNP_pngs")

# Define save_plot function with high-quality settings
save_plot <- function(p, name, width = 7, height = 6) {
  # Save to original directory
  ggsave(
    filename = file.path("PCA_plots_bigsnpr", paste0(name, ".png")),
    plot     = p,
    width    = width,
    height   = height,
    dpi      = 600,
    bg       = "white"
  )
  # Save to BigSNP_pngs directory
  ggsave(
    filename = file.path("BigSNP_pngs", paste0(name, ".png")),
    plot     = p,
    width    = width,
    height   = height,
    dpi      = 600,
    bg       = "white"
  )
  # Also save PDF for publications
  ggsave(
    filename = file.path("BigSNP_pngs", paste0(name, ".pdf")),
    plot     = p,
    width    = width,
    height   = height,
    device   = cairo_pdf
  )
}

# Sanity test
cat("save_plot function exists:", exists("save_plot"), "\n\n")

# -------------------------------------------------------
# DEFINE NATURE-STYLE THEME
# -------------------------------------------------------
theme_nature <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      # Plot elements
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0, margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size, hjust = 0, margin = margin(b = 10)),
      
      # Axes
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      
      # Legend
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.key.size = unit(0.5, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      
      # Panel
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Margins
      plot.margin = margin(10, 10, 10, 10)
    )
}

# -------------------------------------------------------
# DEFINE COLOR PALETTES (Nature-style, colorblind-friendly)
# -------------------------------------------------------
# Custom palette for datasets (2 colors)
colors_dataset <- c(
  "1000G" = "#0072B2",  # Blue
  "Ghana"    = "#D55E00"   # Vermillion/Orange
)

# Custom palette for super populations (5+ colors)
colors_superpop <- c(
  "AFR" = "#E69F00",      # Orange
  "EUR" = "#56B4E9",      # Sky Blue
  "EAS" = "#009E73",      # Bluish Green
  "SAS" = "#F0E442",      # Yellow
  "AMR" = "#CC79A7",      # Reddish Purple
  "Unknown" = "#999999",  # Gray
  "Ghana" = "#A020F0"    # Vermillion for H3
)

# Custom palette for African populations (will use RColorBrewer Set2)
colors_afr_pop <- brewer.pal(8, "Set2")

cat("Nature-style theme and color palettes defined.\n\n")

# -------------------------------------------------------
# 1) Load BED files for bigsnpr (1000G projection)
# -------------------------------------------------------
bed_1000g <- bed("1000G_merged.bed")  # reference (1000G)
bed_h3    <- bed("H3_Africa_order_3220_original_annotation.bed")     # target (Ghana H3)

map_1000g <- bed_1000g$map
map_h3    <- bed_h3$map

fam_1000g <- bed_1000g$fam
fam_h3    <- bed_h3$fam

# Convenient aliases
chr_1000g   <- map_1000g$chromosome
pos_1000g   <- map_1000g$physical.pos
snpid_1000g <- map_1000g$marker.ID
id_1000g    <- fam_1000g$sample.ID

chr_h3   <- map_h3$chromosome
pos_h3   <- map_h3$physical.pos
snpid_h3 <- map_h3$marker.ID
id_h3    <- fam_h3$sample.ID

# SNP counts (raw)
n_snp_1000g_raw <- nrow(map_1000g)
n_snp_h3_raw    <- nrow(map_h3)

cat("Raw SNPs – 1000G:", n_snp_1000g_raw,
    "| H3:", n_snp_h3_raw, "\n\n")

# -------------------------------------------------------
# 1b) Load PLINK PCA results for H3 internal structure (Panels A & B)
# -------------------------------------------------------
cat("=== Loading PLINK PCA for H3 internal structure ===\n")

# Read PCA eigenvectors from PLINK - try multiple possible locations.
# Paths can be overridden via environment variables; defaults are relative to repo root.
plink_eigenvec_file <- Sys.getenv("PLINK_EIGENVEC", "results/pca/H3_QC9_noOutliers_PCA.eigenvec")
plink_eigenval_file <- Sys.getenv("PLINK_EIGENVAL", "results/pca/H3_QC9_noOutliers_PCA.eigenval")

# Also check current directory
if (!file.exists(plink_eigenvec_file)) {
  plink_eigenvec_file <- "H3_QC9_noOutliers_PCA.eigenvec"
}
if (!file.exists(plink_eigenval_file)) {
  plink_eigenval_file <- "H3_QC9_noOutliers_PCA.eigenval"
}

cat("Looking for PLINK PCA files:\n")
cat("  Eigenvec:", plink_eigenvec_file, "- Exists:", file.exists(plink_eigenvec_file), "\n")
cat("  Eigenval:", plink_eigenval_file, "- Exists:", file.exists(plink_eigenval_file), "\n\n")

if (file.exists(plink_eigenvec_file) && file.exists(plink_eigenval_file)) {
  # Read eigenvectors
  eig_plink <- fread(plink_eigenvec_file)
  colnames(eig_plink) <- c("FID", "IID", paste0("PC", 1:(ncol(eig_plink) - 2)))
  
  cat("✓ Read PLINK PCA eigenvec with", nrow(eig_plink), "samples and",
      ncol(eig_plink) - 2, "PCs.\n")
  cat("  First few sample IDs:\n")
  print(head(eig_plink[, c("FID", "IID")], 3))
  
  # Read eigenvalues and compute variance explained
  eigenval_plink <- fread(plink_eigenval_file, header = FALSE)$V1
  variance_explained_plink <- eigenval_plink / sum(eigenval_plink)
  
  pc1_var_plink <- round(variance_explained_plink[1] * 100, 2)
  pc2_var_plink <- round(variance_explained_plink[2] * 100, 2)
  
  pc1_lab_plink <- paste0("PC1 (", pc1_var_plink, "%)")
  pc2_lab_plink <- paste0("PC2 (", pc2_var_plink, "%)")
  
  cat("✓ PLINK PC1 variance:", pc1_var_plink, "%\n")
  cat("✓ PLINK PC2 variance:", pc2_var_plink, "%\n\n")
  
  has_plink_pca <- TRUE
} else {
  cat("✗ Warning: PLINK PCA files not found. Panels A & B will use bigsnpr PCA instead.\n")
  cat("  This means panels A & B will show global structure, not H3 internal structure.\n")
  cat("  Expected files:\n")
  cat("    -", plink_eigenvec_file, "\n")
  cat("    -", plink_eigenval_file, "\n\n")
  has_plink_pca <- FALSE
}

# -------------------------------------------------------
# 2) PCA projection with relaxed matching threshold
# -------------------------------------------------------
set.seed(1)

# Use multiple cores safely (OMP_NUM_THREADS is set at top of script)
res_pca <- bed_projectPCA(
  obj.bed.ref = bed_1000g,
  obj.bed.new = bed_h3,
  k           = 20,
  ncores      = 4,             # Now safe to use multiple cores
  strand_flip = TRUE,          # Allow strand flips
  match.min.prop = 0.1         # Accept if ≥10% variants match (default is 0.2)
)

# Reference SVD object
svd_ref <- res_pca$obj.svd.ref

# Proper 1000G PC scores: U * diag(d)
PC_ref_mat <- sweep(svd_ref$u, 2, svd_ref$d, `*`)

# OADP-projected PCs for H3 (already scaled)
PC_h3_mat  <- res_pca$OADP_proj

# SNP indices that were used in reference PCA
ind.col.ref.used <- attr(svd_ref, "subset")
n_snp_ld_pruned  <- length(ind.col.ref.used)

cat("SNPs used in reference PCA after LD/filters:", n_snp_ld_pruned, "\n\n")

# -------------------------------------------------------
# 3) Create PC dataframes
# -------------------------------------------------------
pc_ref <- data.frame(
  sample.id = id_1000g,
  PC1       = PC_ref_mat[, 1],
  PC2       = PC_ref_mat[, 2],
  PC3       = PC_ref_mat[, 3],
  stringsAsFactors = FALSE
)

pc_h3 <- data.frame(
  sample.id = id_h3,
  PC1       = PC_h3_mat[, 1],
  PC2       = PC_h3_mat[, 2],
  PC3       = PC_h3_mat[, 3],
  stringsAsFactors = FALSE
)

# -------------------------------------------------------
# 4) Variance explained (PVE) and scree plot
# -------------------------------------------------------
d_ref <- svd_ref$d
pve   <- 100 * (d_ref^2 / sum(d_ref^2))

cat("First 10 PCs variance explained:\n")
print(pve[1:10])
cat("Total variance:", sum(pve), "%\n\n")

# Scree plot for 1000G PCA
df_scree <- data.frame(
  PC  = seq_along(pve),
  PVE = pve
)

p_scree <- ggplot(df_scree[1:20, ], aes(x = PC, y = PVE)) +
  geom_line(color = "#0072B2", linewidth = 1) +
  geom_point(color = "#0072B2", size = 3, shape = 21, fill = "white", stroke = 1.5) +
  labs(
    title = "Variance Explained by Principal Components",
    x     = "Principal Component",
    y     = "Variance Explained (%)"
  ) +
  scale_x_continuous(breaks = seq(0, 20, 5)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme_nature()

print(p_scree)
save_plot(p_scree, "Scree_PC_variance_bigsnpr", width = 7, height = 5)

# -------------------------------------------------------
# 5) Load 1000G annotations
# -------------------------------------------------------
annot_1000g <- read.table(
  "1000_Genomes_Datasheet.txt",
  header           = TRUE,
  stringsAsFactors = FALSE
)

# Standardize column names
names(annot_1000g)[names(annot_1000g) == "sample"]    <- "sample.id"
names(annot_1000g)[names(annot_1000g) == "super_pop"] <- "superpop"

# Merge with PC data
pc_ref <- merge(
  pc_ref,
  annot_1000g,
  by    = "sample.id",
  all.x = TRUE
)

cat("Duplicated column names in pc_ref:", anyDuplicated(names(pc_ref)), "\n\n")

# -------------------------------------------------------
# 6) Combine datasets
# -------------------------------------------------------
pc_ref$dataset <- "1000G"
pc_h3$dataset  <- "Ghana"

pc_ref$superpop <- ifelse(is.na(pc_ref$superpop), "Unknown", pc_ref$superpop)
pc_h3$superpop  <- NA_character_

pc_all <- rbind(
  pc_ref[, c("sample.id", "PC1", "PC2", "PC3", "superpop", "dataset")],
  pc_h3[,  c("sample.id", "PC1", "PC2", "PC3", "superpop", "dataset")]
)

cat("Duplicated column names in pc_all:", anyDuplicated(names(pc_all)), "\n\n")

# -------------------------------------------------------
# 7) Interactive plots (plotly)
# -------------------------------------------------------
# 2D interactive: PC1 vs PC2
p_plotly_12 <- plot_ly(
  data = pc_all,
  x = ~PC1,
  y = ~PC2,
  color = ~dataset,
  text = ~sample.id,
  hoverinfo = "text",
  mode = "markers",
  type = "scatter",
  marker = list(size = 4, opacity = 0.7)
) %>%
  layout(
    title = "PCA PC1 vs PC2 (interactive)",
    xaxis = list(title = paste0("PC1 (", round(pve[1], 1), "%)")),
    yaxis = list(title = paste0("PC2 (", round(pve[2], 1), "%)"))
  )

print(p_plotly_12)

# 3D PCA: PC1–PC3
p_plotly_3d <- plot_ly(
  data = pc_all,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~dataset,
  text = ~sample.id,
  hoverinfo = "text",
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, opacity = 0.7)
) %>%
  layout(
    title = "3D PCA (interactive)"
  )

print(p_plotly_3d)

# -------------------------------------------------------
# 8) Static ggplot2 plots - Dataset comparison
# -------------------------------------------------------
p1 <- ggplot(pc_all, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point(alpha = 0.7, size = 2, shape = 16) +
  scale_color_manual(values = colors_dataset) +
  labs(
    title = "PCA: 1000 Genomes vs H3 Ghana",
    subtitle = "PC1 vs PC2",
    x = paste0("PC1 (", round(pve[1], 2), "%)"),
    y = paste0("PC2 (", round(pve[2], 2), "%)"),
    color = "Dataset"
  ) +
  theme_nature()
print(p1)
save_plot(p1, "PC1_PC2_dataset_bigsnpr", width = 7, height = 5)

p2 <- ggplot(pc_all, aes(x = PC1, y = PC3, color = dataset)) +
  geom_point(alpha = 0.7, size = 2, shape = 16) +
  scale_color_manual(values = colors_dataset) +
  labs(
    title = "PCA: 1000 Genomes vs H3 Ghana",
    subtitle = "PC1 vs PC3",
    x = paste0("PC1 (", round(pve[1], 2), "%)"),
    y = paste0("PC3 (", round(pve[3], 2), "%)"),
    color = "Dataset"
  ) +
  theme_nature()
print(p2)
save_plot(p2, "PC1_PC3_dataset_bigsnpr", width = 7, height = 5)

p3 <- ggplot(pc_all, aes(x = PC2, y = PC3, color = dataset)) +
  geom_point(alpha = 0.7, size = 2, shape = 16) +
  scale_color_manual(values = colors_dataset) +
  labs(
    title = "PCA: 1000 Genomes vs H3 Ghana",
    subtitle = "PC2 vs PC3",
    x = paste0("PC2 (", round(pve[2], 2), "%)"),
    y = paste0("PC3 (", round(pve[3], 2), "%)"),
    color = "Dataset"
  ) +
  theme_nature()
print(p3)
save_plot(p3, "PC2_PC3_dataset_bigsnpr", width = 7, height = 5)

# -------------------------------------------------------
# 9) Superpopulation plots with H3 overlay
# -------------------------------------------------------
# PC1 vs PC2
p_ref12 <- ggplot() +
  geom_point(
    data  = pc_ref,
    aes(x = PC1, y = PC2, color = superpop),
    alpha = 0.6,
    size  = 2,
    shape = 16
  ) +
  geom_point(
    data  = pc_h3,
    aes(x = PC1, y = PC2),
    shape = 17,  # Triangle
    size  = 2.5,
    color = "#D55E00",
    alpha = 0.8
  ) +
  scale_color_manual(values = colors_superpop) +
  labs(
    x     = paste0("PC1 (", round(pve[1], 2), "%)"),
    y     = paste0("PC2 (", round(pve[2], 2), "%)"),
    color = "Super Population"
  ) +
  theme_nature() +
  theme(legend.position = "right")
print(p_ref12)
save_plot(p_ref12, "PC1_PC2_superpop_with_H3_bigsnpr", width = 8, height = 6)

# PC1 vs PC3
p_ref13 <- ggplot() +
  geom_point(
    data  = pc_ref,
    aes(x = PC1, y = PC3, color = superpop),
    alpha = 0.6,
    size  = 2,
    shape = 16
  ) +
  geom_point(
    data  = pc_h3,
    aes(x = PC1, y = PC3),
    shape = 17,  # Triangle
    size  = 2.5,
    color = "#D55E00",
    alpha = 0.8
  ) +
  scale_color_manual(values = colors_superpop) +
  labs(
    x     = paste0("PC1 (", round(pve[1], 2), "%)"),
    y     = paste0("PC3 (", round(pve[3], 2), "%)"),
    color = "Super Population"
  ) +
  theme_nature() +
  theme(legend.position = "right")
print(p_ref13)
save_plot(p_ref13, "PC1_PC3_superpop_with_H3_bigsnpr", width = 8, height = 6)

# PC2 vs PC3
p_ref23 <- ggplot() +
  geom_point(
    data  = pc_ref,
    aes(x = PC2, y = PC3, color = superpop),
    alpha = 0.6,
    size  = 2,
    shape = 16
  ) +
  geom_point(
    data  = pc_h3,
    aes(x = PC2, y = PC3),
    shape = 17,  # Triangle
    size  = 2.5,
    color = "#D55E00",
    alpha = 0.8
  ) +
  scale_color_manual(values = colors_superpop) +
  labs(
    x     = paste0("PC2 (", round(pve[2], 2), "%)"),
    y     = paste0("PC3 (", round(pve[3], 2), "%)"),
    color = "Super Population"
  ) +
  theme_nature() +
  theme(legend.position = "right")
print(p_ref23)
save_plot(p_ref23, "PC2_PC3_superpop_with_H3_bigsnpr", width = 8, height = 6)

# -------------------------------------------------------
# 10) Combined dataset plots (superpop + shape by dataset)
# -------------------------------------------------------
# Ensure ancestry labels exist for 1000G
if (!"pop" %in% names(pc_ref)) {
  pc_ref$pop <- NA_character_
}
if (!"superpop" %in% names(pc_ref)) {
  pc_ref$superpop <- NA_character_
}

pc_ref$pop[is.na(pc_ref$pop)]           <- "Unknown_pop"
pc_ref$superpop[is.na(pc_ref$superpop)] <- "Unknown_superpop"

# Give H3 a population + superpop label
pc_h3_for_big <- pc_h3
pc_h3_for_big$pop      <- "Ghana"
pc_h3_for_big$superpop <- "Ghana"

# Combine into one data.frame
pc_big <- rbind(
  pc_ref[,       c("sample.id", "PC1", "PC2", "PC3", "pop", "superpop", "dataset")],
  pc_h3_for_big[, c("sample.id", "PC1", "PC2", "PC3", "pop", "superpop", "dataset")]
)

# PC1 vs PC2
p_big_12 <- ggplot(pc_big,
                   aes(x = PC1, y = PC2, color = superpop, shape = dataset)) +
  geom_point(alpha = 0.6, size = 2.5) +
  scale_shape_manual(values = c("1000G" = 16, "Ghana" = 17)) +
  scale_color_manual(values = colors_superpop) +
  labs(
    x     = paste0("PC1 (", round(pve[1], 2), "%)"),
    y     = paste0("PC2 (", round(pve[2], 2), "%)"),
    color = "Super Population",
    shape = "Dataset"
  ) +
  theme_nature() +
  theme(legend.position = "right")
print(p_big_12)
save_plot(p_big_12, "PC1_PC2_superpop_dataset_bigsnpr", width = 8, height = 6)

# PC1 vs PC3
p_big_13 <- ggplot(pc_big,
                   aes(x = PC1, y = PC3, color = superpop, shape = dataset)) +
  geom_point(alpha = 0.6, size = 2.5) +
  scale_shape_manual(values = c("1000G" = 16, "Ghana" = 17)) +
  scale_color_manual(values = colors_superpop) +
  labs(
    x     = paste0("PC1 (", round(pve[1], 2), "%)"),
    y     = paste0("PC3 (", round(pve[3], 2), "%)"),
    color = "Super Population",
    shape = "Dataset"
  ) +
  theme_nature() +
  theme(legend.position = "right")
print(p_big_13)
save_plot(p_big_13, "PC1_PC3_superpop_dataset_bigsnpr", width = 8, height = 6)

# PC2 vs PC3
p_big_23 <- ggplot(pc_big,
                   aes(x = PC2, y = PC3, color = superpop, shape = dataset)) +
  geom_point(alpha = 0.6, size = 2.5) +
  scale_shape_manual(values = c("1000G" = 16, "Ghana" = 17)) +
  scale_color_manual(values = colors_superpop) +
  labs(
    x     = paste0("PC2 (", round(pve[2], 2), "%)"),
    y     = paste0("PC3 (", round(pve[3], 2), "%)"),
    color = "Super Population",
    shape = "Dataset"
  ) +
  theme_nature() +
  theme(legend.position = "right")
print(p_big_23)
save_plot(p_big_23, "PC2_PC3_superpop_dataset_bigsnpr", width = 8, height = 6)

# -------------------------------------------------------
# 11) AFR zoom
# -------------------------------------------------------
pc_ref_afr <- subset(pc_ref, superpop == "AFR")

pc_h3_afr <- pc_h3
pc_h3_afr$pop      <- "Ghana"
pc_h3_afr$superpop <- "Ghana"
pc_h3_afr$dataset  <- "Ghana"

pc_afr_all <- rbind(
  pc_ref_afr[, c("sample.id", "PC1", "PC2", "PC3", "pop", "superpop", "dataset")],
  pc_h3_afr[,  c("sample.id", "PC1", "PC2", "PC3", "pop", "superpop", "dataset")]
)

# Calculate ranges for consistent zoom
xr_12 <- range(pc_afr_all$PC1)
yr_12 <- range(pc_afr_all$PC2)

xr_13 <- range(pc_afr_all$PC1)
yr_13 <- range(pc_afr_all$PC3)

xr_23 <- range(pc_afr_all$PC2)
yr_23 <- range(pc_afr_all$PC3)

# Get unique populations and assign colors
unique_pops <- unique(pc_afr_all$pop)
n_pops <- length(unique_pops)
afr_colors <- setNames(
  brewer.pal(min(n_pops, 12), "Set3"),
  unique_pops[1:min(n_pops, 12)]
)
# Ensure H3_Ghana has a distinct color
afr_colors["Ghana"] <- "#A020F0"

p_afr_12 <- ggplot(pc_afr_all,
                   aes(x = PC1, y = PC2, color = pop, shape = dataset)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_shape_manual(values = c("1000G" = 16, "Ghana" = 17)) +
  scale_color_manual(values = afr_colors) +
  coord_cartesian(xlim = xr_12, ylim = yr_12) +
  labs(
    x     = paste0("PC1 (", round(pve[1], 2), "%)"),
    y     = paste0("PC2 (", round(pve[2], 2), "%)"),
    color = "Population",
    shape = "Dataset"
  ) +
  theme_nature() +
  theme(legend.position = "right")
print(p_afr_12)
save_plot(p_afr_12, "AFR_zoom_PC1_PC2_pop_bigsnpr", width = 9, height = 6)

p_afr_13 <- ggplot(pc_afr_all,
                   aes(x = PC1, y = PC3, color = pop, shape = dataset)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_shape_manual(values = c("1000G" = 16, "Ghana" = 17)) +
  scale_color_manual(values = afr_colors) +
  coord_cartesian(xlim = xr_13, ylim = yr_13) +
  labs(
    x     = paste0("PC1 (", round(pve[1], 2), "%)"),
    y     = paste0("PC3 (", round(pve[3], 2), "%)"),
    color = "Population",
    shape = "Dataset"
  ) +
  theme_nature() +
  theme(legend.position = "right")
print(p_afr_13)
save_plot(p_afr_13, "AFR_zoom_PC1_PC3_pop_bigsnpr", width = 9, height = 6)

p_afr_23 <- ggplot(pc_afr_all,
                   aes(x = PC2, y = PC3, color = pop, shape = dataset)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_shape_manual(values = c("1000G" = 16, "Ghana" = 17)) +
  scale_color_manual(values = afr_colors) +
  coord_cartesian(xlim = xr_23, ylim = yr_23) +
  labs(
    x     = paste0("PC2 (", round(pve[2], 2), "%)"),
    y     = paste0("PC3 (", round(pve[3], 2), "%)"),
    color = "Population",
    shape = "Dataset"
  ) +
  theme_nature() +
  theme(legend.position = "right")
print(p_afr_23)
save_plot(p_afr_23, "AFR_zoom_PC2_PC3_pop_bigsnpr", width = 9, height = 6)

cat("\n=== SCRIPT COMPLETED SUCCESSFULLY ===\n")
cat("All plots saved to:\n")
cat("  - PCA_plots_bigsnpr/ (original directory)\n")
cat("  - BigSNP_pngs/ (PNG and PDF files)\n")

# -------------------------------------------------------
# 12) Combined Panels
# -------------------------------------------------------
cat("\n=== Creating Combined Panel Plots ===\n")

# Two-panel: Global + AFR zoom (side by side)
combined_2panel <- plot_grid(
  p_big_12, p_afr_12,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 16,
  label_fontface = "bold",
  rel_widths = c(1, 1),
  align = "h",
  axis = "tb"
)

# Save 2-panel combined plot
ggsave(
  filename = "BigSNP_pngs/Combined_Global_AFR_PC1_PC2.png",
  plot = combined_2panel,
  width = 16,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "BigSNP_pngs/Combined_Global_AFR_PC1_PC2.pdf",
  plot = combined_2panel,
  width = 16,
  height = 6,
  device = cairo_pdf
)

cat("Saved 2-panel combined plot (Global + AFR zoom)\n")

# -------------------------------------------------------
# 13) Load H3 ethnicity/strain data and create 4-panel plot
# -------------------------------------------------------
cat("\n=== Checking for H3 phenotype data ===\n")

# Try to load H3 phenotype file if it exists
h3_pheno_file <- "H3_pheno.phe.txt"

if (file.exists(h3_pheno_file)) {
  cat("Found H3 phenotype file. Creating 4-panel combined plot...\n")
  
  # Load phenotype data using fread (handles various formats better)
  h3_pheno <- fread(h3_pheno_file)
  
  cat("Successfully loaded phenotype file with", nrow(h3_pheno), "rows and", ncol(h3_pheno), "columns.\n")
  cat("Column names:", paste(names(h3_pheno), collapse = ", "), "\n")
  
  # -------------------------------------------------------
  # Decide which PCA to use for Panels A & B
  # -------------------------------------------------------
  if (has_plink_pca) {
    cat("\n=== Using PLINK PCA for Panels A & B (H3 internal structure) ===\n")
    
    cat("Phenotype file sample IDs (first 3):\n")
    print(head(h3_pheno[, c("FID", "IID")], 3))
    
    # Merge PLINK PCA with phenotype data
    dat_h3_internal <- eig_plink %>%
      left_join(h3_pheno, by = c("FID", "IID"))
    
    cat("✓ Merged PLINK PCA with phenotype data.\n")
    cat("  Rows after merge:", nrow(dat_h3_internal), "\n")
    cat("  Samples with ETHNICITY:", sum(!is.na(dat_h3_internal$ETHNICITY)), "\n")
    cat("  Samples with STRAIN:", sum(!is.na(dat_h3_internal$STRAIN)), "\n\n")
    
    # Use PLINK variance labels
    pc1_lab_internal <- pc1_lab_plink
    pc2_lab_internal <- pc2_lab_plink
    
  } else {
    cat("\n=== Using bigsnpr PCA for Panels A & B (PLINK PCA not available) ===\n")
    
    # Standardize column names
    if ("FID" %in% names(h3_pheno) && "IID" %in% names(h3_pheno)) {
      h3_pheno$sample.id <- h3_pheno$IID
    } else if ("IID" %in% names(h3_pheno)) {
      h3_pheno$sample.id <- h3_pheno$IID
    } else if ("FID" %in% names(h3_pheno)) {
      h3_pheno$sample.id <- h3_pheno$FID
    }
    
    # Merge bigsnpr PCA with phenotype data
    dat_h3_internal <- merge(pc_h3, h3_pheno, by = "sample.id", all.x = TRUE)
    
    cat("✓ Merged bigsnpr PCA with phenotype data.\n")
    cat("  Rows after merge:", nrow(dat_h3_internal), "\n\n")
    
    # Use bigsnpr variance labels
    pc1_lab_internal <- paste0("PC1 (", round(pve[1], 2), "%)")
    pc2_lab_internal <- paste0("PC2 (", round(pve[2], 2), "%)")
  }
  
  # -------------------------------------------------------
  # Process phenotype data (same for both)
  # -------------------------------------------------------
  
  # Recode ETHNICITY
  dat_h3_internal <- dat_h3_internal %>%
    mutate(
      ETHNICITY_CLEAN = case_when(
        !is.na(ETHNICITY) & tolower(ETHNICITY) == "akan" ~ "Akan",
        !is.na(ETHNICITY) & tolower(ETHNICITY) == "ewe"  ~ "Ewe",
        !is.na(ETHNICITY) & tolower(ETHNICITY) == "ga"   ~ "Ga",
        is.na(ETHNICITY) ~ NA_character_,
        TRUE ~ "Other"
      ),
      ETHNICITY_CLEAN = factor(ETHNICITY_CLEAN, levels = c("Akan", "Ewe", "Ga", "Other"))
    )
  
  # Recode STRAIN - keep only MAF and Mtbss
  dat_h3_internal <- dat_h3_internal %>%
    mutate(
      STRAIN_CLEAN = case_when(
        !is.na(STRAIN) & grepl("maf", STRAIN, ignore.case = TRUE) ~ "MAF",
        !is.na(STRAIN) & grepl("mtb", STRAIN, ignore.case = TRUE) ~ "Mtbss",
        TRUE ~ NA_character_
      ),
      STRAIN_CLEAN = factor(STRAIN_CLEAN, levels = c("MAF", "Mtbss"))
    )
  
  # Filter to keep only MAF and Mtbss
  pc_h3_filtered <- dat_h3_internal %>% filter(!is.na(STRAIN_CLEAN))
  
  # Remove outliers (IQR method)
  remove_outliers <- function(df, pc_col, multiplier = 3) {
    Q1 <- quantile(df[[pc_col]], 0.25, na.rm = TRUE)
    Q3 <- quantile(df[[pc_col]], 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1
    lower_bound <- Q1 - multiplier * IQR_val
    upper_bound <- Q3 + multiplier * IQR_val
    df %>% filter(.data[[pc_col]] >= lower_bound, .data[[pc_col]] <= upper_bound)
  }
  
  pc_h3_filtered <- pc_h3_filtered %>%
    remove_outliers("PC1", multiplier = 3) %>%
    remove_outliers("PC2", multiplier = 3)
  
  cat("H3 samples after filtering (MAF/Mtbss only, outliers removed):", nrow(pc_h3_filtered), "\n")
  
  # Define colors for ethnicity and strain
  ethnicity_colors <- c(
    "Akan"  = "#E69F00",  # Orange
    "Ewe"   = "#56B4E9",  # Sky blue
    "Ga"    = "#009E73",  # Bluish green
    "Other" = "#999999"   # Gray
  )
  
  strain_colors <- c(
    "MAF"   = "#5D4037",  # Dark brown
    "Mtbss" = "#FF0000"   # Bright red
  )
  
  # Plot A: PCA by Ethnicity
  p_eth <- pc_h3_filtered %>%
    filter(!is.na(ETHNICITY_CLEAN)) %>%
    ggplot(aes(x = PC1, y = PC2, color = ETHNICITY_CLEAN)) +
    geom_point(alpha = 0.7, size = 2.5, shape = 16) +
    scale_color_manual(values = ethnicity_colors) +
    labs(
      x = pc1_lab_internal,
      y = pc2_lab_internal,
      color = "Ethnicity"
    ) +
    theme_nature()
  
  # Plot B: PCA by Strain
  p_strain <- pc_h3_filtered %>%
    filter(!is.na(STRAIN_CLEAN)) %>%
    ggplot(aes(x = PC1, y = PC2, color = STRAIN_CLEAN)) +
    geom_point(alpha = 0.7, size = 2.5, shape = 16) +
    scale_color_manual(values = strain_colors) +
    labs(
      x = pc1_lab_internal,
      y = pc2_lab_internal,
      color = "Strain"
    ) +
    theme_nature()
  
  # Four-panel: H3 Ethnicity + H3 Strain + Global + AFR zoom
  combined_4panel <- plot_grid(
    p_eth, p_strain,
    p_big_12, p_afr_12,
    ncol = 2,
    nrow = 2,
    labels = c("A", "B", "C", "D"),
    label_size = 16,
    label_fontface = "bold",
    align = "hv",
    axis = "tblr"
  )
  
  # Save 4-panel combined plot
  ggsave(
    filename = "BigSNP_pngs/Combined_4panel_Ethnicity_Strain_Global_AFR.png",
    plot = combined_4panel,
    width = 16,
    height = 12,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = "BigSNP_pngs/Combined_4panel_Ethnicity_Strain_Global_AFR.pdf",
    plot = combined_4panel,
    width = 16,
    height = 12,
    device = cairo_pdf
  )
  
  cat("Saved 4-panel combined plot (Ethnicity + Strain + Global + AFR)\n")
  
} else {
  cat("H3 phenotype file not found. Skipping 4-panel combined plot.\n")
  cat("To create the 4-panel plot, ensure 'H3_phenotype.txt' exists in the working directory.\n")
}

cat("\n=== ALL COMBINED PANELS COMPLETED ===\n")