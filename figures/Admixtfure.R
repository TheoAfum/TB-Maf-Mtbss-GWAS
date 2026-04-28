#!/usr/bin/env Rscript

################################################################################
# ADMIXTURE Results Visualization - Nature Style
# Creates publication-quality plots of ADMIXTURE results
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(RColorBrewer)
  library(viridis)
})

# -------------------------------------------------------
# Configuration
# -------------------------------------------------------
INPUT_PREFIX <- "H3_QC8_ibdClean"
RESULTS_DIR <- "/home/tafum/H3_imputation/ADMIXTURE/No_merge/ADMIXTURE_results"
OUTPUT_DIR <- "/home/tafum/H3_imputation/ADMIXTURE/No_merge/ADMIXTURE_results/ADMIXTURE_plots"
MIN_K <- 2
MAX_K <- 10

# Optional: Load phenotype data for sorting/labeling
PHENO_FILE <- "/home/tafum/H3_imputation/final/GWAS/H3_pheno.phe.txt"  # Set to NULL if not available

cat("================================\n")
cat("ADMIXTURE Visualization\n")
cat("================================\n")
cat("Results directory:", RESULTS_DIR, "\n")
cat("K range:", MIN_K, "to", MAX_K, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("================================\n\n")

# Create output directory
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# -------------------------------------------------------
# 1) Parse Cross-Validation Errors
# -------------------------------------------------------
cat("=== Parsing Cross-Validation Errors ===\n")

cv_file <- file.path(RESULTS_DIR, "cv_errors.txt")
if (file.exists(cv_file)) {
  cv_data <- readLines(cv_file)
  
  # Extract K and CV error
  cv_df <- data.frame(
    K = as.integer(gsub(".*log(\\d+)\\.out.*", "\\1", cv_data)),
    CV_error = as.numeric(gsub(".*error.*: ([0-9.]+).*", "\\1", cv_data))
  ) %>%
    arrange(K)
  
  cat("Cross-validation errors:\n")
  print(cv_df)
  
  # Find best K
  best_k <- cv_df$K[which.min(cv_df$CV_error)]
  cat("\n✓ Best K (lowest CV error):", best_k, "\n\n")
  
} else {
  cat("Warning: cv_errors.txt not found. Skipping CV error plot.\n\n")
  cv_df <- NULL
  best_k <- NULL
}

# -------------------------------------------------------
# 2) Load FAM file for sample information
# -------------------------------------------------------
cat("=== Loading sample information ===\n")

fam_file <- file.path(RESULTS_DIR, paste0(INPUT_PREFIX, ".fam"))
fam <- fread(fam_file, header = FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")

cat("Loaded", nrow(fam), "samples from FAM file.\n")

# -------------------------------------------------------
# 3) Load phenotype data (if available)
# -------------------------------------------------------
if (!is.null(PHENO_FILE) && file.exists(PHENO_FILE)) {
  cat("Loading phenotype data from:", PHENO_FILE, "\n")
  
  pheno <- fread(PHENO_FILE)
  cat("Phenotype columns:", paste(names(pheno), collapse = ", "), "\n")
  
  # Merge with fam
  samples <- fam %>%
    left_join(pheno, by = c("FID", "IID"))
  
  # Check if we have grouping variables
  has_ethnicity <- "ETHNICITY" %in% names(samples)
  has_strain <- "STRAIN" %in% names(samples)
  
  if (has_ethnicity) {
    cat("✓ Found ETHNICITY column\n")
  }
  if (has_strain) {
    cat("✓ Found STRAIN column\n")
  }
  
} else {
  cat("No phenotype file provided. Using sample order from FAM file.\n")
  samples <- fam
  has_ethnicity <- FALSE
  has_strain <- FALSE
}

# Add sample index
samples$sample_idx <- 1:nrow(samples)

cat("\n")

# -------------------------------------------------------
# 4) Load ADMIXTURE Q matrices for all K
# -------------------------------------------------------
cat("=== Loading ADMIXTURE Q matrices ===\n")

q_list <- list()

for (k in MIN_K:MAX_K) {
  q_file <- file.path(RESULTS_DIR, paste0(INPUT_PREFIX, ".K", k, ".Q"))
  
  if (file.exists(q_file)) {
    q_mat <- fread(q_file, header = FALSE)
    colnames(q_mat) <- paste0("Pop", 1:k)
    
    # Add sample information
    q_df <- cbind(samples, q_mat)
    
    q_list[[paste0("K", k)]] <- q_df
    
    cat("✓ Loaded K =", k, "\n")
  } else {
    cat("✗ Warning: File not found:", q_file, "\n")
  }
}

cat("\n")

# -------------------------------------------------------
# 5) Define color palettes
# -------------------------------------------------------
# Create colorblind-friendly palettes for different K values
get_colors <- function(k) {
  if (k <= 8) {
    # Use Set2 for K <= 8
    brewer.pal(k, "Set2")
  } else if (k <= 12) {
    # Use Set3 for K <= 12
    brewer.pal(k, "Set3")
  } else {
    # Use viridis for larger K
    viridis(k, option = "D")
  }
}

# -------------------------------------------------------
# 6) Define Nature-style theme
# -------------------------------------------------------
theme_admixture <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      # Axis
      axis.title.x = element_text(size = base_size, face = "bold"),
      axis.title.y = element_text(size = base_size, face = "bold"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = base_size - 1, color = "black"),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      
      # Legend
      legend.position = "none",
      
      # Panel
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Title
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      
      # Margins
      plot.margin = margin(5, 5, 5, 5)
    )
}

# -------------------------------------------------------
# 7) Function to create ADMIXTURE bar plot
# -------------------------------------------------------
plot_admixture_k <- function(q_df, k, sort_by = NULL, add_dividers = FALSE) {
  
  # Reshape data for plotting
  q_long <- q_df %>%
    select(sample_idx, FID, IID, starts_with("Pop"), all_of(sort_by)) %>%
    pivot_longer(
      cols = starts_with("Pop"),
      names_to = "Population",
      values_to = "Ancestry"
    )
  
  # Sort samples if requested
  if (!is.null(sort_by) && sort_by %in% names(q_df)) {
    q_long <- q_long %>%
      arrange(.data[[sort_by]], sample_idx)
    
    # Update sample order
    sample_order <- unique(q_long$sample_idx)
    q_long$sample_order <- match(q_long$sample_idx, sample_order)
  } else {
    q_long$sample_order <- q_long$sample_idx
  }
  
  # Get colors
  colors <- get_colors(k)
  
  # Create plot
  p <- ggplot(q_long, aes(x = sample_order, y = Ancestry, fill = Population)) +
    geom_bar(stat = "identity", width = 1, position = "stack") +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(
      title = paste0("K = ", k),
      x = "Individuals",
      y = "Ancestry"
    ) +
    theme_admixture()
  
  # Add dividers between groups if requested
  if (add_dividers && !is.null(sort_by) && sort_by %in% names(q_df)) {
    group_data <- q_long %>%
      group_by(.data[[sort_by]]) %>%
      summarise(
        start = min(sample_order),
        end = max(sample_order),
        .groups = "drop"
      )
    
    # Add vertical lines between groups
    if (nrow(group_data) > 1) {
      for (i in 1:(nrow(group_data) - 1)) {
        p <- p + geom_vline(xintercept = group_data$end[i] + 0.5, 
                           color = "black", linewidth = 0.5)
      }
    }
  }
  
  return(p)
}

# -------------------------------------------------------
# 8) Create individual plots for each K
# -------------------------------------------------------
cat("=== Creating individual K plots ===\n")

for (k in MIN_K:MAX_K) {
  k_name <- paste0("K", k)
  
  if (k_name %in% names(q_list)) {
    q_df <- q_list[[k_name]]
    
    # Create plot without sorting
    p <- plot_admixture_k(q_df, k)
    
    # Save
    ggsave(
      filename = file.path(OUTPUT_DIR, paste0("Admixture_K", k, ".png")),
      plot = p,
      width = 10,
      height = 3,
      dpi = 600,
      bg = "white"
    )
    
    ggsave(
      filename = file.path(OUTPUT_DIR, paste0("Admixture_K", k, ".pdf")),
      plot = p,
      width = 10,
      height = 3,
      device = cairo_pdf
    )
    
    # If we have grouping variables, create sorted versions
    if (has_strain) {
      p_strain <- plot_admixture_k(q_df, k, sort_by = "STRAIN", add_dividers = TRUE)
      
      ggsave(
        filename = file.path(OUTPUT_DIR, paste0("Admixture_K", k, "_by_STRAIN.png")),
        plot = p_strain,
        width = 10,
        height = 3,
        dpi = 600,
        bg = "white"
      )
    }
    
    if (has_ethnicity) {
      p_eth <- plot_admixture_k(q_df, k, sort_by = "ETHNICITY", add_dividers = TRUE)
      
      ggsave(
        filename = file.path(OUTPUT_DIR, paste0("Admixture_K", k, "_by_ETHNICITY.png")),
        plot = p_eth,
        width = 10,
        height = 3,
        dpi = 600,
        bg = "white"
      )
    }
    
    cat("✓ Created plots for K =", k, "\n")
  }
}

cat("\n")

# -------------------------------------------------------
# 9) Create combined multi-panel plot (K=2 to K=10)
# -------------------------------------------------------
cat("=== Creating combined multi-K panel ===\n")

plot_list <- list()

for (k in MIN_K:MAX_K) {
  k_name <- paste0("K", k)
  
  if (k_name %in% names(q_list)) {
    q_df <- q_list[[k_name]]
    
    # Decide how to sort
    sort_by <- NULL
    if (has_strain) {
      sort_by <- "STRAIN"
    } else if (has_ethnicity) {
      sort_by <- "ETHNICITY"
    }
    
    p <- plot_admixture_k(q_df, k, sort_by = sort_by, add_dividers = TRUE)
    plot_list[[k_name]] <- p
  }
}

if (length(plot_list) > 0) {
  # Combine all K plots
  combined <- plot_grid(
    plotlist = plot_list,
    ncol = 1,
    align = "v",
    axis = "lr"
  )
  
  # Save combined plot
  ggsave(
    filename = file.path(OUTPUT_DIR, "Admixture_Combined_K2_K10.png"),
    plot = combined,
    width = 12,
    height = 2.5 * length(plot_list),
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(OUTPUT_DIR, "Admixture_Combined_K2_K10.pdf"),
    plot = combined,
    width = 12,
    height = 2.5 * length(plot_list),
    device = cairo_pdf
  )
  
  cat("✓ Created combined multi-K plot\n\n")
}

# -------------------------------------------------------
# 10) Plot Cross-Validation Error
# -------------------------------------------------------
if (!is.null(cv_df)) {
  cat("=== Creating CV error plot ===\n")
  
  p_cv <- ggplot(cv_df, aes(x = K, y = CV_error)) +
    geom_line(color = "#0072B2", linewidth = 1) +
    geom_point(color = "#0072B2", size = 3, shape = 21, fill = "white", stroke = 1.5) +
    geom_point(
      data = cv_df[cv_df$K == best_k, ],
      aes(x = K, y = CV_error),
      color = "#D55E00",
      size = 4,
      shape = 21,
      fill = "#D55E00"
    ) +
    scale_x_continuous(breaks = MIN_K:MAX_K) +
    labs(
      title = "ADMIXTURE Cross-Validation Error",
      x = "Number of ancestral populations (K)",
      y = "Cross-validation error"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = file.path(OUTPUT_DIR, "CV_error_plot.png"),
    plot = p_cv,
    width = 7,
    height = 5,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(OUTPUT_DIR, "CV_error_plot.pdf"),
    plot = p_cv,
    width = 7,
    height = 5,
    device = cairo_pdf
  )
  
  cat("✓ Created CV error plot\n")
  cat("  Best K (highlighted in red):", best_k, "\n\n")
}

# -------------------------------------------------------
# 11) Create ancestry proportion summary table
# -------------------------------------------------------
cat("=== Creating ancestry proportion summaries ===\n")

if (!is.null(best_k) && paste0("K", best_k) %in% names(q_list)) {
  q_best <- q_list[[paste0("K", best_k)]]
  
  # Calculate mean ancestry proportions
  ancestry_summary <- q_best %>%
    select(starts_with("Pop")) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(everything(), names_to = "Population", values_to = "Mean_Proportion") %>%
    mutate(Percentage = round(Mean_Proportion * 100, 2)) %>%
    arrange(desc(Mean_Proportion))
  
  cat("\nMean ancestry proportions for best K (K =", best_k, "):\n")
  print(ancestry_summary)
  
  # Save to file
  write.csv(
    ancestry_summary,
    file = file.path(OUTPUT_DIR, paste0("Ancestry_proportions_K", best_k, ".csv")),
    row.names = FALSE
  )
  
  cat("\n✓ Saved ancestry proportions summary\n")
}

# -------------------------------------------------------
# Summary
# -------------------------------------------------------
cat("\n================================\n")
cat("ADMIXTURE visualization complete!\n")
cat("================================\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("\nFiles created:\n")
cat("  - Individual K plots (K2 to K10)\n")
cat("  - Combined multi-K panel\n")
if (!is.null(cv_df)) {
  cat("  - Cross-validation error plot\n")
  cat("  - Best K:", best_k, "\n")
}
if (has_strain || has_ethnicity) {
  cat("  - Sorted plots by phenotype\n")
}
cat("\n================================\n")
