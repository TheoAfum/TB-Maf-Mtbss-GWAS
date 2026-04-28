#!/usr/bin/env Rscript

################################################################################
# Complete Package Installer for PCA Analysis
# Installs all dependencies in correct order
################################################################################

cat("================================\n")
cat("Installing All Required Packages\n")
cat("This may take 10-15 minutes\n")
cat("================================\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Function to install with error handling
safe_install <- function(pkg) {
  cat("\n--- Installing", pkg, "---\n")
  
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "already installed\n")
    return(TRUE)
  }
  
  tryCatch({
    install.packages(pkg, dependencies = TRUE, quiet = FALSE)
    library(pkg, character.only = TRUE)
    cat("✓", pkg, "successfully installed\n")
    return(TRUE)
  }, error = function(e) {
    cat("✗ Failed to install", pkg, "\n")
    cat("  Error:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# -------------------------------------------------------
# Phase 1: Core system dependencies (MUST be first!)
# -------------------------------------------------------
cat("\n=== Phase 1: Core Dependencies ===\n")

core_deps <- c(
  "stringi",      # Required by stringr and many others
  "stringr",      # String manipulation
  "S7",           # Required by bigstatsr
  "cli",          # Command line interface
  "glue",         # String interpolation
  "rlang",        # R language tools
  "lifecycle",    # Lifecycle management
  "vctrs",        # Vector tools
  "pillar",       # Table printing
  "tibble",       # Modern data frames
  "fansi",        # ANSI escape codes
  "utf8",         # UTF-8 handling
  "crayon"        # Colored terminal output
)

for (pkg in core_deps) {
  safe_install(pkg)
  Sys.sleep(0.5)  # Brief pause between installs
}

# -------------------------------------------------------
# Phase 2: Data manipulation packages
# -------------------------------------------------------
cat("\n=== Phase 2: Data Manipulation ===\n")

data_pkgs <- c(
  "magrittr",     # Pipe operator
  "dplyr",        # Data manipulation
  "tidyr",        # Data tidying
  "purrr",        # Functional programming
  "data.table",   # Fast data manipulation
  "readr"         # Reading data
)

for (pkg in data_pkgs) {
  safe_install(pkg)
  Sys.sleep(0.5)
}

# -------------------------------------------------------
# Phase 3: Plotting packages
# -------------------------------------------------------
cat("\n=== Phase 3: Plotting Packages ===\n")

plot_pkgs <- c(
  "scales",       # Scaling functions
  "RColorBrewer", # Color palettes
  "viridis",      # Color scales
  "ggplot2",      # Plotting
  "cowplot",      # Plot arrangements
  "patchwork",    # Combining plots
  "plotly",       # Interactive plots
  "htmlwidgets",  # HTML widgets
  "htmltools"     # HTML tools
)

for (pkg in plot_pkgs) {
  safe_install(pkg)
  Sys.sleep(0.5)
}

# -------------------------------------------------------
# Phase 4: Statistical packages
# -------------------------------------------------------
cat("\n=== Phase 4: Statistical Packages ===\n")

stats_pkgs <- c(
  "foreach",      # Parallel loops
  "doParallel",   # Parallel backend
  "Rcpp",         # C++ interface
  "RcppArmadillo",# Armadillo C++ library
  "Matrix",       # Sparse matrices
  "bigmemory",    # Large datasets
  "bigstatsr",    # Big data statistics
  "bigparallelr"  # Parallel processing
)

for (pkg in stats_pkgs) {
  safe_install(pkg)
  Sys.sleep(0.5)
}

# -------------------------------------------------------
# Phase 5: Specialized dependencies for bigsnpr
# -------------------------------------------------------
cat("\n=== Phase 5: bigsnpr Dependencies ===\n")

bigsnpr_deps <- c(
  "checkmate",    # Parameter validation
  "htmlTable",    # HTML tables
  "Hmisc",        # Harrell miscellaneous
  "pcadapt"       # PCA adaptation
)

for (pkg in bigsnpr_deps) {
  safe_install(pkg)
  Sys.sleep(0.5)
}

# -------------------------------------------------------
# Phase 6: Finally install bigsnpr
# -------------------------------------------------------
cat("\n=== Phase 6: Installing bigsnpr ===\n")

safe_install("bigsnpr")

# -------------------------------------------------------
# Verification
# -------------------------------------------------------
cat("\n================================\n")
cat("Verifying Installation\n")
cat("================================\n\n")

test_packages <- c(
  "stringi", "S7", "scales", "bigstatsr", "bigsnpr",
  "ggplot2", "dplyr", "data.table", "cowplot"
)

all_ok <- TRUE

for (pkg in test_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "\n")
  } else {
    cat("✗", pkg, "FAILED\n")
    all_ok <- FALSE
  }
}

cat("\n================================\n")

if (all_ok) {
  cat("✓ ALL PACKAGES INSTALLED SUCCESSFULLY!\n")
  cat("================================\n\n")
  cat("You can now run:\n")
  cat("  Rscript pca_analysis_fixed.R\n\n")
} else {
  cat("✗ SOME PACKAGES FAILED\n")
  cat("================================\n\n")
  cat("Please try installing failed packages manually:\n")
  cat("  install.packages('package_name', dependencies=TRUE)\n\n")
}
