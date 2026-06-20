#!/usr/bin/env Rscript

################################################################################
# Visualize Case-Case GWAS Results
# Creates Manhattan plots, QQ plots, and regional plots
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

# -------------------------------------------------------
# Configuration
# -------------------------------------------------------
# Paths can be overridden via environment variables; defaults are relative to the
# repository root. e.g. GWAS_RESULTS=/path/to/file.txt Rscript figures/GWAS.R
RESULTS_FILE <- Sys.getenv("GWAS_RESULTS", "results/gwas/GWAS_MAF_vs_Mtbss_cleaned.txt")  # Or _logistic_cleaned.txt
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", "results/gwas/GWAS_plots")

cat("================================\n")
cat("GWAS Visualization\n")
cat("================================\n\n")

# Create output directory
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# -------------------------------------------------------
# Load results
# -------------------------------------------------------
cat("Loading GWAS results...\n")

# Try to find results file
if (!file.exists(RESULTS_FILE)) {
  # Try alternative names
  alt_files <- c(
    "GWAS_MAF_vs_Mtbss_firth_cleaned.txt",
    "GWAS_MAF_vs_Mtbss_logistic_cleaned.txt"
  )
  
  for (f in alt_files) {
    if (file.exists(f)) {
      RESULTS_FILE <- f
      break
    }
  }
}

if (!file.exists(RESULTS_FILE)) {
  stop("Cannot find GWAS results file!")
}

gwas <- fread(RESULTS_FILE)

# Standardize column names (PLINK2 vs PLINK1.9)
if ("#CHROM" %in% names(gwas)) {
  # PLINK2 format
  gwas <- gwas %>% 
    rename(
      CHR = `#CHROM`,
      BP = POS,
      SNP = ID
    )
  cat("✓ Detected PLINK2 format\n")
} else if ("CHR" %in% names(gwas) && "BP" %in% names(gwas)) {
  # PLINK 1.9 format - already standardized
  if (!"SNP" %in% names(gwas) && "ID" %in% names(gwas)) {
    gwas <- gwas %>% rename(SNP = ID)
  }
  cat("✓ Detected PLINK 1.9 format\n")
} else {
  stop("Cannot identify CHR/BP columns in results!")
}

# Check for required columns
if (!all(c("CHR", "BP", "P") %in% names(gwas))) {
  cat("Available columns:", paste(names(gwas), collapse = ", "), "\n")
  stop("Missing required columns: CHR, BP, or P")
}

# Add -log10(P)
gwas <- gwas %>%
  mutate(
    CHR = as.numeric(CHR),
    logP = -log10(P)
  ) %>%
  filter(!is.na(CHR) & CHR <= 22)  # Autosomal only

cat("✓ Loaded", nrow(gwas), "variants\n")
cat("  Columns available:", paste(names(gwas), collapse = ", "), "\n\n")

# -------------------------------------------------------
# Calculate genomic inflation
# -------------------------------------------------------
cat("=== GWAS Summary Statistics ===\n")

observed_chisq <- qchisq(gwas$P, df = 1, lower.tail = FALSE)
lambda <- median(observed_chisq, na.rm = TRUE) / qchisq(0.5, df = 1)

cat("Genomic inflation factor (λ):", round(lambda, 4), "\n")

if (lambda > 1.1) {
  cat("⚠️  High inflation detected! Check for population stratification.\n")
} else if (lambda < 0.95) {
  cat("⚠️  Low inflation detected! May be over-corrected.\n")
} else {
  cat("✓ Inflation within acceptable range\n")
}

n_gws <- sum(gwas$P < 5e-8)
n_sug <- sum(gwas$P < 1e-5)

cat("Genome-wide significant (P < 5×10⁻⁸):", n_gws, "\n")
cat("Suggestive (P < 1×10⁻⁵):", n_sug, "\n\n")

# -------------------------------------------------------
# Prepare data for Manhattan plot
# -------------------------------------------------------
cat("Preparing Manhattan plot data...\n")

# Calculate cumulative BP position
gwas_plot <- gwas %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP)) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  left_join(gwas, by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot)

# Axis setup
axis_set <- gwas_plot %>%
  group_by(CHR) %>%
  summarize(center = mean(BPcum))

# -------------------------------------------------------
# Manhattan Plot
# -------------------------------------------------------
cat("Creating Manhattan plot...\n")

sig_threshold <- -log10(5e-8)
sug_threshold <- -log10(1e-5)

p_manhattan <- ggplot(gwas_plot, aes(x = BPcum, y = logP, color = as.factor(CHR))) +
  geom_point(alpha = 0.75, size = 0.8) +
  scale_color_manual(values = rep(c("#0072B2", "#56B4E9"), 11)) +
  scale_x_continuous(
    label = axis_set$CHR,
    breaks = axis_set$center,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(gwas_plot$logP) * 1.05)) +
  geom_hline(yintercept = sig_threshold, color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = sug_threshold, color = "blue", linetype = "dashed", linewidth = 0.5) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(P)))
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

ggsave(
  file.path(OUTPUT_DIR, "Manhattan_plot.png"),
  p_manhattan,
  width = 12,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(OUTPUT_DIR, "Manhattan_plot.pdf"),
  p_manhattan,
  width = 12,
  height = 6,
  device = cairo_pdf
)

cat("✓ Saved Manhattan plot\n")

# -------------------------------------------------------
# QQ Plot
# -------------------------------------------------------
cat("Creating QQ plot...\n")

# Remove any infinite or invalid P-values
qq_data <- gwas %>%
  filter(P > 0 & P <= 1 & !is.na(P)) %>%
  arrange(P) %>%
  mutate(
    observed = -log10(P),
    expected = -log10(ppoints(n()))
  )

# Only add confidence intervals if we have enough data
if (nrow(qq_data) > 100) {
  qq_data <- qq_data %>%
    mutate(
      clower = -log10(qbeta(0.025, 1:n(), n():1)),
      cupper = -log10(qbeta(0.975, 1:n(), n():1))
    )
  
  has_ci <- TRUE
} else {
  has_ci <- FALSE
}

# Create plot
p_qq <- ggplot(qq_data, aes(x = expected, y = observed))

if (has_ci) {
  p_qq <- p_qq + 
    geom_ribbon(aes(ymin = clower, ymax = cupper), fill = "grey80", alpha = 0.5)
}

p_qq <- p_qq +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_point(alpha = 0.5, size = 1) +
  labs(
    x = expression(Expected~~-log[10](italic(P))),
    y = expression(Observed~~-log[10](italic(P)))
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  )

ggsave(
  file.path(OUTPUT_DIR, "QQ_plot.png"),
  p_qq,
  width = 6,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(OUTPUT_DIR, "QQ_plot.pdf"),
  p_qq,
  width = 6,
  height = 6,
  device = cairo_pdf
)

cat("✓ Saved QQ plot\n")

# -------------------------------------------------------
# Top hits table
# -------------------------------------------------------
if (n_sug > 0) {
  cat("\nCreating top hits table...\n")
  
  # Identify which columns are available
  available_cols <- c("CHR", "BP", "SNP", "ID", "A1", "REF", "ALT", 
                      "BETA", "OR", "SE", "P", "L95", "U95", 
                      "OR_95L", "OR_95U", "A1_FREQ", "OBS_CT")
  
  # Map PLINK2 column names to standard names
  if ("#CHROM" %in% names(gwas)) {
    gwas <- gwas %>%
      rename(
        CHR = `#CHROM`,
        BP = POS,
        SNP = ID,
        A1_FREQ = A1_FREQ
      )
  }
  
  # Select columns that exist
  select_cols <- intersect(available_cols, names(gwas))
  
  top_hits <- gwas %>%
    filter(P < 1e-5) %>%
    arrange(P) %>%
    select(all_of(select_cols)) %>%
    head(50)
  
  # Add effect direction if we have BETA or OR
  if ("BETA" %in% names(top_hits)) {
    top_hits <- top_hits %>%
      mutate(Effect_Direction = ifelse(BETA > 0, "MAF_risk", "Mtbss_risk"))
  } else if ("OR" %in% names(top_hits)) {
    top_hits <- top_hits %>%
      mutate(Effect_Direction = ifelse(OR > 1, "MAF_risk", "Mtbss_risk"))
  }
  
  fwrite(top_hits, file.path(OUTPUT_DIR, "top_50_hits.txt"),
         sep = "\t", quote = FALSE)
  
  cat("✓ Saved top hits table\n")
  
  # Print top 10
  cat("\n=== Top 10 Variants ===\n")
  print(head(top_hits, 10))
}

# -------------------------------------------------------
# Summary report
# -------------------------------------------------------
cat("\n================================\n")
cat("Visualization Complete!\n")
cat("================================\n")
cat("\nOutput files in:", OUTPUT_DIR, "/\n")
cat("  - Manhattan_plot.png/pdf\n")
cat("  - QQ_plot.png/pdf\n")
if (n_sug > 0) {
  cat("  - top_50_hits.txt\n")
}
cat("\n================================\n")