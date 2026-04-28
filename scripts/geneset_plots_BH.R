# =============================================================================
# MAGMA Gene Set Analysis — Forest Plot with BH Correction (Nature-level figure)
#
# Reads from MAGMA results with BH correction included.
# Specifically designed for corrected output files from magma_summary_with_BH.sh
#
# Usage:
#   Rscript gsea_forest_plot_BH.R                          # uses default filename
#   Rscript gsea_forest_plot_BH.R geneset_results_with_corrections.txt
#   Rscript gsea_forest_plot_BH.R geneset_results_with_corrections.txt out_prefix
#
# Requires: ggplot2, dplyr
#   install.packages(c("ggplot2", "dplyr"))
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# 0. Config — edit defaults here or pass via command-line arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

GSA_FILE   <- if (length(args) >= 1) args[1] else "geneset_results_with_corrections.txt"
OUT_PREFIX <- if (length(args) >= 2) args[2] else "MAGMA_BH_corrected"
FIG_WIDTH  <- 10     # inches (wider for better visualization)
FIG_HEIGHT <- NULL   # NULL = auto-scale to number of gene sets shown
P_THRESH   <- 0.05   # significance threshold for colour coding
BH_THRESH  <- 0.05   # BH FDR threshold for highlighting
CI_MULT    <- 1.96   # 1.96 = 95% CI; use 2.576 for 99%
TOP_N      <- 10     # gene sets to show in main figure (ranked by BH-corrected p-value)
                     # set to Inf to plot all; rest go to supplementary CSV

# -----------------------------------------------------------------------------
# 1. Read corrected MAGMA results
# -----------------------------------------------------------------------------
if (!file.exists(GSA_FILE)) {
  stop(
    "File not found: ", GSA_FILE,
    "\nPlease provide the path to your corrected MAGMA results file.",
    "\nExpected from: magma_summary_with_BH.sh output"
  )
}

message("Reading: ", GSA_FILE)

gsea_raw <- tryCatch(
  read.table(
    GSA_FILE,
    header           = TRUE,
    sep              = "\t",           # tab-delimited from our script
    stringsAsFactors = FALSE,
    fill             = TRUE,
    quote            = ""
  ),
  error = function(e) stop("Could not parse file: ", e$message)
)

message("  Loaded ", nrow(gsea_raw), " rows")

# -----------------------------------------------------------------------------
# 2. Validate expected columns (including BH corrections)
# -----------------------------------------------------------------------------
required_cols <- c("VARIABLE", "NGENES", "BETA", "SE", "P", "P_BH")
missing_cols  <- setdiff(required_cols, colnames(gsea_raw))

if (length(missing_cols) > 0) {
  stop(
    "Missing expected column(s): ", paste(missing_cols, collapse = ", "),
    "\nColumns found: ", paste(colnames(gsea_raw), collapse = ", "),
    "\nThis script expects output from magma_summary_with_BH.sh"
  )
}

# Keep only SET rows (drops COVAR rows if present in the file)
if ("TYPE" %in% colnames(gsea_raw)) {
  gsea_raw <- filter(gsea_raw, TYPE == "SET")
  message("  After filtering TYPE == SET: ", nrow(gsea_raw), " gene sets")
}

# -----------------------------------------------------------------------------
# 3. Derive plot variables with BH correction emphasis
# -----------------------------------------------------------------------------
gsea <- gsea_raw %>%
  mutate(
    NGENES = as.integer(NGENES),
    BETA   = as.numeric(BETA),
    SE     = as.numeric(SE),
    P      = as.numeric(P),
    P_BH   = as.numeric(P_BH),

    # Human-readable y-axis labels with proper capitalization
    label = gsub("_", " ", VARIABLE),
    label = case_when(
      # Fix specific pathway names that need custom formatting
      grepl("JAK.*STAT", label, ignore.case = TRUE) ~ "JAK-STAT Signaling",
      grepl("MHC.*CLASS", label, ignore.case = TRUE) ~ "MHC Class I",
      grepl("^B.*CELL", label, ignore.case = TRUE) ~ "B Cell Function",
      grepl("^T.*CELL", label, ignore.case = TRUE) ~ "T Cell Activation", 
      grepl("TH1.*RESPONSE", label, ignore.case = TRUE) ~ "Th1 Response",
      grepl("TH2.*RESPONSE", label, ignore.case = TRUE) ~ "Th2 Response",
      grepl("TREG.*RESPONSE", label, ignore.case = TRUE) ~ "Treg Response",
      grepl("VITAMIN.*D", label, ignore.case = TRUE) ~ "Vitamin D Signaling",
      grepl("INFLAMMASOME", label, ignore.case = TRUE) ~ "Inflammasome",
      grepl("CELL.*ADHESION", label, ignore.case = TRUE) ~ "Cell Adhesion",
      # Default formatting for other pathways
      TRUE ~ {
        temp_label <- gsub("SIGNALING", "Signaling", label)
        temp_label <- gsub("PATHWAY", "Pathway", temp_label)
        temp_label <- gsub("RECEPTOR", "Receptor", temp_label)
        temp_label <- gsub("RESPONSE", "Response", temp_label)
        temp_label <- gsub("FUNCTION", "Function", temp_label)
        temp_label <- gsub("ACTIVATION", "Activation", temp_label)
        tools::toTitleCase(tolower(temp_label))
      }
    ),

    # Confidence intervals
    ci_lo = BETA - CI_MULT * SE,
    ci_hi = BETA + CI_MULT * SE,

    # Enhanced significance tier (prioritizes BH correction)
    sig_level = case_when(
      P_BH < 0.001              ~ "BH q < 0.001",
      P_BH < BH_THRESH          ~ paste0("BH q < ", BH_THRESH),
      P < 0.001                 ~ "p < 0.001",
      P < P_THRESH              ~ paste0("p < ", P_THRESH),
      TRUE                      ~ "ns"
    ),
    sig_level = factor(
      sig_level,
      levels = c("BH q < 0.001", paste0("BH q < ", BH_THRESH), 
                 "p < 0.001", paste0("p < ", P_THRESH), "ns")
    ),
    
    # Flag for BH significant results
    BH_significant = P_BH < BH_THRESH
  ) %>%
  arrange(P_BH) %>%  # Sort by BH-corrected p-value
  slice_head(n = min(TOP_N, nrow(.))) %>%  # Take top N by BH p-value
  arrange(BETA) %>%  # Then sort by effect size for forest plot
  mutate(label = factor(label, levels = unique(label)))   # lock sort order

# -----------------------------------------------------------------------------
# 4. Enhanced supplementary table with BH correction
# -----------------------------------------------------------------------------
supp_table <- gsea_raw %>%
  arrange(P_BH) %>%  # Sort by BH-corrected p-value
  transmute(
    `Gene set`     = gsub("_", " ", VARIABLE),
    `N genes`      = NGENES,
    `Beta`         = round(BETA, 4),
    `SE`           = round(SE,   4),
    `95% CI lower` = round(BETA - CI_MULT * SE, 4),
    `95% CI upper` = round(BETA + CI_MULT * SE, 4),
    `P nominal`    = signif(P,    4),
    `P BH (FDR)`   = signif(P_BH, 4),
    `BH significant` = P_BH < BH_THRESH
  )

supp_out <- paste0(OUT_PREFIX, "_supplementary_all_genesets_BH.csv")
write.csv(supp_table, supp_out, row.names = FALSE)
message("Supplementary table -> ", supp_out, "  (", nrow(supp_table), " gene sets)")

message("Main figure: showing top ", nrow(gsea), " gene sets by BH-corrected p-value")

# Report significant findings
n_BH_sig <- sum(gsea$BH_significant)
n_nominal_sig <- sum(gsea$P < P_THRESH & !gsea$BH_significant)
message("  BH FDR significant (q < ", BH_THRESH, "): ", sum(gsea_raw$P_BH < BH_THRESH), " gene set(s)")
message("  Nominal significant (p < ", P_THRESH, "): ", sum(gsea_raw$P < P_THRESH), " gene set(s)")

# -----------------------------------------------------------------------------
# 5. Enhanced colour + alpha palettes (emphasizing BH correction)
# -----------------------------------------------------------------------------
sig_colours <- setNames(
  c("#8B0000", "#C0392B", "#E67E22", "#F39C12", "#95A5A6"),
  c("BH q < 0.001", paste0("BH q < ", BH_THRESH), "p < 0.001", 
    paste0("p < ", P_THRESH), "ns")
)

sig_alpha <- setNames(
  c(1.00, 1.00, 0.85, 0.75, 0.45),
  c("BH q < 0.001", paste0("BH q < ", BH_THRESH), "p < 0.001", 
    paste0("p < ", P_THRESH), "ns")
)

# Annotations for significant rows (prioritize BH correction)
ann <- gsea %>%
  filter(P_BH < BH_THRESH | P < P_THRESH) %>%
  mutate(
    p_label = case_when(
      P_BH < 0.001 ~ "q < 0.001",
      P_BH < BH_THRESH ~ paste0("q = ", formatC(P_BH, digits = 3, format = "f")),
      P < 0.001 ~ "p < 0.001",
      TRUE ~ paste0("p = ", formatC(P, digits = 3, format = "f"))
    )
  )

n_sig <- nrow(ann)
message("  Showing significance labels: ", n_sig, " gene set(s)")

# -----------------------------------------------------------------------------
# 6. Axis limits with room for BH annotations (adjusted for wider figure)
# -----------------------------------------------------------------------------
x_range <- range(c(gsea$ci_lo, gsea$ci_hi), na.rm = TRUE)
x_pad   <- diff(x_range) * 0.06    # slightly less padding since figure is wider
ann_pad <- if (n_sig > 0) 0.30 else 0.04    # adjusted for wider format
x_lims  <- c(x_range[1] - x_pad, x_range[2] + x_pad + ann_pad)
x_breaks <- pretty(c(gsea$ci_lo, gsea$ci_hi), n = 8)  # more breaks for wider figure

# Auto height: longer scale for better screen fit
if (is.null(FIG_HEIGHT)) {
  FIG_HEIGHT <- max(5.5, nrow(gsea) * 0.5 + 2.5)
}

# -----------------------------------------------------------------------------
# 7. Build enhanced plot with BH emphasis
# -----------------------------------------------------------------------------
p <- ggplot(gsea, aes(x = BETA, y = label)) +

  # Alternating row shading with subtle BH highlight -----------------------
  geom_hline(
    data      = filter(gsea, as.integer(label) %% 2 == 0),
    aes(yintercept = as.integer(label)),
    colour    = "#F8F8F8",
    linewidth = 9.5
  ) +
  
  # Highlight BH significant rows with subtle background
  geom_hline(
    data      = filter(gsea, BH_significant),
    aes(yintercept = as.integer(label)),
    colour    = "#FFF9E6",
    linewidth = 9.5,
    alpha     = 0.8
  ) +

  # Null-effect reference line -----------------------------------------------
  geom_vline(
    xintercept = 0, linetype = "dashed",
    colour = "#AAAAAA", linewidth = 0.5
  ) +

  # Confidence interval whiskers ---------------------------------------------
  geom_errorbarh(
    aes(xmin = ci_lo, xmax = ci_hi,
        colour = sig_level, alpha = sig_level),
    height    = 0.4,
    linewidth = 0.6
  ) +

  # Points (size proportional to gene-set size) ------------------------------
  geom_point(
    aes(size = NGENES, fill = sig_level, alpha = sig_level),
    shape  = 21,
    colour = "white",
    stroke = 0.7
  ) +

  # Enhanced inline labels for significant hits ------------------------------
  {if (n_sig > 0)
    geom_text(
      data    = ann,
      aes(x = ci_hi + 0.05, label = p_label, colour = sig_level),
      hjust   = 0,
      size    = 2.8,
      fontface = "bold",
      family  = "Arial"
    )
  } +

  # Scales -------------------------------------------------------------------
  scale_colour_manual(values = sig_colours, name = "Significance") +
  scale_fill_manual  (values = sig_colours, name = "Significance") +
  scale_alpha_manual (values = sig_alpha,   name = "Significance") +
  scale_size_continuous(
    range  = c(2.5, 8),
    breaks = c(5, 15, 25, 35, 50),
    name   = "Gene set size"
  ) +
  scale_x_continuous(
    limits = x_lims,
    breaks = x_breaks,
    expand = c(0, 0)
  ) +

  # Labels -------------------------------------------------------------------
  labs(
    x = expression("Effect size (" * beta * ")"),
    y = NULL
  ) +

  # Enhanced theme (Nature style) -------------------------------------------
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white", colour = NA),
    panel.grid.major.x = element_line(colour = "#EEEEEE", linewidth = 0.35),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),

    axis.line.x   = element_line(colour = "#333333", linewidth = 0.6),
    axis.line.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.ticks.x  = element_line(colour = "#555555", linewidth = 0.45),
    axis.text.y   = element_text(size = 9,  colour = "#222222", hjust = 1,
                                  family = "Arial"),
    axis.text.x   = element_text(size = 9,  colour = "#333333", family = "Arial"),
    axis.title.x  = element_text(size = 10, colour = "#111111", family = "Arial",
                                  margin = margin(t = 8)),

    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin   = margin(15, 120, 12, 10),  # increased right margin for legend positioning

    legend.position    = "right",
    legend.justification = c(0, 0.5),  # position legend further right
    legend.title       = element_text(size = 9, face = "bold", colour = "#222222"),
    legend.text        = element_text(size = 8, colour = "#333333"),
    legend.key.size    = unit(0.9, "lines"),
    legend.background  = element_blank(),
    legend.box.spacing = unit(8, "pt"),  # more spacing from plot
    legend.margin      = margin(l = 20)  # push legend further right
  ) +

  guides(
    colour = guide_legend(order = 1, override.aes = list(size = 4)),
    fill   = guide_legend(order = 1, override.aes = list(size = 4)),
    alpha  = "none",
    size   = guide_legend(
      order = 2,
      override.aes = list(fill = "#888888", colour = "white", alpha = 0.8)
    )
  )

# -----------------------------------------------------------------------------
# 8. Enhanced export with BH correction notation
# -----------------------------------------------------------------------------
pdf_out <- paste0(OUT_PREFIX, "_forest_plot_BH_corrected.pdf")
png_out <- paste0(OUT_PREFIX, "_forest_plot_BH_corrected.png")

ggsave(pdf_out, plot = p, width = FIG_WIDTH, height = FIG_HEIGHT,
       dpi = 300, device = cairo_pdf)
message("Saved PDF  → ", pdf_out)

ggsave(png_out, plot = p, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = 300)
message("Saved PNG  → ", png_out)

# -----------------------------------------------------------------------------
# 9. Summary statistics
# -----------------------------------------------------------------------------
message("\n=== SUMMARY ===")
message("Total gene sets tested: ", nrow(gsea_raw))
message("BH FDR significant (q < ", BH_THRESH, "): ", sum(gsea_raw$P_BH < BH_THRESH))
message("Nominal significant (p < ", P_THRESH, "): ", sum(gsea_raw$P < P_THRESH))

if(sum(gsea_raw$P_BH < BH_THRESH) > 0) {
  message("\nBH FDR SIGNIFICANT GENE SETS:")
  bh_sig <- gsea_raw %>%
    filter(P_BH < BH_THRESH) %>%
    arrange(P_BH) %>%
    select(VARIABLE, BETA, P, P_BH)
  print(bh_sig)
}

message("\nFiles created:")
message("  • ", pdf_out)
message("  • ", png_out) 
message("  • ", supp_out)