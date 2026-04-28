#!/usr/bin/env Rscript

################################################################################
# Pathway Enrichment Analysis - With Comprehensive Visualizations
# Converts MAGMA gene IDs to symbols using biomaRt
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# -------------------------------------------------------
# Configuration
# -------------------------------------------------------
MAGMA_GENE_FILE <- "/home/tafum/H3_imputation/final/GWAS/Final_GWAS/GWAS_plots/MAGMA_analysis/genome_wide_gene_analysis.genes.out"
OUTPUT_DIR <- "/home/tafum/H3_imputation/final/GWAS/Final_GWAS/GWAS_plots/Pathway_analysis"

cat("================================\n")
cat("Pathway Enrichment Analysis\n")
cat("Using MAGMA Gene Results\n")
cat("================================\n\n")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# -------------------------------------------------------
# Install required packages
# -------------------------------------------------------
cat("Checking required packages...\n")

if (!require("biomaRt", quietly = TRUE)) {
  cat("Installing biomaRt...\n")
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("biomaRt", update = FALSE, ask = FALSE)
}

if (!require("enrichR", quietly = TRUE)) {
  cat("Installing enrichR...\n")
  install.packages("enrichR")
}

library(biomaRt)
library(enrichR)
library(ggplot2)
library(cowplot)
library(scales)
library(RColorBrewer)
library(forcats)

cat("✓ Packages loaded\n\n")

# -------------------------------------------------------
# Load MAGMA results
# -------------------------------------------------------
cat("Loading MAGMA gene results...\n")

genes_all <- fread(MAGMA_GENE_FILE)
genes_all <- as.data.frame(genes_all)

# Get significant genes
sig_genes <- genes_all %>%
  filter(P < 0.05) %>%
  arrange(P)

sig_genes <- as.data.frame(sig_genes)

cat("✓ Found", nrow(sig_genes), "significant genes (P<0.05)\n\n")

# -------------------------------------------------------
# Convert Entrez IDs to gene symbols
# -------------------------------------------------------
cat("Converting gene IDs to symbols using biomaRt...\n")
cat("This may take a few minutes...\n\n")

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

cat("✓ Connected to Ensembl\n")

gene_ids <- sig_genes$GENE
batch_size <- 500
n_batches <- ceiling(length(gene_ids) / batch_size)

gene_symbols_list <- list()

for (i in 1:n_batches) {
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, length(gene_ids))

  batch_ids <- gene_ids[start_idx:end_idx]

  cat("  Processing batch", i, "of", n_batches,
      "(", length(batch_ids), "genes)...\n")

  tryCatch({
    gene_conversion <- getBM(
      attributes = c("entrezgene_id", "hgnc_symbol"),
      filters = "entrezgene_id",
      values = batch_ids,
      mart = ensembl
    )

    gene_symbols_list[[i]] <- gene_conversion
    Sys.sleep(1)

  }, error = function(e) {
    cat("    Warning: Batch", i, "failed:", conditionMessage(e), "\n")
  })
}

gene_conversion_all <- bind_rows(gene_symbols_list)
gene_conversion_all <- as.data.frame(gene_conversion_all)

gene_symbols <- gene_conversion_all %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
  pull(hgnc_symbol) %>%
  unique()

cat("\n✓ Converted", length(gene_symbols), "gene symbols\n\n")

# Save conversion table
gene_conversion_table <- sig_genes %>%
  left_join(gene_conversion_all, by = c("GENE" = "entrezgene_id")) %>%
  dplyr::select(EntrezID = GENE, Symbol = hgnc_symbol, CHR, P)

fwrite(
  gene_conversion_table,
  file.path(OUTPUT_DIR, "gene_id_conversion.txt"),
  sep = "\t",
  quote = FALSE
)

cat("Conversion summary:\n")
cat("  Input genes:", nrow(sig_genes), "\n")
cat("  Converted to symbols:", length(gene_symbols), "\n")
cat(
  "  Conversion rate:",
  round(length(gene_symbols) / nrow(sig_genes) * 100, 1),
  "%\n\n"
)

if (length(gene_symbols) < 10) {
  cat("⚠️  Warning: Very few genes converted!\n")
  cat("Proceeding anyway, but results may be limited.\n\n")
}

# -------------------------------------------------------
# Run EnrichR
# -------------------------------------------------------
cat("Running EnrichR pathway analysis...\n")

dbs <- c(
  "KEGG_2021_Human",
  "GO_Biological_Process_2023",
  "GO_Molecular_Function_2023",
  "Reactome_2022",
  "WikiPathway_2023_Human",
  "MSigDB_Hallmark_2020"
)

cat("Databases:", paste(dbs, collapse = ", "), "\n")
cat("Input genes:", length(gene_symbols), "\n\n")

enrichr_results <- enrichr(gene_symbols, dbs)

cat("✓ EnrichR complete\n\n")

# -------------------------------------------------------
# Process results
# -------------------------------------------------------
cat("Processing results...\n")

all_enrichments <- list()

for (db_name in names(enrichr_results)) {
  result <- enrichr_results[[db_name]]

  if (nrow(result) > 0) {
    result$Database <- db_name
    all_enrichments[[db_name]] <- result
  }
}

enrichment_df <- bind_rows(all_enrichments)

enrichment_df <- enrichment_df %>%
  mutate(
    log_pval = -log10(P.value),
    log_adj_pval = -log10(Adjusted.P.value),
    n_genes = sapply(strsplit(Genes, ";"), length),
    enrichment_score = log_adj_pval * sqrt(n_genes)
  ) %>%
  arrange(Adjusted.P.value)

fwrite(
  enrichment_df,
  file.path(OUTPUT_DIR, "pathway_enrichment_results_all.txt"),
  sep = "\t",
  quote = FALSE
)

enrichment_sig <- enrichment_df %>%
  filter(Adjusted.P.value < 0.1)

fwrite(
  enrichment_sig,
  file.path(OUTPUT_DIR, "pathway_enrichment_results.txt"),
  sep = "\t",
  quote = FALSE
)

cat("✓ Results processed\n")
cat("  Total pathways tested:", nrow(enrichment_df), "\n")
cat("  Significant pathways (FDR<0.1):", nrow(enrichment_sig), "\n\n")

# -------------------------------------------------------
# CREATE VISUALIZATIONS
# -------------------------------------------------------

cat("Creating visualizations...\n\n")

# Database colors
db_colors <- c(
  "KEGG_2021_Human" = "#e74c3c",
  "GO_Biological_Process_2023" = "#3498db",
  "GO_Molecular_Function_2023" = "#2ecc71",
  "Reactome_2022" = "#f39c12",
  "WikiPathway_2023_Human" = "#9b59b6",
  "MSigDB_Hallmark_2020" = "#1abc9c"
)

# -------------------------------------------------------
# PLOT 1: Top Enriched Pathways (Bar Plot)
# -------------------------------------------------------
cat("Creating Plot 1: Top enriched pathways...\n")

top_pathways <- enrichment_df %>%
  filter(Adjusted.P.value < 0.1) %>%
  group_by(Database) %>%
  slice_min(order_by = Adjusted.P.value, n = 5) %>%
  ungroup() %>%
  arrange(desc(log_adj_pval)) %>%
  head(30)

if(nrow(top_pathways) > 0) {
  top_pathways <- top_pathways %>%
    mutate(
      Term = as.character(Term),  # Convert to character first
      Term_short = ifelse(nchar(Term) > 60, 
                         paste0(substr(Term, 1, 57), "..."), 
                         Term),
      Term_short = factor(Term_short, levels = unique(Term_short))
    )
  
  p1 <- ggplot(top_pathways, aes(x = log_adj_pval, y = fct_reorder(Term_short, log_adj_pval), 
                                  fill = Database)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = n_genes), hjust = -0.2, size = 2.5) +
    scale_fill_manual(values = db_colors, name = "Database") +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 8),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    ) +
    labs(
      title = "Top Enriched Pathways",
      x = "-log10(Adjusted P-value)",
      y = "Pathway",
      subtitle = paste("Showing top 30 pathways (FDR < 0.1) |", length(gene_symbols), "input genes")
    ) +
    xlim(0, max(top_pathways$log_adj_pval) * 1.15)
  
  ggsave(
    file.path(OUTPUT_DIR, "enrichment_top_pathways.pdf"),
    p1,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  ggsave(
    file.path(OUTPUT_DIR, "enrichment_top_pathways.png"),
    p1,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  cat("  ✓ Saved: enrichment_top_pathways.pdf/png\n")
} else {
  cat("  ⚠ No significant pathways to plot\n")
}

# -------------------------------------------------------
# PLOT 2: Database Comparison (Dot Plot)
# -------------------------------------------------------
cat("Creating Plot 2: Database comparison...\n")

top_per_db <- enrichment_df %>%
  filter(Adjusted.P.value < 0.1) %>%
  group_by(Database) %>%
  slice_min(order_by = Adjusted.P.value, n = 10) %>%
  ungroup()

if(nrow(top_per_db) > 0) {
  top_per_db <- top_per_db %>%
    mutate(
      Term = as.character(Term),  # Convert to character
      Term_short = ifelse(nchar(Term) > 50, 
                         paste0(substr(Term, 1, 47), "..."), 
                         Term)
    )
  
  p2 <- ggplot(top_per_db, aes(x = log_adj_pval, y = fct_reorder(Term_short, log_adj_pval),
                                size = n_genes, color = Database)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = db_colors, name = "Database") +
    scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 7),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    ) +
    labs(
      title = "Pathway Enrichment Across Databases",
      x = "-log10(Adjusted P-value)",
      y = "Pathway",
      subtitle = "Top 10 pathways per database (FDR < 0.1)"
    )
  
  ggsave(
    file.path(OUTPUT_DIR, "enrichment_by_database.pdf"),
    p2,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  ggsave(
    file.path(OUTPUT_DIR, "enrichment_by_database.png"),
    p2,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  cat("  ✓ Saved: enrichment_by_database.pdf/png\n")
}

# -------------------------------------------------------
# PLOT 3: Enrichment Bubble Plot
# -------------------------------------------------------
cat("Creating Plot 3: Enrichment bubble plot...\n")

bubble_data <- enrichment_df %>%
  filter(Adjusted.P.value < 0.05) %>%
  group_by(Database) %>%
  slice_min(order_by = Adjusted.P.value, n = 8) %>%
  ungroup() %>%
  mutate(
    Term = as.character(Term),  # Convert to character
    Term_short = ifelse(nchar(Term) > 45, 
                       paste0(substr(Term, 1, 42), "..."), 
                       Term)
  )

if(nrow(bubble_data) > 0) {
  p3 <- ggplot(bubble_data, aes(x = Database, y = fct_reorder(Term_short, -log_adj_pval),
                                 size = n_genes, color = log_adj_pval)) +
    geom_point(alpha = 0.7) +
    scale_color_gradient(low = "#fee5d9", high = "#a50f15", name = "-log10(FDR)") +
    scale_size_continuous(name = "Gene Count", range = c(3, 10)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      panel.grid.major = element_line(color = "grey90")
    ) +
    labs(
      title = "Pathway Enrichment Bubble Plot",
      x = "Database",
      y = "Pathway",
      subtitle = "Top 8 pathways per database (FDR < 0.05)"
    )
  
  ggsave(
    file.path(OUTPUT_DIR, "enrichment_bubble_plot.pdf"),
    p3,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  ggsave(
    file.path(OUTPUT_DIR, "enrichment_bubble_plot.png"),
    p3,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  cat("  ✓ Saved: enrichment_bubble_plot.pdf/png\n")
}

# -------------------------------------------------------
# PLOT 4: Summary Statistics per Database
# -------------------------------------------------------
cat("Creating Plot 4: Database summary...\n")

db_summary <- enrichment_df %>%
  group_by(Database) %>%
  summarise(
    Total = n(),
    Significant = sum(Adjusted.P.value < 0.1),
    Mean_log_pval = mean(log_adj_pval[Adjusted.P.value < 0.1], na.rm = TRUE),
    .groups = "drop"
  )

p4a <- ggplot(db_summary, aes(x = fct_reorder(Database, Significant), y = Significant, fill = Database)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Significant), hjust = -0.2, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = db_colors) +
  coord_flip() +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  ) +
  labs(
    title = "Significant Pathways per Database",
    x = "Database",
    y = "Number of Significant Pathways (FDR < 0.1)"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

p4b <- ggplot(db_summary %>% filter(!is.na(Mean_log_pval)), 
              aes(x = fct_reorder(Database, Mean_log_pval), y = Mean_log_pval, fill = Database)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = round(Mean_log_pval, 2)), hjust = -0.2, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = db_colors) +
  coord_flip() +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  ) +
  labs(
    title = "Mean Enrichment Strength",
    x = "Database",
    y = "Mean -log10(Adjusted P-value)"
  )

p4 <- plot_grid(p4a, p4b, ncol = 1, align = "v")

ggsave(
  file.path(OUTPUT_DIR, "enrichment_database_summary.pdf"),
  p4,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(OUTPUT_DIR, "enrichment_database_summary.png"),
  p4,
  width = 10,
  height = 8,
  dpi = 300
)

cat("  ✓ Saved: enrichment_database_summary.pdf/png\n")

# -------------------------------------------------------
# PLOT 5: Gene Ratio vs P-value
# -------------------------------------------------------
cat("Creating Plot 5: Gene ratio plot...\n")

ratio_data <- enrichment_df %>%
  filter(Adjusted.P.value < 0.1) %>%
  mutate(
    gene_ratio = n_genes / sapply(strsplit(Overlap, "/"), function(x) as.numeric(x[2]))
  ) %>%
  group_by(Database) %>%
  slice_min(order_by = Adjusted.P.value, n = 15) %>%
  ungroup()

if(nrow(ratio_data) > 0) {
  p5 <- ggplot(ratio_data, aes(x = gene_ratio, y = log_adj_pval, 
                                color = Database, size = n_genes)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = db_colors, name = "Database") +
    scale_size_continuous(name = "Gene Count", range = c(2, 10)) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    ) +
    labs(
      title = "Gene Ratio vs Significance",
      x = "Gene Ratio (Enriched Genes / Total Genes in Pathway)",
      y = "-log10(Adjusted P-value)",
      subtitle = "Top 15 pathways per database"
    )
  
  ggsave(
    file.path(OUTPUT_DIR, "enrichment_gene_ratio.pdf"),
    p5,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  ggsave(
    file.path(OUTPUT_DIR, "enrichment_gene_ratio.png"),
    p5,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  cat("  ✓ Saved: enrichment_gene_ratio.pdf/png\n")
}

# -------------------------------------------------------
# Final Summary
# -------------------------------------------------------

cat("\n================================\n")
cat("Pathway Analysis Complete!\n")
cat("================================\n\n")

cat("Summary Statistics:\n")
cat("  Total genes analyzed:", length(gene_symbols), "\n")
cat("  Total pathways tested:", nrow(enrichment_df), "\n")
cat("  Significant pathways (FDR<0.1):", nrow(enrichment_sig), "\n")
cat("  Significant pathways (FDR<0.05):", sum(enrichment_df$Adjusted.P.value < 0.05), "\n\n")

cat("Output files saved to:", OUTPUT_DIR, "\n")
cat("  - gene_id_conversion.txt\n")
cat("  - pathway_enrichment_results_all.txt\n")
cat("  - pathway_enrichment_results.txt\n")
cat("  - enrichment_top_pathways.pdf/png\n")
cat("  - enrichment_by_database.pdf/png\n")
cat("  - enrichment_bubble_plot.pdf/png\n")
cat("  - enrichment_database_summary.pdf/png\n")
cat("  - enrichment_gene_ratio.pdf/png\n\n")

cat("✓ Analysis complete!\n")