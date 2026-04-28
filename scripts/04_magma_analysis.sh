#!/bin/bash

################################################################################
# MAGMA Genome-Wide Gene-Based Analysis
# Analyzes ALL chromosomes for unbiased discovery
# Nature publication quality
################################################################################

set -e

# -------------------------------------------------------
# Configuration
# -------------------------------------------------------
GWAS_FILE="/home/tafum/H3_imputation/final/GWAS/Final_GWAS/GWAS_MAF_vs_Mtbss_cleaned.txt"
GENO_PREFIX="/home/tafum/H3_imputation/final/GWAS/Final_GWAS/imputed_all_chr_QC"
OUTPUT_DIR="/home/tafum/H3_imputation/final/GWAS/Final_GWAS/GWAS_plots/MAGMA_analysis"
GENE_LOC_FILE="$OUTPUT_DIR/NCBI38.gene.loc"
GENESETS_FILE="/home/tafum/H3_imputation/final/GWAS/Final_GWAS/GWAS_plots/tb_immune_genesets.txt"

echo "================================"
echo "MAGMA Genome-Wide Analysis"
echo "ALL Chromosomes"
echo "================================"
echo ""

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# -------------------------------------------------------
# Step 1: Prepare genome-wide GWAS data
# -------------------------------------------------------
echo "Step 1: Preparing genome-wide GWAS data..."

Rscript - <<'EOF'
library(data.table)
library(dplyr)

gwas <- fread("/home/tafum/H3_imputation/final/GWAS/Final_GWAS/GWAS_MAF_vs_Mtbss_cleaned.txt")

# Standardize columns
if ("#CHROM" %in% names(gwas)) {
  gwas <- gwas %>% rename(CHR = `#CHROM`, BP = POS, SNP = ID)
}

# Genome-wide file
gwas_full <- gwas %>%
  filter(CHR >= 1 & CHR <= 22) %>%  # Autosomal only
  select(SNP, CHR, BP, P, N = OBS_CT)

fwrite(gwas_full, "gwas_genomewide_magma_input.txt", sep="\t", quote=FALSE)

cat("? Prepared", format(nrow(gwas_full), big.mark=","), "variants\n")
cat("  Chromosomes:", paste(sort(unique(gwas_full$CHR)), collapse=", "), "\n")
EOF

echo ""

# -------------------------------------------------------
# Step 2: Annotate SNPs to genes (genome-wide)
# -------------------------------------------------------
echo "Step 2: Annotating SNPs to genes (genome-wide)..."

magma \
  --annotate window=10,10 \
  --snp-loc gwas_genomewide_magma_input.txt \
  --gene-loc $GENE_LOC_FILE \
  --out genomewide_snp_gene_annotation

echo "? Genome-wide SNP annotation complete"
echo ""

# -------------------------------------------------------
# Step 3: Gene-based analysis (genome-wide)
# -------------------------------------------------------
echo "Step 3: Running genome-wide gene-based analysis..."

if [ -f "${GENO_PREFIX}.bed" ]; then
    echo "Using PLINK files for LD-aware analysis"
    magma \
      --bfile ${GENO_PREFIX} \
      --pval gwas_genomewide_magma_input.txt ncol=N \
      --gene-annot genomewide_snp_gene_annotation.genes.annot \
      --out genome_wide_gene_analysis
else
    echo "Running summary-statistics-only analysis"
    magma \
      --pval gwas_genomewide_magma_input.txt ncol=N \
      --gene-annot genomewide_snp_gene_annotation.genes.annot \
      --out genome_wide_gene_analysis
fi

echo "? Genome-wide gene analysis complete"
echo ""

# -------------------------------------------------------
# Step 4: Gene-set analysis
# -------------------------------------------------------
echo "Step 4: Running gene-set analysis..."

if [ -f "${GENESETS_FILE}" ]; then
    magma \
      --gene-results genome_wide_gene_analysis.genes.raw \
      --set-annot ${GENESETS_FILE} \
      --out genome_wide_geneset_analysis
    
    echo "? Gene-set analysis complete"
else
    echo "??  Gene-set file not found: ${GENESETS_FILE}"
    echo "  Skipping gene-set analysis"
fi

echo ""

# -------------------------------------------------------
# Step 5: Summary statistics
# -------------------------------------------------------
echo "Step 5: Generating summary..."

Rscript - <<'EOF'
library(data.table)
library(dplyr)

genes <- fread("genome_wide_gene_analysis.genes.out")

cat("\n=== GENOME-WIDE RESULTS ===\n")
cat("Total genes tested:", format(nrow(genes), big.mark=","), "\n")
cat("Bonferroni threshold: P <", format(0.05/nrow(genes), scientific=TRUE), "\n")
cat("Bonferroni significant:", sum(genes$P < 0.05/nrow(genes)), "\n")
cat("Nominally significant (P<0.05):", sum(genes$P < 0.05), "\n")
cat("Suggestive (P<0.001):", sum(genes$P < 0.001), "\n\n")

# Top genes
cat("Top 20 genes:\n")
top <- genes %>%
  arrange(P) %>%
  select(GENE, CHR, START, STOP, NSNPS, ZSTAT, P) %>%
  head(20)
print(top)

# By chromosome
cat("\n\nSignificant genes by chromosome:\n")
chr_summary <- genes %>%
  filter(P < 0.05) %>%
  group_by(CHR) %>%
  summarise(
    n_genes = n(),
    top_gene = GENE[which.min(P)],
    top_P = min(P),
    .groups = "drop"
  ) %>%
  arrange(CHR)

if (nrow(chr_summary) > 0) {
  print(chr_summary)
} else {
  cat("  No significant genes\n")
}

# Save summary
fwrite(top, "top_20_genes_genomewide.txt", sep="\t", quote=FALSE)

# -------------------------------------------------------
# Gene-set results (if available)
# -------------------------------------------------------
cat("\n\n=== GENE-SET RESULTS ===\n")

if (file.exists("genome_wide_geneset_analysis.gsa.out")) {
  genesets <- fread("genome_wide_geneset_analysis.gsa.out")
  
  cat("Gene sets tested:", nrow(genesets), "\n")
  cat("Significant (P<0.05):", sum(genesets$P < 0.05), "\n\n")
  
  if (nrow(genesets) > 0) {
    cat("Top 10 gene sets:\n")
    top_sets <- genesets %>%
      arrange(P) %>%
      select(VARIABLE, NGENES, BETA, BETA_STD, SE, P) %>%
      head(10)
    print(top_sets)
    
    # Save
    fwrite(genesets, "genome_wide_geneset_results.txt", sep="\t", quote=FALSE)
  }
} else {
  cat("Gene-set results file not found (this may be expected)\n")
}
EOF

echo ""

# -------------------------------------------------------
# Summary
# -------------------------------------------------------
echo "================================"
echo "Genome-Wide MAGMA Complete!"
echo "================================"
echo ""
echo "Output files:"
echo "  - genome_wide_gene_analysis.genes.out"
echo "  - genome_wide_gene_analysis.genes.raw"
echo "  - genome_wide_geneset_analysis.gsa.out"
echo "  - top_20_genes_genomewide.txt"
echo ""
echo "Next steps:"
echo "  1. Rscript visualize_magma_genomewide.R"
echo "  2. Rscript run_pathway_enrichment.R"
echo "  3. Rscript create_integrated_figure.R"
echo ""
echo "================================"

cd ..