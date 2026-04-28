#!/bin/bash

################################################################################
# Case-Case GWAS: MAF vs Mtbss - CORRECTED VERSION
# Fixed phenotype coding for PLINK2
################################################################################

set -e

# -------------------------------------------------------
# Configuration
# -------------------------------------------------------
GENO_PREFIX="imputed_rsID"
PHENO_FILE="H3_pheno.phe.txt"
OUTPUT_PREFIX="GWAS_MAF_vs_Mtbss"

echo "================================"
echo "Case-Case GWAS: MAF vs Mtbss"
echo "CORRECTED VERSION"
echo "================================"
echo ""

# -------------------------------------------------------
# Step 1: Fix phenotype coding
# -------------------------------------------------------
# echo "Step 1: Preparing phenotype with correct coding..."

# Rscript fix_phenotype.R

# if [ ! -f "phenotype_MAF_vs_Mtbss_fixed.txt" ]; then
 #    echo "Error: Failed to create phenotype file!"
#    exit 1
# fi

# echo ""

# -------------------------------------------------------
# Step 2: Extract cases from genotype data
# -------------------------------------------------------
echo "Step 2: Extracting MAF and Mtbss cases..."

plink --bfile ${GENO_PREFIX} \
      --keep samples_MAF_Mtbss.txt \
      --make-bed \
      --out ${GENO_PREFIX}_cases_only

echo "✓ Extracted $(wc -l < samples_MAF_Mtbss.txt) cases"
echo ""

# -------------------------------------------------------
# Step 3: QC on cases
# -------------------------------------------------------
echo "Step 3: Quality control..."

plink --bfile ${GENO_PREFIX}_cases_only \
      --geno 0.05 \
      --maf 0.01 \
      --hwe 1e-6 \
      --make-bed \
      --out ${GENO_PREFIX}_cases_QC

echo "✓ QC complete"
echo ""

# -------------------------------------------------------
# Step 4: Check if we already have PCs
# -------------------------------------------------------
echo "Step 4: Preparing principal components..."

if [ ! -f "covariates_PCs.txt" ]; then
    echo "Computing PCA on cases..."
    
    # LD pruning
    plink --bfile ${GENO_PREFIX}_cases_QC \
          --indep-pairwise 50 5 0.2 \
          --out ld_pruned
    
    # PCA
    plink --bfile ${GENO_PREFIX}_cases_QC \
          --extract ld_pruned.prune.in \
          --pca 10 \
          --out pca_cases
    
    # Format for covariates
    awk 'BEGIN {print "FID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10"} \
         {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' \
         pca_cases.eigenvec > covariates_PCs.txt
    
    echo "✓ Computed PCs"
else
    echo "✓ Using existing PCs: covariates_PCs.txt"
fi

echo ""

# -------------------------------------------------------
# Step 5: Prepare covariates
# -------------------------------------------------------
echo "Step 5: Merging covariates..."

Rscript - <<'EOF'
library(data.table)
library(dplyr)

# Load PCs
pcs <- fread("covariates_PCs.txt")

# Load phenotype to get AGE
pheno <- fread("phenotype_MAF_vs_Mtbss_fixed.txt")
pheno_full <- fread("H3_pheno.phe.txt")

# Merge to get AGE
pheno_with_age <- pheno %>%
  left_join(pheno_full %>% select(FID, IID, AGE), by = c("FID", "IID"))

# Get SEX from FAM
fam <- fread("imputed_rsID_cases_QC.fam", header = FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO_FAM")

# Merge all
covariates <- pcs %>%
  left_join(pheno_with_age %>% select(FID, IID, AGE), by = c("FID", "IID")) %>%
  left_join(fam %>% select(FID, IID, SEX), by = c("FID", "IID"))

# Recode
covariates <- covariates %>%
  mutate(
    AGE = ifelse(is.na(AGE) | AGE == "", -9, AGE),
    SEX = ifelse(SEX == 2, 0, ifelse(SEX == 1, 1, -9))  # 1=male, 0=female
  )

# Write
fwrite(covariates, "covariates_final.txt",
       sep = "\t", col.names = TRUE, quote = FALSE)

cat("\n=== Covariate Summary ===\n")
cat("Samples:", nrow(covariates), "\n")
cat("Covariates: PC1-PC10, AGE, SEX\n")
cat("AGE missing:", sum(covariates$AGE == -9), "\n")
cat("SEX missing:", sum(covariates$SEX == -9), "\n")
cat("✓ Covariates ready\n")
EOF

echo ""

# -------------------------------------------------------
# Step 6: Run GWAS with PLINK2
# -------------------------------------------------------
echo "Step 6: Running GWAS..."

# Check which PLINK version we have
if command -v plink2 &> /dev/null; then
    echo "Using PLINK2 with Firth fallback..."
    
    plink2 --bfile ${GENO_PREFIX}_cases_QC \
           --pheno phenotype_MAF_vs_Mtbss_fixed.txt \
           --covar covariates_final.txt \
           --covar-variance-standardize \
           --glm firth-fallback hide-covar \
           --ci 0.95 \
           --out ${OUTPUT_PREFIX}
    
    echo "✓ PLINK2 GWAS complete"
    
elif plink --version 2>&1 | grep -q "PLINK v2"; then
    echo "Using PLINK v2 with Firth fallback..."
    
    plink --bfile ${GENO_PREFIX}_cases_QC \
          --pheno phenotype_MAF_vs_Mtbss_fixed.txt \
          --covar covariates_final.txt \
          --covar-variance-standardize \
          --glm firth-fallback hide-covar \
          --ci 0.95 \
          --out ${OUTPUT_PREFIX}
    
    echo "✓ PLINK v2 GWAS complete"
    
else
    echo "Using PLINK 1.9 (no Firth correction available)..."
    echo "⚠️  Warning: Consider installing PLINK2 for Firth regression"
    
    plink --bfile ${GENO_PREFIX}_cases_QC \
          --pheno phenotype_MAF_vs_Mtbss_fixed.txt \
          --covar covariates_final.txt \
          --logistic hide-covar \
          --ci 0.95 \
          --out ${OUTPUT_PREFIX}
    
    echo "✓ PLINK 1.9 logistic regression complete"
fi

echo ""

# -------------------------------------------------------
# Step 7: Check and process results
# -------------------------------------------------------
echo "Step 7: Processing results..."

Rscript - <<'EOF'
library(data.table)
library(dplyr)

# Find the results file
results_files <- list.files(pattern = "GWAS_MAF_vs_Mtbss.*\\.(glm|assoc)")

if (length(results_files) == 0) {
  stop("No GWAS results files found!")
}

cat("Found results files:\n")
for (f in results_files) {
  cat(" ", f, "\n")
}

# Load the main results file
main_file <- results_files[grepl("glm.logistic.hybrid|glm.firth|assoc.logistic", results_files)][1]

if (is.na(main_file)) {
  main_file <- results_files[1]
}

cat("\nProcessing:", main_file, "\n")

gwas <- fread(main_file)

# Check structure
cat("\nColumns in results:\n")
print(names(gwas))

# Filter
if ("TEST" %in% names(gwas)) {
  gwas_clean <- gwas %>% filter(TEST == "ADD" | is.na(TEST))
} else {
  gwas_clean <- gwas
}

gwas_clean <- gwas_clean %>%
  filter(!is.na(P) & P > 0 & P <= 1)

# Calculate lambda
observed_chisq <- qchisq(gwas_clean$P, df = 1, lower.tail = FALSE)
lambda <- median(observed_chisq, na.rm = TRUE) / qchisq(0.5, df = 1)

cat("\n=== GWAS Summary ===\n")
cat("Total variants:", nrow(gwas_clean), "\n")
cat("Genomic inflation (λ):", round(lambda, 4), "\n")
cat("Genome-wide sig (P<5e-8):", sum(gwas_clean$P < 5e-8), "\n")
cat("Suggestive (P<1e-5):", sum(gwas_clean$P < 1e-5), "\n")

# Save cleaned results
fwrite(gwas_clean, "GWAS_MAF_vs_Mtbss_cleaned.txt",
       sep = "\t", quote = FALSE)

# Save top hits
if (sum(gwas_clean$P < 1e-5) > 0) {
  top_hits <- gwas_clean %>%
    filter(P < 1e-5) %>%
    arrange(P) %>%
    head(100)
  
  fwrite(top_hits, "GWAS_MAF_vs_Mtbss_top_hits.txt",
         sep = "\t", quote = FALSE)
  
  cat("\nTop 10 hits:\n")
  print(head(top_hits %>% select(any_of(c("ID", "SNP", "CHR", "BP", "P"))), 10))
}

cat("\n✓ Results processed\n")
EOF

echo ""

# -------------------------------------------------------
# Final summary
# -------------------------------------------------------
echo "================================"
echo "GWAS Complete!"
echo "================================"
echo ""
echo "Output files:"
echo "  - GWAS_MAF_vs_Mtbss_cleaned.txt"
echo "  - GWAS_MAF_vs_Mtbss_top_hits.txt (if P<1e-5)"
echo "  - phenotype_MAF_vs_Mtbss_fixed.txt"
echo "  - covariates_final.txt"
echo ""
echo "Next: Rscript visualize_gwas.R"
echo "================================"