#!/bin/bash

################################################################################
# Case-Case GWAS Analysis Pipeline
# Compares infection by two MTBC lineages using logistic regression
#
# USAGE:
#   bash scripts/03_case_case_gwas.sh
#
# REQUIREMENTS:
#   - PLINK v1.9 or v2.0 (v2.0 recommended for Firth correction)
#   - R 4.0+
#   - R packages: data.table, dplyr
#   - Input files: genotypes (PLINK format), phenotype, covariates (PCs, age, sex)
#
# INPUT FILES EXPECTED:
#   - data/genotypes.bed/bim/fam (imputed PLINK format)
#   - data/samples_case_case.txt (list of samples: FID IID)
#   - data/phenotype.txt (FID IID PHENO)
#   - data/covariates.txt (FID IID PC1-PC10 AGE SEX)
#
# OUTPUT FILES:
#   - results/gwas_case_case_cleaned.txt
#   - results/gwas_top_hits.txt
#   - results/gwas_qc_summary.txt
#
# AUTHOR: Theophilus Afum
# DATE: 2024
# LICENSE: MIT
################################################################################

set -euo pipefail

# Enable error reporting
trap 'echo "Error on line $LINENO"; exit 1' ERR

################################################################################
# CONFIGURATION - EDIT THESE
################################################################################

# Input/Output directories
DATA_DIR="data"
RESULTS_DIR="results"
LOGS_DIR="logs"

# Input files
GENO_PREFIX="${DATA_DIR}/genotypes"
SAMPLES_FILE="${DATA_DIR}/samples_case_case.txt"
PHENOTYPE_FILE="${DATA_DIR}/phenotype.txt"
COVARIATES_FILE="${DATA_DIR}/covariates.txt"

# Output prefix
OUTPUT_PREFIX="${RESULTS_DIR}/gwas_case_case"

# GWAS parameters
MAF_THRESHOLD=0.01
HWE_THRESHOLD=1e-6
GENO_MISSING=0.05
GWS_THRESHOLD=5e-8
SUGGESTIVE_THRESHOLD=1e-5
NUM_PCs=10

################################################################################
# SETUP
################################################################################

# Create output directories
mkdir -p "${RESULTS_DIR}" "${LOGS_DIR}"

echo "================================================================================"
echo "Case-Case GWAS: Logistic Regression Analysis"
echo "================================================================================"
echo ""
echo "Configuration:"
echo "  Input genotypes: ${GENO_PREFIX}"
echo "  Phenotype file: ${PHENOTYPE_FILE}"
echo "  Covariates file: ${COVARIATES_FILE}"
echo "  Output directory: ${RESULTS_DIR}"
echo "  MAF threshold: ${MAF_THRESHOLD}"
echo "  HWE threshold: ${HWE_THRESHOLD}"
echo "================================================================================"
echo ""

# Check input files exist
for file in "${PHENOTYPE_FILE}" "${COVARIATES_FILE}" "${SAMPLES_FILE}"; do
    if [ ! -f "${file}" ]; then
        echo "ERROR: Input file not found: ${file}"
        exit 1
    fi
done

# Check genotype files
if [ ! -f "${GENO_PREFIX}.bed" ] || [ ! -f "${GENO_PREFIX}.bim" ] || [ ! -f "${GENO_PREFIX}.fam" ]; then
    echo "ERROR: Genotype files not found (${GENO_PREFIX}.{bed,bim,fam})"
    exit 1
fi

################################################################################
# STEP 1: Extract case-case samples
################################################################################

echo "Step 1: Extracting case-case samples..."

plink --bfile "${GENO_PREFIX}" \
      --keep "${SAMPLES_FILE}" \
      --make-bed \
      --out "${RESULTS_DIR}/genotypes_cases_only" \
      2>&1 | tee "${LOGS_DIR}/step1_sample_extraction.log"

N_SAMPLES=$(wc -l < "${SAMPLES_FILE}")
echo "✓ Extracted ${N_SAMPLES} case-case samples"
echo ""

################################################################################
# STEP 2: Quality Control
################################################################################

echo "Step 2: Quality control on case-case samples..."

plink --bfile "${RESULTS_DIR}/genotypes_cases_only" \
      --geno "${GENO_MISSING}" \
      --maf "${MAF_THRESHOLD}" \
      --hwe "${HWE_THRESHOLD}" \
      --make-bed \
      --out "${RESULTS_DIR}/genotypes_qc" \
      2>&1 | tee "${LOGS_DIR}/step2_qc.log"

echo "✓ QC complete"
echo ""

################################################################################
# STEP 3: Principal Component Analysis (PCA)
################################################################################

echo "Step 3: Computing principal components..."

# LD pruning for PCA
plink --bfile "${RESULTS_DIR}/genotypes_qc" \
      --indep-pairwise 50 5 0.2 \
      --out "${RESULTS_DIR}/ld_pruned" \
      2>&1 | tee "${LOGS_DIR}/step3a_ld_pruning.log"

# Compute PCs
plink --bfile "${RESULTS_DIR}/genotypes_qc" \
      --extract "${RESULTS_DIR}/ld_pruned.prune.in" \
      --pca "${NUM_PCs}" \
      --out "${RESULTS_DIR}/pca_results" \
      2>&1 | tee "${LOGS_DIR}/step3b_pca.log"

# Format PCs for use as covariates
awk 'BEGIN {
    printf "FID\tIID"
    for (i = 1; i <= 10; i++) printf "\tPC%d", i
    printf "\n"
}
NR > 1 {
    printf "%s\t%s", $1, $2
    for (i = 3; i <= 12; i++) printf "\t%s", $i
    printf "\n"
}' "${RESULTS_DIR}/pca_results.eigenvec" > "${RESULTS_DIR}/pca_formatted.txt"

echo "✓ PCA complete (${NUM_PCs} components)"
echo ""

################################################################################
# STEP 4: Prepare Final Covariates
################################################################################

echo "Step 4: Merging covariates (PCs, age, sex)..."

Rscript --vanilla - <<'REND'
library(data.table)
library(dplyr)

# Detect script directory for relative imports
script_dir <- dirname(sys.frame(1)$ofile)
if (is.null(script_dir)) script_dir <- getwd()

# Load input files
message("Loading phenotype and covariate files...")
phenotype <- fread("data/phenotype.txt")
covariates_input <- fread("data/covariates.txt")

message("Phenotype samples: ", nrow(phenotype))
message("Covariate samples: ", nrow(covariates_input))

# Ensure FID and IID exist
if (!all(c("FID", "IID") %in% names(phenotype))) {
    stop("Phenotype file must have columns: FID IID PHENO")
}
if (!all(c("FID", "IID") %in% names(covariates_input))) {
    stop("Covariates file must have columns: FID IID [covariates]")
}

# Merge phenotype and covariates
final_covariates <- covariates_input %>%
    left_join(phenotype %>% select(FID, IID), by = c("FID", "IID")) %>%
    select(FID, IID, everything())

# Check for missing values
missing_fid <- sum(is.na(final_covariates$FID))
missing_iid <- sum(is.na(final_covariates$IID))

if (missing_fid > 0 || missing_iid > 0) {
    warning("Found missing FID/IID values")
}

# Write merged covariates
fwrite(final_covariates, "results/covariates_final.txt",
       sep = "\t", quote = FALSE)

message("\n=== Covariate Summary ===")
message("Total samples: ", nrow(final_covariates))
message("Columns: ", paste(names(final_covariates), collapse = ", "))
message("✓ Covariates ready for GWAS")
REND

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to prepare covariates"
    exit 1
fi

echo "✓ Covariates prepared"
echo ""

################################################################################
# STEP 5: Run GWAS
################################################################################

echo "Step 5: Running case-case GWAS with logistic regression..."
echo ""

# Detect PLINK version
if command -v plink2 &> /dev/null; then
    PLINK_VERSION="plink2"
    echo "Using PLINK2 with Firth fallback correction..."
    echo ""
    
    plink2 --bfile "${RESULTS_DIR}/genotypes_qc" \
           --pheno "${PHENOTYPE_FILE}" \
           --covar "${RESULTS_DIR}/covariates_final.txt" \
           --covar-variance-standardize \
           --glm firth-fallback hide-covar \
           --ci 0.95 \
           --out "${OUTPUT_PREFIX}" \
           2>&1 | tee "${LOGS_DIR}/step5_gwas.log"
    
    RESULTS_FILE="${OUTPUT_PREFIX}.*.glm.firth"
    
elif plink --version 2>&1 | grep -q "v2"; then
    PLINK_VERSION="plink"
    echo "Using PLINK v2 with Firth fallback correction..."
    echo ""
    
    plink --bfile "${RESULTS_DIR}/genotypes_qc" \
          --pheno "${PHENOTYPE_FILE}" \
          --covar "${RESULTS_DIR}/covariates_final.txt" \
          --covar-variance-standardize \
          --glm firth-fallback hide-covar \
          --ci 0.95 \
          --out "${OUTPUT_PREFIX}" \
          2>&1 | tee "${LOGS_DIR}/step5_gwas.log"
    
    RESULTS_FILE="${OUTPUT_PREFIX}.*.glm.firth"
    
else
    PLINK_VERSION="plink1.9"
    echo "⚠️  Warning: Using PLINK v1.9 (Firth correction unavailable)"
    echo "Recommendation: Install PLINK2 for improved handling of rare variants"
    echo ""
    
    plink --bfile "${RESULTS_DIR}/genotypes_qc" \
          --pheno "${PHENOTYPE_FILE}" \
          --covar "${RESULTS_DIR}/covariates_final.txt" \
          --logistic hide-covar \
          --ci 0.95 \
          --out "${OUTPUT_PREFIX}" \
          2>&1 | tee "${LOGS_DIR}/step5_gwas.log"
    
    RESULTS_FILE="${OUTPUT_PREFIX}.assoc.logistic"
fi

echo "✓ GWAS complete (${PLINK_VERSION})"
echo ""

################################################################################
# STEP 6: Process and Summarize Results
################################################################################

echo "Step 6: Processing GWAS results..."

Rscript --vanilla - <<'REND'
library(data.table)
library(dplyr)

# Find results file
results_files <- list.files("results", pattern = "gwas_case_case.*\\.(glm|assoc)", full.names = TRUE)

if (length(results_files) == 0) {
    stop("No GWAS results files found in results/ directory")
}

main_file <- results_files[1]
message("Loading results from: ", main_file)

# Load results
gwas <- tryCatch({
    fread(main_file)
}, error = function(e) {
    stop("Failed to read GWAS results: ", e$message)
})

message("Loaded ", nrow(gwas), " variants")

# Filter results
if ("TEST" %in% names(gwas)) {
    gwas_clean <- gwas %>% filter(TEST == "ADD" | is.na(TEST))
    message("Filtered to ADD test variants: ", nrow(gwas_clean))
} else {
    gwas_clean <- gwas
}

# Remove invalid p-values
gwas_clean <- gwas_clean %>%
    filter(!is.na(P) & P > 0 & P <= 1)

# Calculate genomic inflation
observed_chisq <- qchisq(gwas_clean$P, df = 1, lower.tail = FALSE)
observed_chisq <- observed_chisq[!is.na(observed_chisq)]
lambda <- median(observed_chisq, na.rm = TRUE) / qchisq(0.5, df = 1)

# Count significant variants
n_gws <- sum(gwas_clean$P < 5e-8, na.rm = TRUE)
n_suggestive <- sum(gwas_clean$P < 1e-5, na.rm = TRUE)

message("\n=== GWAS Summary ===")
message("Total variants tested: ", nrow(gwas_clean))
message("Genomic inflation (λ): ", round(lambda, 4))
message("Genome-wide significant (P < 5e-8): ", n_gws)
message("Suggestive (P < 1e-5): ", n_suggestive)

# Save cleaned results
fwrite(gwas_clean, "results/gwas_case_case_cleaned.txt",
       sep = "\t", quote = FALSE)

# Save summary
summary_stats <- data.frame(
    metric = c("total_variants", "lambda", "gws_count", "suggestive_count"),
    value = c(nrow(gwas_clean), round(lambda, 4), n_gws, n_suggestive)
)
fwrite(summary_stats, "results/gwas_qc_summary.txt",
       sep = "\t", quote = FALSE)

# Save top hits if any
if (n_suggestive > 0) {
    top_hits <- gwas_clean %>%
        filter(P < 1e-5) %>%
        arrange(P) %>%
        select(any_of(c("ID", "SNP", "CHR", "BP", "A1", "A2", "OR", "SE", "P", "BETA")))
    
    fwrite(top_hits, "results/gwas_top_hits.txt",
           sep = "\t", quote = FALSE)
    
    message("\n=== Top 10 Hits ===")
    print(head(top_hits, 10))
} else {
    message("\nNo variants with P < 1e-5 found.")
    file.create("results/gwas_top_hits.txt")  # Create empty file
}

message("\n✓ Results processed and saved")
REND

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to process GWAS results"
    exit 1
fi

echo "✓ Results processed"
echo ""

################################################################################
# FINAL SUMMARY
################################################################################

echo "================================================================================"
echo "GWAS Analysis Complete!"
echo "================================================================================"
echo ""
echo "Output files:"
echo "  ✓ results/gwas_case_case_cleaned.txt"
echo "  ✓ results/gwas_top_hits.txt"
echo "  ✓ results/gwas_qc_summary.txt"
echo ""
echo "Intermediate files (for reference):"
echo "  - results/genotypes_qc.bed/bim/fam"
echo "  - results/pca_results.eigenvec"
echo "  - results/covariates_final.txt"
echo ""
echo "Log files:"
echo "  - logs/step*.log"
echo ""
echo "Next steps:"
echo "  1. Review gwas_qc_summary.txt for genomic inflation (λ)"
echo "  2. Run visualization script: bash scripts/06_visualization.sh"
echo "  3. Perform gene-based analysis: bash scripts/04_magma_analysis.sh"
echo "================================================================================"
echo ""
