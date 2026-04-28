#!/bin/bash

################################################################################
# ADMIXTURE Analysis Pipeline
# Runs ADMIXTURE for K=2 to K=10 and prepares data for R plotting
################################################################################

set -e  # Exit on error

# -------------------------------------------------------
# Configuration
# -------------------------------------------------------
INPUT_PREFIX="H3_QC8_ibdClean"
OUTPUT_DIR="ADMIXTURE_results"
MIN_K=2
MAX_K=10
THREADS=4

echo "=================================="
echo "ADMIXTURE Analysis Pipeline"
echo "=================================="
echo "Input files: ${INPUT_PREFIX}.bed/bim/fam"
echo "K range: ${MIN_K} to ${MAX_K}"
echo "Output directory: ${OUTPUT_DIR}"
echo "=================================="
echo ""

# -------------------------------------------------------
# Create output directory
# -------------------------------------------------------
mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# Copy input files to output directory
echo "Copying input files..."
cp ../${INPUT_PREFIX}.bed .
cp ../${INPUT_PREFIX}.bim .
cp ../${INPUT_PREFIX}.fam .

# -------------------------------------------------------
# Run ADMIXTURE for each K
# -------------------------------------------------------
echo ""
echo "Running ADMIXTURE analysis..."
echo ""

for K in $(seq ${MIN_K} ${MAX_K}); do
    echo "========================================"
    echo "Running ADMIXTURE for K=${K}"
    echo "========================================"
    
    # Run ADMIXTURE with cross-validation
    admixture --cv ${INPUT_PREFIX}.bed ${K} -j${THREADS} | tee log${K}.out
    
    # Rename output files to include K value
    mv ${INPUT_PREFIX}.${K}.Q ${INPUT_PREFIX}.K${K}.Q
    mv ${INPUT_PREFIX}.${K}.P ${INPUT_PREFIX}.K${K}.P
    
    echo "✓ Completed K=${K}"
    echo ""
done

# -------------------------------------------------------
# Extract cross-validation errors
# -------------------------------------------------------
echo ""
echo "Extracting cross-validation errors..."

grep "CV error" log*.out > cv_errors.txt

echo ""
echo "Cross-validation errors:"
cat cv_errors.txt

echo ""
echo "=================================="
echo "ADMIXTURE analysis completed!"
echo "Results saved in: ${OUTPUT_DIR}/"
echo "=================================="
echo ""
echo "Next step: Run the R plotting script"
echo "  Rscript plot_admixture.R"
echo ""

cd ..
