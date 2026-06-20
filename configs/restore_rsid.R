#!/usr/bin/env Rscript

################################################################################
# Restore rsIDs to PLINK files from annotation file
# Maps CHR:POS:REF:ALT to rsIDs
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# -------------------------------------------------------
# Configuration
# -------------------------------------------------------
# Paths can be overridden via environment variables; defaults are relative to repo root.
ANNOTATION_FILE <- Sys.getenv("ANNOTATION_FILE", "data/H3_Africa_order_3220_original_annotation.anno.txt")  # Annotation file
INPUT_PREFIX <- Sys.getenv("INPUT_PREFIX", "data/imputed_all_chr_QC")  # Input PLINK files (without extension)
OUTPUT_PREFIX <- Sys.getenv("OUTPUT_PREFIX", "imputed_rsID") # Output PLINK files (without extension)

cat("================================\n")
cat("Restore rsIDs to PLINK files\n")
cat("================================\n")
cat("Annotation file:", ANNOTATION_FILE, "\n")
cat("Input PLINK:", INPUT_PREFIX, "\n")
cat("Output PLINK:", OUTPUT_PREFIX, "\n")
cat("================================\n\n")

# -------------------------------------------------------
# 1) Load annotation file
# -------------------------------------------------------
cat("=== Loading annotation file ===\n")

# Read annotation file
# Format: rsID rsID CHR POS STRAND CONTEXT REF ALT
annot <- fread(ANNOTATION_FILE, header = FALSE)

# Check if file has header
if (grepl("^rs|^kgp|^snp", annot$V1[1])) {
  # No header, assign column names
  colnames(annot) <- c("rsID_1", "rsID_2", "CHR", "POS", "STRAND", "CONTEXT", "REF", "ALT")
} else {
  # Has header, reload with header
  annot <- fread(ANNOTATION_FILE, header = TRUE)
  colnames(annot) <- c("rsID_1", "rsID_2", "CHR", "POS", "STRAND", "CONTEXT", "REF", "ALT")
}

# Use the first rsID column (they appear to be duplicates)
annot <- annot %>%
  select(rsID = rsID_1, CHR, POS, REF, ALT) %>%
  mutate(CHR = as.character(CHR))

cat("Loaded", nrow(annot), "variants from annotation file.\n")
cat("Chromosomes:", paste(unique(annot$CHR), collapse = ", "), "\n")
cat("First few entries:\n")
print(head(annot, 3))
cat("\n")

# -------------------------------------------------------
# 2) Load PLINK BIM file
# -------------------------------------------------------
cat("=== Loading PLINK BIM file ===\n")

bim_file <- paste0(INPUT_PREFIX, ".bim")
bim <- fread(bim_file, header = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "POS", "A1", "A2")

cat("Loaded", nrow(bim), "variants from BIM file.\n")
cat("Chromosomes:", paste(unique(bim$CHR), collapse = ", "), "\n")
cat("First few entries:\n")
print(head(bim, 3))
cat("\n")

# -------------------------------------------------------
# 3) Create matching keys
# -------------------------------------------------------
cat("=== Creating matching keys ===\n")

# For annotation: CHR:POS:REF:ALT
annot <- annot %>%
  mutate(
    CHR = as.character(CHR),
    match_key = paste(CHR, POS, REF, ALT, sep = ":")
  )

# For BIM: need to extract CHR:POS:REF:ALT from SNP ID
# Common formats after imputation:
# CHR:POS:REF:ALT
# CHR_POS_REF_ALT
# CHR-POS-REF-ALT

bim <- bim %>%
  mutate(
    CHR = as.character(CHR),
    # Try to extract components from SNP ID
    SNP_original = SNP,
    # Create match key using BIM position and alleles
    match_key_1 = paste(CHR, POS, A1, A2, sep = ":"),
    match_key_2 = paste(CHR, POS, A2, A1, sep = ":")  # Try flipped alleles
  )

cat("Created matching keys.\n")
cat("Example BIM match keys:\n")
print(head(bim[, c("SNP", "match_key_1", "match_key_2")], 3))
cat("\n")

# -------------------------------------------------------
# 4) Match and update rsIDs
# -------------------------------------------------------
cat("=== Matching variants ===\n")

# Create lookup table from annotation
rsid_lookup <- setNames(annot$rsID, annot$match_key)

# Try matching with both orientations
bim$rsID_matched <- rsid_lookup[bim$match_key_1]
bim$rsID_matched[is.na(bim$rsID_matched)] <- rsid_lookup[bim$match_key_2[is.na(bim$rsID_matched)]]

# Count matches
n_matched <- sum(!is.na(bim$rsID_matched))
pct_matched <- round(n_matched / nrow(bim) * 100, 2)

cat("Matched", n_matched, "variants (", pct_matched, "%).\n")
cat("Unmatched:", nrow(bim) - n_matched, "variants.\n\n")

# -------------------------------------------------------
# 5) Update BIM file
# -------------------------------------------------------
cat("=== Updating BIM file ===\n")

# Replace SNP IDs where we found a match
bim_updated <- bim %>%
  mutate(
    SNP_new = ifelse(!is.na(rsID_matched), rsID_matched, SNP_original)
  ) %>%
  select(CHR, SNP_new, CM, POS, A1, A2)

# Prepare for writing (no column names for BIM)
colnames(bim_updated) <- NULL

# Write updated BIM file
output_bim <- paste0(OUTPUT_PREFIX, ".bim")
fwrite(bim_updated, output_bim, sep = "\t", col.names = FALSE, quote = FALSE)

cat("✓ Wrote updated BIM file:", output_bim, "\n\n")

# -------------------------------------------------------
# 6) Copy BED and FAM files
# -------------------------------------------------------
cat("=== Copying BED and FAM files ===\n")

# Copy BED file (binary, unchanged)
system(paste("cp", paste0(INPUT_PREFIX, ".bed"), paste0(OUTPUT_PREFIX, ".bed")))
cat("✓ Copied BED file\n")

# Copy FAM file (unchanged)
system(paste("cp", paste0(INPUT_PREFIX, ".fam"), paste0(OUTPUT_PREFIX, ".fam")))
cat("✓ Copied FAM file\n\n")

# -------------------------------------------------------
# 7) Summary statistics
# -------------------------------------------------------
cat("=== Summary Statistics ===\n\n")

# Count different types of IDs in updated file
bim_final <- fread(output_bim, header = FALSE)
colnames(bim_final) <- c("CHR", "SNP", "CM", "POS", "A1", "A2")

n_rs <- sum(grepl("^rs", bim_final$SNP))
n_kgp <- sum(grepl("^kgp", bim_final$SNP))
n_snp_known <- sum(grepl("^snp-known", bim_final$SNP))
n_other <- nrow(bim_final) - n_rs - n_kgp - n_snp_known

cat("Updated BIM file summary:\n")
cat("  Total variants:", nrow(bim_final), "\n")
cat("  rsIDs (rs*):", n_rs, "\n")
cat("  KGP IDs (kgp*):", n_kgp, "\n")
cat("  Known SNPs (snp-known*):", n_snp_known, "\n")
cat("  Other IDs:", n_other, "\n\n")

# Show examples of updated variants
cat("Examples of updated variants:\n")
matched_examples <- bim %>%
  filter(!is.na(rsID_matched)) %>%
  select(SNP_original, SNP_new = rsID_matched, CHR, POS) %>%
  head(10)
print(matched_examples)

cat("\n")

# Show examples that couldn't be matched
unmatched_examples <- bim %>%
  filter(is.na(rsID_matched)) %>%
  select(SNP_original, CHR, POS, A1, A2) %>%
  head(5)

if (nrow(unmatched_examples) > 0) {
  cat("Examples of unmatched variants (kept original IDs):\n")
  print(unmatched_examples)
  cat("\n")
}

# -------------------------------------------------------
# 8) Create matching report
# -------------------------------------------------------
report_file <- paste0(OUTPUT_PREFIX, "_rsid_matching_report.txt")

sink(report_file)
cat("rsID Restoration Report\n")
cat("======================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Input files:\n")
cat("  Annotation:", ANNOTATION_FILE, "\n")
cat("  PLINK prefix:", INPUT_PREFIX, "\n\n")
cat("Output files:\n")
cat("  PLINK prefix:", OUTPUT_PREFIX, "\n\n")
cat("Matching statistics:\n")
cat("  Total variants in BIM:", nrow(bim), "\n")
cat("  Variants in annotation:", nrow(annot), "\n")
cat("  Successfully matched:", n_matched, "(", pct_matched, "%)\n")
cat("  Unmatched:", nrow(bim) - n_matched, "\n\n")
cat("Updated BIM composition:\n")
cat("  rsIDs (rs*):", n_rs, "\n")
cat("  KGP IDs (kgp*):", n_kgp, "\n")
cat("  Known SNPs (snp-known*):", n_snp_known, "\n")
cat("  Other IDs:", n_other, "\n\n")
sink()

cat("✓ Saved matching report:", report_file, "\n\n")

# -------------------------------------------------------
# Final message
# -------------------------------------------------------
cat("================================\n")
cat("rsID restoration complete!\n")
cat("================================\n")
cat("Output files:\n")
cat("  ", paste0(OUTPUT_PREFIX, ".bed"), "\n")
cat("  ", paste0(OUTPUT_PREFIX, ".bim"), "\n")
cat("  ", paste0(OUTPUT_PREFIX, ".fam"), "\n")
cat("  ", report_file, "\n\n")
cat("You can now use these files for GWAS analysis.\n")
cat("================================\n")
