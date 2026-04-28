#!/usr/bin/env Rscript

############################################################
# QC PROGRESS TRACKER (with filter metadata)
# ------------------------------------------
# Usage:
#   Rscript qc_tracker.R <plink_prefix> "<step_name>" [qc_csv] [options]
#
# Options (all optional, free text allowed):
#   --snp-miss <value>       # e.g. 0.05 for --geno 0.05
#   --sample-miss <value>    # e.g. 0.05 for --mind 0.05
#   --maf <value>            # e.g. 0.01 for --maf 0.01
#   --hwe <value>            # e.g. 1e-6 for --hwe 1e-6
#   --info <value>           # e.g. 0.8 for post-imputation info filter
#   --other "<free text>"    # any notes, multiple words allowed
#
# Examples:
#   Rscript qc_tracker.R H3_noAmbig "Removed ambiguous SNPs"
#
#   Rscript qc_tracker.R H3_mind0.05 "Sample missingness filter" QC_progress_tracker.csv \
#       --sample-miss 0.05
#
#   Rscript qc_tracker.R H3_maf0.01_hwe1e-6 "MAF + HWE filter" QC_progress_tracker.csv \
#       --maf 0.01 --hwe 1e-6 --snp-miss 0.05 \
#       --other "Autosomes only; removed monomorphic SNPs"
#
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
})

# -------------------------
# 1. Parse command-line args
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage:\n",
    "  Rscript qc_tracker.R <plink_prefix> \"<step_name>\" [qc_csv] [options]\n\n",
    "Options:\n",
    "  --snp-miss <value>    (SNP missingness filter, e.g. 0.05)\n",
    "  --sample-miss <value> (Sample missingness filter, e.g. 0.05)\n",
    "  --maf <value>         (MAF threshold, e.g. 0.01)\n",
    "  --hwe <value>         (HWE p-value threshold, e.g. 1e-6)\n",
    "  --info <value>        (INFO threshold, e.g. 0.8)\n",
    "  --other \"text\"       (Free-text description of other filters)\n"
  )
}

plink_prefix <- args[1]
step_name    <- args[2]

# Default CSV name if not provided
if (length(args) >= 3 && !startsWith(args[3], "--")) {
  qc_csv <- args[3]
  opt_start <- 4
} else {
  qc_csv <- "QC_progress_tracker.csv"
  opt_start <- 3
}

cat(">>> PLINK prefix  :", plink_prefix, "\n")
cat(">>> QC step name  :", step_name, "\n")
cat(">>> Output CSV    :", qc_csv, "\n")

# -------------------------
# 1b. Parse optional filter flags
# -------------------------
snp_missingness_filter    <- NA_character_
sample_missingness_filter <- NA_character_
maf_filter                <- NA_character_
hwe_filter                <- NA_character_
info_filter               <- NA_character_
other_filters             <- NA_character_

i <- opt_start
while (i <= length(args)) {
  arg <- args[i]
  if (!startsWith(arg, "--")) {
    i <- i + 1
    next
  }
  key <- sub("^--", "", arg)

  # Grab next argument as value if present and not another flag
  val <- NA_character_
  if (i + 1 <= length(args) && !startsWith(args[i + 1], "--")) {
    val <- args[i + 1]
    i <- i + 1
  }

  if (key == "snp-miss") {
    snp_missingness_filter <- val
  } else if (key == "sample-miss") {
    sample_missingness_filter <- val
  } else if (key == "maf") {
    maf_filter <- val
  } else if (key == "hwe") {
    hwe_filter <- val
  } else if (key == "info") {
    info_filter <- val
  } else if (key == "other") {
    # if multiple --other are used, concatenate
    if (is.na(other_filters) || other_filters == "") {
      other_filters <- val
    } else {
      other_filters <- paste(other_filters, val, sep = " | ")
    }
  } else {
    warning("Unrecognized option: --", key, " (ignored)")
  }

  i <- i + 1
}

cat("\n>>> Parsed filter metadata:\n")
cat("   SNP missingness  (geno): ", snp_missingness_filter, "\n")
cat("   Sample missingness(mind): ", sample_missingness_filter, "\n")
cat("   MAF threshold          : ", maf_filter, "\n")
cat("   HWE threshold          : ", hwe_filter, "\n")
cat("   INFO threshold         : ", info_filter, "\n")
cat("   Other filters          : ", other_filters, "\n\n")

# -------------------------
# 2. Read BIM and FAM
# -------------------------
bim_file   <- paste0(plink_prefix, ".bim")
fam_file   <- paste0(plink_prefix, ".fam")
imiss_file <- paste0(plink_prefix, ".imiss")

if (!file.exists(bim_file)) stop("Cannot find BIM file: ", bim_file)
if (!file.exists(fam_file)) stop("Cannot find FAM file: ", fam_file)

bim <- data.table::fread(
  bim_file,
  header   = FALSE,
  col.names = c("CHR", "SNP", "GD", "BP", "A1", "A2")
)

fam <- data.table::fread(
  fam_file,
  header   = FALSE,
  col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")
)

n_variants <- nrow(bim)
n_samples  <- nrow(fam)

# -------------------------
# 3. Cases / controls & sex
# -------------------------
fam <- fam %>%
  mutate(
    SEX        = as.numeric(SEX),
    PHENO      = as.numeric(PHENO),
    is_male    = ifelse(SEX == 1, 1L, 0L),
    is_female  = ifelse(SEX == 2, 1L, 0L),
    is_case    = ifelse(PHENO == 2, 1L, 0L),
    is_control = ifelse(PHENO == 1, 1L, 0L),
    sex_missing   = ifelse(is.na(SEX) | !(SEX %in% c(1, 2)), 1L, 0L),
    pheno_missing = ifelse(is.na(PHENO) | !(PHENO %in% c(1, 2)), 1L, 0L)
  )

n_male          <- sum(fam$is_male, na.rm = TRUE)
n_female        <- sum(fam$is_female, na.rm = TRUE)
n_sex_missing   <- sum(fam$sex_missing, na.rm = TRUE)
n_cases         <- sum(fam$is_case, na.rm = TRUE)
n_controls      <- sum(fam$is_control, na.rm = TRUE)
n_pheno_missing <- sum(fam$pheno_missing, na.rm = TRUE)

cat("Samples:", n_samples,
    " | Cases:", n_cases,
    " | Controls:", n_controls, "\n")
cat("Males:", n_male,
    " | Females:", n_female,
    " | Sex missing:", n_sex_missing, "\n")
cat("Variants:", n_variants, "\n\n")

# -------------------------
# 4. Genotyping rate (from .imiss)
# -------------------------
genotyping_rate <- NA_real_

if (file.exists(imiss_file)) {
  imiss <- data.table::fread(imiss_file, header = TRUE)
  if (all(c("N_MISS", "N_GENO") %in% names(imiss))) {
    total_miss <- sum(imiss$N_MISS, na.rm = TRUE)
    total_geno <- sum(imiss$N_GENO, na.rm = TRUE)
    if (total_geno > 0) {
      genotyping_rate <- 1 - (total_miss / total_geno)
    }
  } else {
    warning("File found but missing N_MISS/N_GENO columns: ", imiss_file)
  }
} else {
  warning("Missingness file not found (no genotyping rate): ", imiss_file)
}

cat("Genotyping rate (overall call rate): ",
    ifelse(is.na(genotyping_rate),
           "NA (no .imiss)",
           sprintf("%.6f", genotyping_rate)),
    "\n\n")

# -------------------------
# 5. Read previous QC CSV (if exists) & compute exclusions
# -------------------------
if (file.exists(qc_csv)) {
  prev <- readr::read_csv(qc_csv, show_col_types = FALSE)

  last_step_index  <- tail(prev$step_index, 1)
  last_n_variants  <- tail(prev$n_variants, 1)
  last_total_excl  <- tail(prev$total_variants_excluded, 1)

  step_index <- last_step_index + 1

  variants_excluded_since_last <- last_n_variants - n_variants
  if (is.na(variants_excluded_since_last)) {
    variants_excluded_since_last <- NA_real_
  }

  if (is.na(last_total_excl)) last_total_excl <- 0
  total_variants_excluded <- last_total_excl +
    max(variants_excluded_since_last, 0, na.rm = TRUE)

} else {
  step_index <- 1
  variants_excluded_since_last <- NA_real_
  total_variants_excluded      <- 0
}

cat("Step index:", step_index, "\n")
cat("Variants excluded since last step:",
    ifelse(is.na(variants_excluded_since_last),
           "NA (first step)",
           variants_excluded_since_last),
    "\n")
cat("Total variants excluded (cumulative):",
    total_variants_excluded, "\n\n")

# -------------------------
# 6. Build summary row
# -------------------------
qc_row <- tibble::tibble(
  date_time                   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  step_index                  = step_index,
  step_name                   = step_name,
  plink_prefix                = plink_prefix,
  n_variants                  = n_variants,
  n_samples                   = n_samples,
  n_cases                     = n_cases,
  n_controls                  = n_controls,
  n_pheno_missing             = n_pheno_missing,
  n_male                      = n_male,
  n_female                    = n_female,
  n_sex_missing               = n_sex_missing,
  genotyping_rate             = genotyping_rate,
  variants_excluded_since_last = variants_excluded_since_last,
  total_variants_excluded     = total_variants_excluded,
  snp_missingness_filter      = snp_missingness_filter,
  sample_missingness_filter   = sample_missingness_filter,
  maf_filter                  = maf_filter,
  hwe_filter                  = hwe_filter,
  info_filter                 = info_filter,
  other_filters               = other_filters
)

qc_row_df <- as.data.frame(qc_row)

# -------------------------
# 7. Append / create CSV
# -------------------------
if (file.exists(qc_csv)) {
  write.table(
    qc_row_df,
    file      = qc_csv,
    sep       = ",",
    row.names = FALSE,
    col.names = FALSE,
    append    = TRUE,
    quote     = TRUE
  )
} else {
  write.table(
    qc_row_df,
    file      = qc_csv,
    sep       = ",",
    row.names = FALSE,
    col.names = TRUE,
    append    = FALSE,
    quote     = TRUE
  )
}

cat("✅ QC summary row written to:", qc_csv, "\n")
