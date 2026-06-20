# Host Genetic Determinants of Lineage-Specific TB Susceptibility

## Overview
This repository contains analysis scripts and code for the case-case genome-wide association study (GWAS) comparing *M. africanum* and *M. tuberculosis sensu stricto* susceptibility in a Ghanaian cohort (N=232 TB cases).

## Citation
If you use this code, please cite:


GitHub: https://github.com/TheoAfum/TB-Maf-Mtbss-GWAS

## Contents

### Scripts (`scripts/`)
- **01_qc_pipeline.qmd** - Genotype QC, liftover, imputation prep and post-imputation steps (PLINK/bcftools)
- **02_pca_BigSNP.R** - Population structure / PCA (bigsnpr projection onto 1000 Genomes)
- **02.5_admixture.sh** - ADMIXTURE runs for K=2–10 with cross-validation
- **03_case_case_gwas.sh** - Case-case GWAS (logistic regression, Firth-fallback)
- **04_magma_analysis.sh** - Genome-wide gene-based and gene-set association (MAGMA)
- **05_pathway_with_geneset_plot.R** - Pathway enrichment (EnrichR) with plots
- **geneset_plots_BH.R** - GSEA forest plots from BH-corrected MAGMA gene-set results
- **qc_tracker.R** - Per-step QC tracking table (counts, filters, exclusions)

### Figures (`figures/`)
- **GWAS.R** - Manhattan and QQ plots
- **locuszoom.R** - Regional (LocusZoom-style) plots
- **Admixtfure.R** - ADMIXTURE bar plots and CV-error plot
- **genomewide_magma.R** - Integrated multi-panel gene/pathway figure

### Configs (`configs/`)
- **install_dependencies.R** - Install required R packages
- **post_impute_QC_and_merge.sh** - Concatenate, filter and merge imputed VCFs to a PLINK dataset
- **restore_rsid.R** - Restore rsIDs to PLINK files from an annotation file
- **covert_genesets_to_entrez.R** - Convert gene-set symbols to Entrez IDs
- **tb_immune_genesets_symbols.txt** - Curated TB immune gene sets

### Configuration / paths
Scripts use **relative paths by default** (rooted at the repository, e.g. `data/`, `results/`, `configs/`) and create output folders as needed. No paths are hardcoded to a specific machine. Any input or output location can be overridden with an environment variable without editing the code, for example:

```bash
# R scripts read overrides via Sys.getenv(); shell scripts via ${VAR:-default}
GWAS_RESULTS=/my/path/gwas.txt Rscript figures/GWAS.R
GENO_PREFIX=/my/path/genotypes OUTPUT_DIR=/my/results bash scripts/04_magma_analysis.sh
```

The variable name is shown next to each path in the script's configuration block. Place input data under `data/` (or point the relevant variable at your own location) to run end-to-end.

### Data
Individual-level genotype data are available through the European Genome-phenotype Archive (EGA) under managed access. 
- **EGA Accession:** [INSERT WHEN ASSIGNED]
- Data Access Request: https://ega-archive.org/request-access/

Summary statistics and metadata are provided in the `results/` folder.

## Requirements
- R 4.5.1 or later
- PLINK 1.9 & 2.0
- MAGMA v1.10
- Python 3.x
- See individual scripts for package dependencies

## Installation & Usage

### Quick Start
```bash
# Clone this repository
git clone https://github.com/TheoAfum/TB-Maf-Mtbss-GWAS.git
cd TB-Maf-Mtbss-GWAS

# Run analysis pipeline (requires input data under data/)
# See "Configuration / paths" below to override locations.
Rscript configs/install_dependencies.R
# QC / imputation steps are documented in scripts/01_qc_pipeline.qmd
Rscript scripts/02_pca_BigSNP.R
bash    scripts/03_case_case_gwas.sh
bash    scripts/04_magma_analysis.sh
Rscript scripts/05_pathway_with_geneset_plot.R
# ... then figures in figures/
```

### Detailed Instructions
See each script's header comments for:
- Input file formats
- Required software versions
- Expected runtime
- Output descriptions

## Key Findings
- No genome-wide significant SNPs (GWAS)
- B9D1 significant by gene-based analysis (MAGMA, P=2.01×10⁻⁶)
- Three FDR-significant pathways (ORA):
  - Purinergic nucleotide receptor activity
  - Deoxycytidine deaminase activity
  - JAK-STAT signaling (MAGMA GSA)

## Ethics & Data Sharing
This study was approved by:
- Ghana Health Service Ethical Review Committee (GHS-ERC: 010/08/23)
- Noguchi Memorial Institute Institutional Review Board (NMIMR-IRB: 072/19-20)

Individual-level data are managed-access through EGA to protect participant confidentiality.

## Contact
For questions or collaboration:
- **Dorothy Yeboah-Manu** - dyeboah-manu@noguchi.ug.edu.gh
- **Theophilus Afum** - tafum@noguchi.ug.edu.gh

## License
MIT License - See LICENSE file for details
