# Host Genetic Determinants of Lineage-Specific TB Susceptibility

## Overview
This repository contains analysis scripts and code for the case-case genome-wide association study (GWAS) comparing *M. africanum* and *M. tuberculosis sensu stricto* susceptibility in a Ghanaian cohort (N=232 TB cases).

## Citation
If you use this code, please cite:


GitHub: https://github.com/TheoAfum/TB-Maf-Mtbss-GWAS

## Contents

### Scripts
- **01_qc_pipeline.R** - Genotype quality control (PLINK)
- **02_pca_admixture.R** - Population structure analysis
- **03_gwas_analysis.sh** - Case-case GWAS (logistic regression, Firth-fallback)
- **04_magma_analysis.R** - Gene-based association testing
- **05_pathway_enrichment.R** - Pathway enrichment using EnrichR & MAGMA GSA
- **06_visualization.R** - Generate publication figures

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
git clone https://github.com/YOUR-USERNAME/TB-Lineage-GWAS.git
cd TB-Lineage-GWAS

# Run analysis pipeline (requires input data)
Rscript scripts/01_qc_pipeline.R
Rscript scripts/02_pca_admixture.R
# ... etc
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
- **Theophilus Afum** - tafum@noguchi.ug.edugh

## License
MIT License - See LICENSE file for details
