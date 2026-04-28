# ADMIXTURE Analysis Pipeline - Complete Guide

## 📋 Overview
This pipeline runs ADMIXTURE analysis for K=2 to K=10 and creates publication-quality Nature-style plots.

## 📁 Required Files
```
H3_QC8_ibdClean.bed
H3_QC8_ibdClean.bim
H3_QC8_ibdClean.fam
H3_pheno.phe.txt (optional - for sorted plots)
```

## 🚀 Quick Start

### Step 1: Run ADMIXTURE Analysis
```bash
bash run_admixture_analysis.sh
```

This will:
- Run ADMIXTURE for K=2 to K=10
- Use cross-validation to find optimal K
- Create `ADMIXTURE_results/` directory with all output files

**Time estimate**: 10-60 minutes depending on data size and number of SNPs

### Step 2: Create Visualizations
```bash
Rscript plot_admixture.R
```

This will create:
- Individual plots for each K (K=2 to K=10)
- Combined multi-K panel plot
- Cross-validation error plot
- Sorted plots by ethnicity/strain (if phenotype file exists)

## 📊 Output Files

### ADMIXTURE Results (`ADMIXTURE_results/`)
```
H3_QC8_ibdClean.K2.Q     - Ancestry proportions for K=2
H3_QC8_ibdClean.K2.P     - Allele frequencies for K=2
log2.out                 - ADMIXTURE log for K=2
cv_errors.txt            - Cross-validation errors for all K
```

### Plots (`ADMIXTURE_plots/`)
```
Admixture_K2.png/pdf                           - Individual K=2 plot
Admixture_K3.png/pdf                           - Individual K=3 plot
...
Admixture_K10.png/pdf                          - Individual K=10 plot
Admixture_Combined_K2_K10.png/pdf              - All K values stacked
CV_error_plot.png/pdf                          - Cross-validation plot
Ancestry_proportions_K[best].csv               - Summary statistics

# If phenotype file exists:
Admixture_K2_by_STRAIN.png                     - K=2 sorted by strain
Admixture_K2_by_ETHNICITY.png                  - K=2 sorted by ethnicity
```

## ⚙️ Customization

### Modify K range
Edit `run_admixture_analysis.sh`:
```bash
MIN_K=2
MAX_K=10
```

### Change input files
Edit both scripts:
```bash
# In run_admixture_analysis.sh
INPUT_PREFIX="H3_QC8_ibdClean"

# In plot_admixture.R
INPUT_PREFIX <- "H3_QC8_ibdClean"
PHENO_FILE <- "H3_pheno.phe.txt"
```

### Adjust plot appearance
Edit `plot_admixture.R`:
```r
# Plot dimensions
width = 10, height = 3

# DPI (resolution)
dpi = 600

# Colors
get_colors <- function(k) {
  # Modify color palettes here
}
```

## 🎨 Plot Features

### Nature-Style Design
- Clean, professional appearance
- Colorblind-friendly palettes
- High resolution (600 DPI)
- Both PNG and PDF formats
- No unnecessary elements

### Sorting Options
If phenotype file is provided, samples are automatically sorted by:
- **Strain** (MAF, Mtbss, Control)
- **Ethnicity** (Akan, Ewe, Ga, Other)
- Vertical dividers separate groups

### Cross-Validation Plot
- Shows CV error for each K
- Best K highlighted in red
- Helps determine optimal number of populations

## 📝 Requirements

### Software
```bash
# ADMIXTURE (must be in PATH)
admixture --version

# R packages
install.packages(c("data.table", "dplyr", "tidyr", "ggplot2", 
                   "cowplot", "RColorBrewer", "viridis"))
```

### System Requirements
- **Memory**: ~4-8 GB (depends on data size)
- **CPU**: Multi-core recommended (script uses 4 threads)
- **Disk**: ~1 GB for results

## 🔍 Interpreting Results

### Cross-Validation Error
- **Lower is better**
- Look for "elbow" in CV error plot
- Best K is where CV error is minimized

### Ancestry Plots
- Each color represents one ancestral population
- Bar height shows proportion of ancestry
- Individuals sorted for easier interpretation

### Common K Values
- **K=2**: Major population splits
- **K=3-5**: Fine-scale population structure
- **K>5**: May show overfitting or very fine structure

## 🐛 Troubleshooting

### "admixture: command not found"
```bash
# Install ADMIXTURE
wget http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
tar -xvf admixture_linux-1.3.0.tar.gz
sudo mv admixture /usr/local/bin/
```

### "Error: Cannot find .Q files"
- Check that ADMIXTURE ran successfully
- Look in `ADMIXTURE_results/` for .Q files
- Check log files for errors

### Plots look wrong
- Ensure phenotype file has FID and IID columns
- Check that sample order matches FAM file
- Verify phenotype column names (ETHNICITY, STRAIN)

## 📚 References

- Alexander, D.H., Novembre, J. & Lange, K. (2009) Fast model-based estimation of ancestry in unrelated individuals. *Genome Research* 19:1655-1664.

## 💡 Tips

1. **Pruning**: Use LD-pruned data for ADMIXTURE (you already have this!)
2. **CV error**: Run multiple times if unstable
3. **Best K**: Don't just trust CV - look at biological meaning
4. **Large K**: K>10 may be computationally expensive
5. **Publication**: Use PDF outputs for manuscripts

## 🎯 Next Steps

After running ADMIXTURE:
1. Examine CV error plot to determine best K
2. Look at ancestry patterns for biological interpretation
3. Compare with PCA results
4. Include best K plot in publication
5. Report CV errors in methods section

## 📧 Need Help?

Check the log files in `ADMIXTURE_results/log*.out` for detailed error messages.
