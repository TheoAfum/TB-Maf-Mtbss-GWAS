#!/bin/bash
set -euo pipefail

cd ~/H3_imputation/final/imputated_data

# ---------------- CONFIG ----------------
INFO_FIELD="R2"      # from VCF header: Estimated Imputation Accuracy (R-square)
INFO_THRESH=0.3      # you can increase to 0.5 or 0.8 if you want stricter QC
MAF_PLINK=0.01       # MAF threshold for PLINK QC
GENO_PLINK=0.02      # SNP missingness
MIND_PLINK=0.02      # sample missingness
# ----------------------------------------

echo "=== 1. Concatenate chunk VCFs per chromosome ==="

for chr in {1..22}; do
  echo ""
  echo ">>> Chromosome ${chr}: concatenating chunk VCFs"

  ls chr_${chr}/*.vcf.gz > chr_${chr}_vcf_list.txt

  bcftools concat -f chr_${chr}_vcf_list.txt -Oz -o imputed_chr${chr}.vcf.gz
  tabix -p vcf imputed_chr${chr}.vcf.gz
done

echo ""
echo "=== 2. VCF-level post-imputation QC (INFO/R2 + SNPs only) ==="

for chr in {1..22}; do
  IN_VCF="imputed_chr${chr}.vcf.gz"
  OUT_VCF="imputed_chr${chr}_filt.vcf.gz"

  echo ""
  echo ">>> Chromosome ${chr}: filtering ${IN_VCF} -> ${OUT_VCF}"

  bcftools view \
    -i "INFO/${INFO_FIELD}>=${INFO_THRESH}" \
    -m2 -M2 -v snps \
    -Oz -o "${OUT_VCF}" \
    "${IN_VCF}"

  tabix -p vcf "${OUT_VCF}"
done

echo ""
echo "=== 3. Convert filtered VCFs to PLINK per chromosome ==="

for chr in {1..22}; do
  VCF="imputed_chr${chr}_filt.vcf.gz"
  OUT="imputed_chr${chr}_plink"

  echo ""
  echo ">>> Chromosome ${chr}: converting ${VCF} -> ${OUT}.bed/bim/fam"

  plink \
    --vcf "${VCF}" \
    --double-id \
    --keep-allele-order \
    --make-bed \
    --out "${OUT}"
done

echo ""
echo "=== 4. Merge all chromosomes into one PLINK dataset ==="

BASE_PREFIX="imputed_chr1_plink"
MERGE_LIST="merge_list.txt"

rm -f "${MERGE_LIST}"

for chr in {2..22}; do
  echo "imputed_chr${chr}_plink.bed imputed_chr${chr}_plink.bim imputed_chr${chr}_plink.fam" >> "${MERGE_LIST}"
done

plink \
  --bfile "${BASE_PREFIX}" \
  --merge-list "${MERGE_LIST}" \
  --make-bed \
  --out imputed_all_chr_merged

echo ""
echo "=== 5. Final PLINK QC on merged imputed dataset ==="

plink \
  --bfile imputed_all_chr_merged \
  --geno ${GENO_PLINK} \
  --mind ${MIND_PLINK} \
  --maf ${MAF_PLINK} \
  --hwe 1e-6 midp \
  --make-bed \
  --out imputed_all_chr_QC

echo ""
echo "All done ✅"
echo "Final GWAS-ready dataset: imputed_all_chr_QC.{bed,bim,fam}"
