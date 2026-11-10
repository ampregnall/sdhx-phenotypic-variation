#!/bin/bash
#BSUB -J query-sdhx-variants
#BSUB -o logs/query_variants.%J.out
#BSUB -e logs/query_variants.%J.err
#BSUB -q voltron_normal
#BSUB -M 10000
#BSUB -n 1

cd /project/damrauer_shared/Users/ampregnall/SDHx-Phenotypic-Variation
OUTDIR=data/pmbb-variants/work
mkdir -p $OUTDIR

module load bcftools

for GENE in SDHA SDHB SDHC SDHD SDHAF2 TMEM127 MAX
do
    bcftools query \
        -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\n]' \
        "${OUTDIR}/PMBB_${GENE}.vcf.gz" \
        -o "${OUTDIR}/PMBB_${GENE}_variants.txt"
done