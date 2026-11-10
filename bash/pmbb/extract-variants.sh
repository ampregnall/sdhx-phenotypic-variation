#!/bin/bash
#BSUB -J extract-sdhx-variants
#BSUB -o logs/variant_extraction.%J.out
#BSUB -e logs/variant_extraction.%J.err
#BSUB -q voltron_normal
#BSUB -M 10000
#BSUB -n 1

cd /project/damrauer_shared/Users/ampregnall/SDHx-Phenotypic-Variation
OUTDIR=data/pmbb-variants/work
mkdir -p $OUTDIR

module load plink/2.0

while IFS=$' \t\r' read -r CHR START END GENE
do
    echo "Running extraction for: $GENE"
    plink2 --bfile "/static/PMBB/PMBB-Release-2026-4.0/Exome/by_chromosome_GL/PMBB-Release-2026-4.0_genetic_exome.gl.chr${CHR}" \
      --chr $CHR \
      --from-bp $START \
      --to-bp $END \
      --export vcf bgz \
      --out ${OUTDIR}/PMBB_${GENE}
done < assets/gene-locations.txt
