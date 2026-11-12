library(data.table)
source("src/utils/classify-variants.R")

GENES <- c(
    "SDHA",
    "SDHB",
    "SDHC",
    "SDHD",
    "SDHAF2",
    "TMEM127",
    "MAX"
)

# Load variants annotated with VEP/ClinVar/AutoGVP
dt_variants <- fread(
    file = "data/pmbb-variants/PMBB_SDHx_annotated_11-03-2026.txt_autoGVP_all.txt",
    colClasses = c("exon_vep" = "character")
)

dt_variants <- classify_variants(
    dt = dt_variants,
    allele_freq = 0.01,
    clinvar_annotation = "clinvar_clinsig_autogvp",
    gnomad_frequency = "gnomad41_genome_fafmax_faf99_max",
    autogvp_annotation = "autogvp_call",
    vep_annotation = "variant_classification_vep",
    exon_col = "exon_vep",
    am_annotation = "AlphaMissense_score"
)

dt_variants <- dt_variants[,
    .(
        CHROM,
        POS,
        REF,
        ALT,
        ID,
        Gene = gene_symbol_vep,
        VEP = variant_classification_vep,
        ClinVar = clinvar_clinsig_autogvp,
        AutoGVP = autogvp_call,
        AlphaMissense = AlphaMissense_score,
        Consequence
    )
]

dt_variants <- dt_variants[Gene %in% GENES]
dt_variants <- unique(dt_variants, by = "ID")

dt_summary <- dt_variants[, .N, by = .(Gene, Consequence)]

fwrite(
    x = dt_variants,
    file = "results/variant-classification/pmbb-variants-classified.txt",
    sep = "\t"
)
fwrite(
    x = dt_summary,
    file = "results/variant-classification/pmbb-variant-classification-summary.txt",
    sep = "\t"
)
