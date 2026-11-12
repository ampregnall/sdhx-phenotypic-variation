library(data.table)

GENES <- c(
    "SDHA",
    "SDHB",
    "SDHC",
    "SDHD",
    "SDHAF2",
    "TMEM127",
    "MAX"
)

# Load demographic data
dt_demographics <- fread(
    "/static/PMBB/PMBB-Release-2026-4.0/Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_covariates.txt"
)
PMBB_SAMPLE_SIZE <- nrow(dt_demographics)

# Load genotypes
dt_genotypes <- rbindlist(lapply(GENES, function(x) {
    fread(
        file = sprintf("data/pmbb-variants/work/PMBB_%s_variants.txt", x),
        header = FALSE,
        col.names = c("person_id", "chr", "pos", "ref", "alt", "genotype")
    )[genotype %in% c("0/1", "1/1")][, gene := x]
}))

# Create variant identifier in format of classification table
dt_genotypes[, ID := sprintf(fmt = "chr%s_%s_%s_%s", chr, pos, ref, alt)]

# load classified variants
dt_variants <- fread(
    "results/variant-classification/pmbb-variants-classified.txt"
)

dt_genotypes <- dt_genotypes[ID %in% dt_variants[Consequence == 'pLoF']$ID]

dt_prevalence <- dt_genotypes[,
    .(
        n_carriers = length(unique(person_id)),
        prevalence = (length(unique(person_id)) / PMBB_SAMPLE_SIZE) * 100
    ),
    by = gene
]

dt_genotypes <- dcast(
    data = dt_genotypes,
    formula = person_id ~ gene,
    fun.aggregate = function(x) return(1),
    fill = 0
)

dt_genotypes[, person_id := tstrsplit(person_id, "_", fixed = TRUE, keep = 1L)]

dt_demographics <- dt_demographics[,
    .(
        person_id,
        age = sample_age,
        sex = sequenced_gender,
        sample_date,
        population = Class,
        crep_highrisk_flag
    )
]

dt_demographics <- merge(
    dt_demographics,
    dt_genotypes,
    by = "person_id",
    all.x = TRUE
)
setnafill(
    dt_demographics,
    type = "const",
    fill = 0,
    cols = c("SDHA", "SDHB", "SDHC", "SDHD", "SDHAF2", "MAX")
)

fwrite(
    dt_demographics,
    "data/pmbb-variants/pmbb-carrier-status-annotated.txt",
    sep = "\t"
)

fwrite(
    dt_prevalence,
    "results/prevalence/pmbb-prevalence.txt",
    sep = "\t"
)
