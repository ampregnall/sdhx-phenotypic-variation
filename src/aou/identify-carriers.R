library(data.table)

# Load demographic data
dt_demographics <- fread(
  "data/all-of-us/phenotypes/aou-demographics.txt"
)

N_GENOMICS <- nrow(dt_demographics[has_whole_genome_variant == TRUE])

dt <- rbindlist(lapply(list.files(path = "data/all-of-us/variants/", full.names = TRUE), fread))
dt[, ID := sprintf("chr%s_%s_%s_%s", chrom, pos, ref, alt)]
dt <- unique(dt, by = c("sample", "chrom", "pos", "ref", "alt"))

# load classified variants
dt_variants <- fread(
  "results/variant-classification/aou-variants-classified.txt"
)

dt <- dt[ID %in% dt_variants[Consequence == 'pLoF']$ID]
dt[dt_variants, gene := i.Gene, on = "ID"]


dt_prevalence <- dt[,
                    .(
                      n_carriers = length(unique(sample)),
                      prevalence = (length(unique(sample)) / N_GENOMICS) * 100
                      ),
                    by = gene
]


dt_genotypes <- dcast(
  data = dt,
  formula = sample ~ gene,
  fun.aggregate = function(x) return(1),
  fill = 0
)

dt_demographics <- dt_demographics[has_whole_genome_variant & has_ehr_data]
dt_demographics <- dt_demographics[,
                                   .(
                                     person_id,
                                     age,
                                     sex = T_DISP_sex_at_birth,
                                     year_of_birth
                                   )
]

dt_demographics <- merge(
  dt_demographics,
  dt_genotypes,
  by.x = "person_id",
  by.y = "sample",
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
  "data/all-of-us/aou-carrier-status-annotated.txt",
  sep = "\t"
)

fwrite(
  dt_prevalence,
  "results/prevalence/aou-prevalence.txt",
  sep = "\t"
)
