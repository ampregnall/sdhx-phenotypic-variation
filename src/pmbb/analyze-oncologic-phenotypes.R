library(data.table)
library(glue)
source("src/utils/firth-regression.R")

dt_model <- fread("data/pmbb-variants/pmbb-carrier-status-annotated.txt")

dt_diagnoses <- fread(
    "/static/PMBB/PMBB-Release-2026-4.0/Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_conditions_icd10.txt"
)

dt_diagnoses <- dt_diagnoses[,
    .(
        any_cancer = as.numeric(any(grepl(
            pattern = "^C",
            condition_source_value
        ))),
        head_neck = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 0:14)}"), collapse = "|"),
            condition_source_value
        ))),
        digestive = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 15:26)}"), collapse = "|"),
            condition_source_value
        ))),
        respiratory = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 30:29)}"), collapse = "|"),
            condition_source_value
        ))),
        bone = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 40:41)}"), collapse = "|"),
            condition_source_value
        ))),
        melanoma = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 43)}"), collapse = "|"),
            condition_source_value
        ))),
        skin = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 44)}"), collapse = "|"),
            condition_source_value
        ))),
        soft_tissue = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 45:49)}"), collapse = "|"),
            condition_source_value
        ))),
        female_genitals = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 51:58)}"), collapse = "|"),
            condition_source_value
        ))),
        male_genitals = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 60:63)}"), collapse = "|"),
            condition_source_value
        ))),
        urinary = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 64:68)}"), collapse = "|"),
            condition_source_value
        ))),
        cns = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 69:72)}"), collapse = "|"),
            condition_source_value
        ))),
        endocrine = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 73:75)}"), collapse = "|"),
            condition_source_value
        ))),
        illdefined = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 76:80)}"), collapse = "|"),
            condition_source_value
        ))),
        lymphoid = as.numeric(any(grepl(
            pattern = paste(glue("C{sprintf('%02d', 81:96)}"), collapse = "|"),
            condition_source_value
        )))
    ),
    by = person_id
]

dt <- merge(dt_model, dt_diagnoses, by = "person_id")
dt <- dt[crep_highrisk_flag == 0]

# Load smoking information
dt_smoking <- fread(
    file = "/static/PMBB/PMBB-Release-2026-4.0/Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_substance_use.txt"
)

dt_smoking <- dt_smoking[,
    .(
        smoking = as.numeric(any((grepl("Yes|Quit", TOBACCO_USER))))
    ),
    by = person_id
]

dt <- merge(dt, dt_smoking, by = "person_id", all.x = TRUE)
dt <- dt[!is.na(smoking)]

GENES <- c("SDHA", "SDHB", "SDHC", "SDHD", "SDHAF2", "MAX")
COVARS <- c("age", "sex", "smoking")

PHENOS <- c(
    any_cancer = "Any cancer",
    head_neck = "Head and neck",
    digestive = "Digestive",
    respiratory = "Respiratory",
    bone = "Bone",
    skin = "Basal and squamous cell carcinoma",
    melanoma = "Melanoma",
    soft_tissue = "Soft tissue",
    female_genitals = "Female genital",
    male_genitals = "Male genital",
    urinary = "Urinary",
    cns = "Central nervous system",
    endocrine = "Endocrine",
    illdefined = "Illdefined",
    lymphoid = "Lymphoid and myeloid"
)

BIOBANK <- "pmbb"
out_dir <- "results/oncologic-associations"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

grid <- CJ(phenotype = names(PHENOS), gene = GENES, sorted = FALSE)

res <- rbindlist(Map(
    function(p, g) fit_one(dt, p, g),
    grid$phenotype,
    grid$gene
))

res[, biobank := BIOBANK]
setcolorder(res, "biobank")

fwrite(
    res,
    file.path(out_dir, sprintf("%s-oncologic-associations.txt", BIOBANK)),
    sep = "\t",
    na = "NA"
)
