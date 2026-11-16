library(data.table)
library(logistf)
source("src/utils/firth-regression.R")

dt_model <- fread("data/all-of-us/aou-carrier-status-annotated.txt")
dt_diagnoses <- fread(
  "data/all-of-us/phenotypes/aou-icd10-diagnoses.txt"
)

dt_diagnoses <- dt_diagnoses[,
                             .(
                               breast = any(grepl(pattern = "^C50", x = condition_source_value)),
                               renal_cell = any(
                                 condition_source_value %in% c("C64.1", "C64.2", "C64.9", "C64")
                               ), # AOU contains C64, C64.9, C64.2, & C64.1
                               gist = any(grepl(pattern = "^C49", x = condition_source_value)),
                               pcc = any(condition_source_value %in% c("C74.10", "C74.11", "C74.12")),
                               pgl = any(condition_source_value %in% c("C75.4", "C75.5")),
                               pituitary = any(condition_source_value %in% c("D35.2")),
                               adrenal = any(
                                 condition_source_value %in% c("D35.00", "D35.01", "D35.02")
                               ),
                               benign_pgl = any(
                                 condition_source_value %in% c("D35.6", "D44.6", "D44.7")
                               )
                             ),
                             by = person_id
]

dt <- merge(dt_model, dt_diagnoses, by = "person_id")

# Load smoking information
dt_smoking <- fread("data/all-of-us/phenotypes/aou-smoking.txt")

dt <- merge(dt, dt_smoking, by = "person_id", all.x = TRUE)
dt <- dt[!is.na(smoking) & sex %in% c("Male", "Female")]
dt <- dt[, age := 2025 - year_of_birth]

GENES <- c("SDHA", "SDHB", "SDHC", "SDHD", "SDHAF2", "MAX")
COVARS <- c("age", "sex", "smoking")

PHENOS <- c(
  breast = "Breast",
  renal_cell = "Renal cell carcinoma",
  gist = "GIST",
  pcc = "Pheochromocytoma",
  pgl = "Paraganglioma",
  pituitary = "Pituitary adenoma",
  adrenal = "Adrenal adenoma",
  benign_pgl = "Benign paraganglioma"
)

BIOBANK <- "aou"
out_dir <- "results/known-associations"
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
  file.path(out_dir, sprintf("%s-known-associations.txt", BIOBANK)),
  sep = "\t",
  na = "NA"
)
