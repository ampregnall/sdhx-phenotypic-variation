library(data.table)

dt <- fread("data/pmbb-variants/pmbb-carrier-status-annotated.txt")
dt[, dob := sample_date - round(age * 365.25)]

# load date of death date and merge into existing data.table
dt_death <- fread(
  file = "/static/PMBB/PMBB-Release-2026-4.0/Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_death.txt"
)

dt <- merge(dt, dt_death, by = "person_id", all.x = TRUE)
dt <- dt[, age_death := as.numeric(death_date - dob) / 365.25]

# load diagnosis data and find earliest date of cancer diagnosis
dt_diagnoses <- fread(
  file = "/static/PMBB/PMBB-Release-2026-4.0/Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_conditions_icd10.txt"
)

dt_diagnoses <- dt_diagnoses[
  grepl(pattern = "^C", condition_source_value)
][,
  .(
    diagnosis_date = min(condition_start_date)
  ),
  by = person_id
]

dt <- merge(dt, dt_diagnoses, by = "person_id", all.x = TRUE)
dt <- dt[, diagnosis_age := as.numeric(diagnosis_date - dob) / 365.25]

dt_visits <- fread(
  file = "/static/PMBB/PMBB-Release-2026-4.0/Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_visit_occurrence.txt"
)

dt_visits <- dt_visits[,
  .(
    first_visit_date = min(visit_start_date),
    last_visit_date = max(visit_start_date)
  ),
  by = person_id
]

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
dt <- merge(dt, dt_visits, by = "person_id", all.x = TRUE)

fwrite(
  x = dt,
  file = "data/pmbb/pmbb-time-to-event.txt",
  sep = "\t"
)
