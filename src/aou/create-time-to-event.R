library(data.table)
library(lubridate)

dt_model <- fread("data/all-of-us/aou-carrier-status-annotated.txt")
dt_enrollment <- fread("data/all-of-us/phenotypes/aou-enrollment-date.txt")
dt_death <- fread("data/all-of-us/phenotypes/aou-death.txt")

# Add enrollment date
dt_model <- merge(dt_model, dt_enrollment, by = "person_id", all.x = TRUE)
dt_model <- merge(dt_model, dt_death, by = "person_id", all.x = TRUE)

# Get diagnosis dates
dt_visits <- fread("data/all-of-us/phenotypes/aou-icd10-diagnoses.txt")

cap <- as.IDate("2025-01-01")

dt_visits <- dt_visits[, {
  dx_dates <- condition_start_date[grepl("^C", condition_source_value)]
  .(
    first_visit_date = min(condition_start_date, na.rm = TRUE),
    last_visit_date  = min(max(condition_start_date, na.rm = TRUE), cap),
    diagnosis_date   = if (length(dx_dates)) min(dx_dates, na.rm = TRUE) else dx_dates[NA_integer_]
  )
}, by = person_id]


dt_model <- merge(dt_model, dt_visits, by = "person_id", all.x = TRUE)

dt_model <- dt_model[, dob := as.IDate(age)]
dt_model <- dt_model[, age_death := year(death_date) - year_of_birth]
dt_model <- dt_model[, diagnosis_age := as.numeric(diagnosis_date - dob) / 365.25]

# Load smoking information
dt_smoking <- fread("data/all-of-us/phenotypes/aou-smoking.txt")

dt_model <- merge(dt_model, dt_smoking, by = "person_id", all.x = TRUE)

dt_model <- dt_model[!is.na(smoking) & sex %in% c("Male", "Female") & !is.na(first_visit_date)]

fwrite(dt_model, "data/all-of-us/phenotypes/aou-time-to-event.txt", sep = "\t")
