library(data.table)
library(coxphf)
library(survival)
source("src/utils/all-cause-mortality.R")

dt <- fread("data/pmbb/pmbb-time-to-event.txt")

dt <- dt[, age_last_followup := as.numeric(last_visit_date - dob) / 365.25]
dt <- dt[, age_first_visit := as.numeric(first_visit_date - dob) / 365.25]
dt <- dt[
    crep_highrisk_flag == 0 &
        !is.na(age_last_followup) &
        !is.na(dob) &
        !is.na(smoking)
]


dt[, `:=`(
    time = fifelse(!is.na(age_death), age_death, age_last_followup),
    status = as.integer(!is.na(age_death))
)]

GENES <- c("SDHA", "SDHB", "SDHC", "SDHD", "SDHAF2", "MAX")
COVARS <- c("sex", "smoking")

BIOBANK <- "pmbb"
out_dir <- "results/all-cause-mortality"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res_mortality <- rbindlist(Map(
    function(g) fit_one_mortality(dt, g),
    GENES
))

res_mortality[, biobank := BIOBANK]
setcolorder(res_mortality, "biobank")

fwrite(
    res_mortality,
    file.path(out_dir, sprintf("%s-all-cause-mortality.txt", BIOBANK)),
    sep = "\t",
    na = "NA"
)

dt[,
    status := fcase(
        !is.na(diagnosis_age) , 1 ,
        default = 0
    )
]
dt[,
    time := fcase(
        !is.na(diagnosis_age) , diagnosis_age ,
        default = age_last_followup
    )
]
dt <- dt[time > 0]

out_dir <- "results/time-to-cancer"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res_time <- rbindlist(Map(
    function(g) fit_one_mortality(dt, g),
    GENES
))

res_time[, biobank := BIOBANK]
setcolorder(res_time, "biobank")

fwrite(
    res_time,
    file.path(out_dir, sprintf("%s-time-to-cancer.txt", BIOBANK)),
    sep = "\t",
    na = "NA"
)
