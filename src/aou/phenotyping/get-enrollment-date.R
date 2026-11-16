library(bigrquery)
library(tidyverse)

download_data <- function(query) {
  tb <- bq_project_query(Sys.getenv('GOOGLE_PROJECT'), query = query, default_dataset = Sys.getenv('WORKSPACE_CDR'))
  bq_table_download(tb)
}

DATASET <- Sys.getenv('WORKSPACE_CDR')

df = download_data(str_glue("SELECT DISTINCT
person_id,
MIN(observation_date) AS primary_consent_date
FROM `{DATASET}.concept`
JOIN `{DATASET}.concept_ancestor` on concept_id = ancestor_concept_id
JOIN `{DATASET}.observation` on descendant_concept_id = observation_concept_id
WHERE concept_name = 'Consent PII' AND concept_class_id = 'Module'
GROUP BY 1"))

data.table::fwrite(df, "data/all-of-us/phenotypes/aou-enrollment-date.txt", sep = "\t")
