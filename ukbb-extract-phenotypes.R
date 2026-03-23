system("pip3 install pandas", intern = TRUE)
install.packages(c("readr", "stringr", "dplyr", "glue", "reactable", "janitor", "remotes"))
remotes::install_github("laderast/xvhelper")

library(dplyr)
library(readr)
library(stringr)
library(glue)
library(reactable)
library(xvhelper)

# Assign dataset record-id
# record id can be either from a dataset or from a cohort
projectid <- "project-J62G6z8JjQG4YJ0JKzKXzKzZ"
rid <- "record-J63J398J8686PKpYBg39xQ9j"
# Assign joint dataset project-id:record-id
dataset <- glue::glue("{projectid}:{rid}")

cmd <- glue::glue("dx extract_dataset {dataset} -ddd")
system(cmd)
df_codings <- readr::read_csv("app69929_20260207035043.dataset.codings.csv")
df_entity <- readr::read_csv("app69929_20260207035043.dataset.entity_dictionary.csv")
df_data_dict <- readr::read_csv("app69929_20260207035043.dataset.data_dictionary.csv")

df_data_dict <- df_data_dict %>%
    relocate(name, title) %>%
    mutate(ent_field = glue::glue("{entity}.{name}"))

df_data_dict <- df_data_dict |>
  select(
    title, 
    name, 
    entity, 
    ent_field, 
    coding_name, 
    is_multi_select, 
    is_sparse_coding
  )

filtered_dict <- df_data_dict |> 
  filter(
      grepl("Diagnoses - ICD10", title) & entity == "participant"
  ) |> 
  arrange(title)

field_list <- paste(c("participant.p41270", "participant.eid"), collapse = ",")
cohort_template <- "dx extract_dataset {dataset} --fields {field_list} -o ukbb-diagnoses.csv"
cmd <- glue::glue(cohort_template)
system(cmd)

filtered_dict <- df_data_dict |> 
  filter(
    grepl("Sex", title) | 
      grepl("Smoking", title) | 
      grepl("Age at recruitment", title) |
      grepl("Diagnoses - ICD10", title) |
      grepl("Body mass index|bmi", title)
  ) |> 
  arrange(title)

field_list <- paste(c(
  "participant.eid",
  "participant.p21022", 
  "participant.p21001_i0",
  "participant.p23104_i0",
  "participant.p21001_i1",
  "participant.p23104_i1",
  "participant.p21001_i2",
  "participant.p23104_i2",
  "participant.p21001_i3",
  "participant.p23104_i3",
  "participant.p31",
  "participant.p20116_i0",
  "participant.p20116_i1",
  "participant.p20116_i2",
  "participant.p20116_i3"), 
  collapse = ","
)

cohort_template <- "dx extract_dataset {dataset} --fields {field_list} -o ukbb-demographics.csv"
cmd <- glue::glue(cohort_template)
system(cmd)
