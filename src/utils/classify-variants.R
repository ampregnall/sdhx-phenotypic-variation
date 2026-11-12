classify_variants <- function(
  dt,
  allele_freq = 0.01,
  clinvar_annotation = "ClinVar.Sig",
  gnomad_frequency = "gnomAD.MAX_AF",
  autogvp_annotation = "AutoGVP",
  vep_annotation = "Variant.Consequence",
  exon_col = "EXON",
  am_annotation = "AlphaMissense",
  exon_sep = "/"
) {
  cols <- c(
    clinvar_annotation,
    gnomad_frequency,
    autogvp_annotation,
    vep_annotation,
    exon_col,
    am_annotation
  )
  missing <- setdiff(cols, names(dt))
  if (length(missing)) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }

  pathogenic_filter <- c(
    "Likely_pathogenic",
    "Pathogenic",
    "Pathogenic/Likely_pathogenic"
  )

  # Subset returns a new table, so no copy() needed
  dt <- dt[
    !grepl("benign", clin, ignore.case = TRUE),
    env = list(clin = clinvar_annotation)
  ]

  dt[,
    `:=`(
      .af = suppressWarnings(as.numeric(af)),
      .am = suppressWarnings(as.numeric(am))
    ),
    env = list(af = gnomad_frequency, am = am_annotation)
  ]

  dt[,
    Consequence := fcase(
      clin %chin% pathogenic_filter , "pLoF" ,
      agvp %chin% pathogenic_filter , "pLoF" ,
      default = NA_character_
    ),
    env = list(clin = clinvar_annotation, agvp = autogvp_annotation)
  ]

  dt <- dt[!is.na(Consequence) | is.na(.af) | .af < allele_freq]

  dt[,
    .last_exon := vapply(
      strsplit(as.character(ex), exon_sep, fixed = TRUE),
      function(x) length(x) >= 2 && x[1] == x[2],
      logical(1)
    ),
    env = list(ex = exon_col)
  ]

  dt[,
    Consequence := fcase(
      !is.na(Consequence)                                                   , Consequence          ,
      grepl("frameshift_variant|stop_gained|start_lost", vep) & !.last_exon , "pLoF"               ,
      grepl("splice_acceptor_variant|splice_donor_variant", vep)            , "pLoF"               ,
      grepl("missense", vep) & !is.na(.am) & .am >= 0.8                     , "Damaging Missense"  ,
      grepl("missense", vep)                                                , "Other Missense"     ,
      grepl("synonymous_variant", vep)                                      , "Synonymous Variant" ,
      default = NA_character_
    ),
    env = list(vep = vep_annotation)
  ]

  dt[, c(".af", ".am", ".last_exon") := NULL]
  dt[!is.na(Consequence)]
}
