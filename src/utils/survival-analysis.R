fit_one_mortality <- function(dt, gene, genes = GENES, covars = COVARS) {
  # Drop carriers of every *other* gene

  others <- setdiff(genes, gene)
  d <- dt[rowSums(dt[, ..others]) == 0]

  out <- data.table(
    gene = gene,
    median_carrier = NA_real_,
    median_noncarrier = NA_real_,
    HR = NA_real_,
    CI_lower = NA_real_,
    CI_upper = NA_real_,
    P = NA_real_
  )

  g <- d[[gene]]
  y <- d[["status"]]

  # Nothing to estimate if there are no cases or no carriers
  if (sum(g == 1 & y == 1) == 0) {
    return(out)
  }

  message(sprintf("Fitting survfit object for: %s", gene))
  f <- reformulate(termlabels = c(gene), response = "Surv(time, status)")
  fit <- do.call(survfit, list(formula = f, data = quote(d)))
  fit.table <- summary(fit)$table

  out[, `:=`(
    median_carrier = fit.table[sprintf("%s=1", gene), ]["median"],
    median_noncarrier = fit.table[sprintf("%s=0", gene), ]["median"]
  )]

  # Drop covariates with no variation (e.g. sex once breast is restricted to females)
  covars_use <- covars[vapply(
    covars,
    function(v) uniqueN(d[[v]]) > 1,
    logical(1)
  )]

  message(sprintf("Fitting coxphf object for: %s", gene))
  cox <- coxphf(
    reformulate(c(gene, covars_use), response = "Surv(time, status)"),
    data = d
  )

  out <- out[, `:=`(
    HR = exp(coef(cox)[gene]),
    CI_lower = cox$ci.lower[gene],
    CI_upper = cox$ci.upper[gene],
    P = cox$prob[gene]
  )]
}
