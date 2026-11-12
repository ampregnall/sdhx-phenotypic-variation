fit_one <- function(dt, pheno, gene, genes = GENES, covars = COVARS) {
    # Drop carriers of every *other* gene
    others <- setdiff(genes, gene)
    d <- dt[rowSums(dt[, ..others]) == 0]

    d <- d[complete.cases(d[, c(pheno, gene, covars), with = FALSE])]
    d[, (pheno) := as.integer(get(pheno))]

    y <- d[[pheno]]
    g <- d[[gene]]
    out <- data.table(
        phenotype = pheno,
        gene = gene,
        n_case_carrier = sum(y == 1 & g == 1),
        n_case = sum(y == 1),
        n_ctrl_carrier = sum(y == 0 & g == 1),
        n_ctrl = sum(y == 0),
        OR = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        P = NA_real_
    )

    # Nothing to estimate if there are no cases or no carriers
    if (out$n_case == 0 || sum(g == 1) == 0) {
        return(out)
    }

    # Drop covariates with no variation (e.g. sex once breast is restricted to females)
    covars_use <- covars[vapply(
        covars,
        function(v) uniqueN(d[[v]]) > 1,
        logical(1)
    )]

    fit <- tryCatch(
        logistf(reformulate(c(gene, covars_use), response = pheno), data = d),
        error = function(e) {
            message(sprintf("%s / %s: %s", pheno, gene, conditionMessage(e)))
            NULL
        }
    )
    if (is.null(fit)) {
        return(out)
    }

    i <- match(gene, names(coef(fit))) # assumes gene columns are numeric 0/1
    out[, `:=`(
        OR = exp(coef(fit)[i]),
        CI_lower = exp(fit$ci.lower[i]),
        CI_upper = exp(fit$ci.upper[i]),
        P = fit$prob[i]
    )]
    out
}
