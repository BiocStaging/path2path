#' Whole-matrix permutation: one global row and one global
#' column permutation.
#'
#' Preserves the marginal value distribution and within-row gene-gene
#' correlation (the same column permutation is applied to every row) while
#' breaking pathway membership, because the original labels stay attached to
#' the permuted positions.
#'
#' @param z perturbations x genes numeric matrix.
#' @return the permuted matrix (labels retained in place, so every label
#'   sits atop another row/column's data).
#' @examples
#' data(demo_z)
#' set.seed(1)
#' zp <- permute_input(demo_z)
#' identical(dimnames(zp), dimnames(demo_z))
#' @export
permute_input <- function(z) {
    .check_matrix_input(z, "z")
    Z <- z[sample.int(nrow(z)), sample.int(ncol(z))]
    dimnames(Z) <- dimnames(z)
    Z
}

#' Whole-pipeline permutation null calibration.
#'
#' Runs the same path2path() pipeline on B permuted matrices and
#' pools |M'| entries into a sorted null pool. Each null M' is saved to
#' out_dir for reuse (per-pathway regrouping, tail fits).
#'
#' @param ... additional path2path() arguments; these override entries of
#'   `config` with the same name.
#' @param z,pathways the inputs of the observed run: `z` a numeric matrix
#'   (perturbations x genes, rownames = perturbed gene symbols) and
#'   `pathways` a named list of character vectors of gene symbols.
#' @param B integer(1): number of outer permutations.
#' @param out_dir `NULL` or character(1): optional directory persisting
#'   each null map
#'   (M_null_bNNN.rds), per-map BH discovery counts, and each null map's
#'   sorted nominal p-values; completed maps found there are reloaded
#'   instead of recomputed (resume). Use a separate directory for each input
#'   dataset, pathway collection, and pipeline configuration.
#' @param config `NULL` or a list, a `path2path()` result's `$config`
#'   component, so null runs repeat the observed configuration exactly.
#' @param seed_base `NULL` or integer(1): if non-NULL, uses `seed_base + b`
#'   for the input permutation in run b. For exact reproducibility of the
#'   GSEA calculations, also use seeded `BPPARAM` objects in `config` or
#'   `...`.
#' @param verbose logical(1): report per-permutation progress messages
#'   (default \code{FALSE}).
#' @return a list with two components: `null_M`, a length-`B` list of
#'   pathway x pathway numeric matrices (one null map per permutation), and
#'   `pool_abs`, a sorted numeric vector pooling the absolute values of all
#'   non-missing null-map entries.
#' @examples
#' data(demo_z)
#' data(demo_pathways)
#' cfg <- list(stage1 = list(nperm = 100, min_size = 5),
#'     stage2 = list(min_pu = 5),
#'     BPPARAM = BiocParallel::SerialParam())
#' pool <- build_null_pool(demo_z, demo_pathways, B = 2, seed_base = 1,
#'     stage1 = cfg$stage1, stage2 = cfg$stage2,
#'     BPPARAM = cfg$BPPARAM)
#' length(pool$pool_abs)
#' @export
build_null_pool <- function(z, pathways, B = 10L, out_dir = NULL,
                            config = NULL, seed_base = NULL,
                            verbose = FALSE, ...) {
    .check_matrix_input(z, "z")
    .check_pathways(pathways)
    .check_positive_integer(B, "B")
    if (!is.null(seed_base)) .check_integer(seed_base, "seed_base")
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
        stop("`verbose` must be TRUE or FALSE")
    if (!is.null(out_dir) &&
        (!is.character(out_dir) || length(out_dir) != 1L || is.na(out_dir))) {
        stop("`out_dir` must be NULL or one non-missing path")
    }
    args <- if (is.null(config)) list() else config
    dots <- list(...)
    if (length(dots)) args <- utils::modifyList(args, dots)
    if (is.null(args$verbose)) args$verbose <- FALSE
    if (!is.null(out_dir))
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    bh_file <- if (!is.null(out_dir))
        file.path(out_dir, "null_bh_counts.csv") else NULL
    alphas <- c(0.05, 0.02, 0.01, 0.005, 0.001)
    null_M <- vector("list", B)
    for (b in seq_len(B)) {
        m_file <- if (!is.null(out_dir))
            file.path(out_dir, sprintf("M_null_b%03d.rds", b)) else NULL
        ## A per-run seed makes the input permutation reproducible in isolation.
        if (!is.null(m_file) && file.exists(m_file)) {
            if (verbose)
                message(sprintf("[%s] null permutation %d / %d (cached)",
                    format(Sys.time(), "%H:%M:%S"), b, B))
            cached <- readRDS(m_file)
            expected <- c(length(pathways), length(pathways))
            if (!is.matrix(cached) || !identical(dim(cached), expected) ||
                !identical(rownames(cached), names(pathways)) ||
                !identical(colnames(cached), names(pathways))) {
                stop(
                    "cached null map does not match the pathway collection: ",
                    m_file)
            }
            null_M[[b]] <- cached
            next
        }
        if (verbose)
            message(sprintf("[%s] null permutation %d / %d",
                format(Sys.time(), "%H:%M:%S"), b, B))
        z_null <- if (is.null(seed_base)) {
            permute_input(z)
        } else {
            withr::with_seed(seed_base + b, permute_input(z))
        }
        res <- do.call(path2path,
            c(list(z = z_null, pathways = pathways), args))
        if (is.null(res$M))
            stop("permutation calibration requires `score_type = 'std'`")
        null_M[[b]] <- res$M
        if (!is.null(m_file)) saveRDS(res$M, m_file)
        ## Per-null-map Route A BH discovery counts: how often BH fires on
        ## all-null data (its family-wise entry behaviour).
        if (!is.null(bh_file) && !is.null(res$padj)) {
            cnts <- vapply(
                alphas, function(a) sum(res$padj < a, na.rm = TRUE), 0L
            )
            row <- data.frame(b = b, populated = sum(!is.na(res$M)),
                t(stats::setNames(cnts, paste0("bh_", alphas))))
            utils::write.table(
                row, bh_file, sep = ",", append = file.exists(bh_file),
                col.names = !file.exists(bh_file), row.names = FALSE)
        }
        ## Full sorted null nominal-p vector: lets E[V] for ANY p-threshold rule
        ## (in particular the observed map's BH cut, unknown in advance) be
        ## computed exactly afterwards. A few MB gzipped per map.
        if (!is.null(out_dir) && !is.null(res$pval)) {
            pv <- sort(res$pval[!is.na(res$pval)])
            data.table::fwrite(data.table::data.table(p = pv),
                file.path(
                    out_dir, sprintf("pvals_null_b%03d.csv.gz", b)
            ))
        }
    }
    pool <- sort(unlist(lapply(null_M, function(M) abs(M)[!is.na(M)]),
        use.names = FALSE))
    list(null_M = null_M, pool_abs = pool)
}

#' Pooled empirical p-values and BH-FDR from a sorted |null| pool.
#'
#' @noRd
pooled_pvalues <- function(M_obs, pool_abs_sorted) {
    L <- length(pool_abs_sorted)
    exceed <- L - findInterval(abs(M_obs), pool_abs_sorted, left.open = TRUE)
    pval <- (1 + exceed) / (L + 1)
    dim(pval) <- dim(M_obs)
    dimnames(pval) <- dimnames(M_obs)
    pval[is.na(M_obs)] <- NA
    list(pval_pool = pval, fdr_pool = bh_populated(pval))
}

#' Generalized Pareto survival function above a threshold.
#' @noRd
.gpd_survival <- function(x, u, xi, beta, tail_fraction = 1) {
    if (abs(xi) < 1e-8) {
        return(tail_fraction * exp(-(x - u) / beta))
    }
    tail_fraction * pmax(1 + xi * (x - u) / beta, 0)^(-1 / xi)
}

#' Generalized Pareto tail fit (internal): for p-value extrapolation
#' beyond the pool resolution (Knijnenburg et al. 2009). Fits exceedances over
#' the `q_thresh` quantile of the null pool by maximum likelihood and returns a
#' function mapping |M| to a tail-extrapolated p-value.
#' @noRd
gpd_tail_fit <- function(pool_abs_sorted, q_thresh = 0.99) {
    u <- stats::quantile(pool_abs_sorted, q_thresh, names = FALSE)
    exc <- pool_abs_sorted[pool_abs_sorted > u] - u
    n_exc <- length(exc)
    n_tot <- length(pool_abs_sorted)
    if (n_exc < 30L) stop("Too few exceedances for a stable GPD fit")
    if (!is.finite(stats::sd(exc)) || stats::sd(exc) <= 0)
        stop("GPD exceedances have no usable variation")
    nll <- function(par) {
        xi <- par[1]
        beta <- exp(par[2])
        y <- 1 + xi * exc / beta
        if (any(y <= 0)) return(1e10)
        if (abs(xi) < 1e-8) {
            return(sum(log(beta) + exc / beta))
        }
        sum(log(beta) + (1 / xi + 1) * log(y))
    }
    # Nelder-Mead at default tolerances can stop short of the optimum on this
    # likelihood (verified against an independent implementation); tighten and
    # polish with BFGS from the NM solution
    fit <- stats::optim(
        c(0.1, log(stats::sd(exc))), nll, method = "Nelder-Mead",
        control = list(reltol = 1e-12, maxit = 5000))
    polished <- tryCatch(
        stats::optim(fit$par, nll, method = "BFGS",
            control = list(reltol = 1e-12, maxit = 500)),
        error = function(e) NULL)
    if (!is.null(polished) && polished$convergence == 0L &&
        is.finite(polished$value) && polished$value <= fit$value)
        fit <- polished
    if (fit$convergence != 0L || !is.finite(fit$value))
        stop("GPD optimization did not converge")
    xi <- fit$par[1]
    beta <- exp(fit$par[2])
    # For xi < 0 the fitted GPD has a finite support endpoint; beyond it the
    # fitted null assigns probability zero and no finite p-value is estimable.
    support_end <- if (xi < 0) u + beta / abs(xi) else Inf
    p_of <- function(x) {
        # empirical branch below u counts >= x, matching pooled_pvalues();
        # Above u, use the GPD estimate scaled by the observed tail fraction.
        tail <- .gpd_survival(x, u, xi, beta, n_exc / n_tot)
        p <- ifelse(
            x <= u,
            (1 + (n_tot - findInterval(
                x, pool_abs_sorted, left.open = TRUE
            ))) / (n_tot + 1),
            tail
        )
        beyond <- !is.na(x) & x >= support_end
        if (any(beyond)) {
            warning(
                "gpd_tail_fit: ", sum(beyond),
                " value(s) lie at or beyond the ",
                "fitted GPD support endpoint (", signif(support_end, 4), "); ",
                "their p-values are not estimable from this fit and are ",
                "returned as NA. Consider a larger B or a higher q_thresh.")
            p[beyond] <- NA_real_
        }
        p
    }
    list(u = u, xi = xi, beta = beta, n_exc = n_exc, n_tot = n_tot,
        support_end = support_end, pval = p_of)
}

#' One-call calibrated p-values for an existing map.
#'
#' Runs (or reuses) the outer-permutation calibration and returns the map's
#' entries with pooled-empirical and GPD-calibrated p-values attached.
#' Expensive: B full pipeline runs. `pool` may be supplied from a previous
#' build_null_pool() to skip recomputation.
#'
#' @param result a list as returned by path2path() (uses result$M, and reuses
#'   result$config so the null runs repeat the observed run's
#'   configuration exactly; there is deliberately no way to override it
#'   here, since mismatched null runs would invalidate the calibration.
#'   For full control, call build_null_pool() directly and pass its pool
#'   via `pool`).
#' @param z,pathways the inputs that produced `result`: `z` a numeric
#'   matrix (perturbations x genes) and `pathways` a named list of
#'   character vectors of gene symbols. Both may be `NULL` when a
#'   precomputed `pool` is supplied.
#' @param B integer(1): number of outer permutations (ignored if `pool`
#'   given).
#' @param pool `NULL` or numeric vector: optional sorted |M'| pool from
#'   build_null_pool()$pool_abs.
#' @param null_dir `NULL` or character(1), optional directory: each null
#'   map is saved there as
#'   M_null_bNNN.rds (persistence only; forwarded to build_null_pool's
#'   out_dir).
#' @return list with K x K matrices `pval_gpd` (the calibrated p-values:
#'   identical to `pval_pool` up to the pool's `q_thresh` quantile, GPD
#'   tail extrapolation beyond it when an adequate fit is available),
#'   `padj_gpd` (single BH family over the
#'   populated entries of `pval_gpd`, the calibrated analogue of
#'   result$padj), and `pval_pool` (pooled empirical p-values, kept for
#'   diagnostics); plus `gpd` (fit parameters) and `pool_abs` (the pool,
#'   reusable via `pool =`).
#' @param q_thresh numeric(1) in (0, 1): GPD threshold quantile of the
#'   pool (default 0.99); a
#'   fit-adequacy guard lowers it stepwise (0.975, 0.95, 0.90) if the fitted
#'   tail's finite endpoint falls below the largest observed |M|.
#' @param seed_base `NULL` or integer(1): forwarded to
#'   \code{build_null_pool()}.
#' @param verbose logical(1): forwarded to \code{build_null_pool()}
#'   (per-permutation progress messages; default \code{FALSE}).
#' @references
#' Knijnenburg TA, Wessels LFA, Reinders MJT, Shmulevich I (2009). Fewer
#' permutations, more accurate p-values. Bioinformatics 25:i161-i168.
#' doi:10.1093/bioinformatics/btp211.
#' @examples
#' data(demo_z)
#' data(demo_pathways)
#' res <- path2path(demo_z, demo_pathways,
#'     stage1 = list(nperm = 100, min_size = 5),
#'     stage2 = list(min_pu = 5),
#'     BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
#' cal <- calibrate_pvalues(res, demo_z, demo_pathways, B = 2, seed_base = 1)
#' summary(as.vector(cal$pval_gpd))
#' @export
calibrate_pvalues <- function(result, z = NULL, pathways = NULL, B = 50L,
                              pool = NULL, q_thresh = 0.99, null_dir = NULL,
                              seed_base = NULL, verbose = FALSE) {
    if (is.null(result$M) || !is.matrix(result$M)) {
        stop("`result` must be a standard path2path result containing `$M`")
    }
    .check_probability(q_thresh, "q_thresh", open = TRUE)
    ## Fit-adequacy guard: a fitted GPD whose finite support endpoint lies at
    ## or below the largest observed |M| contradicts the data (it assigns the
    ## strongest entries probability zero). In that case the threshold is
    ## lowered stepwise -- more exceedances, stabler fit -- until adequate.
    q_ladder <- unique(c(q_thresh, c(0.975, 0.95, 0.90)[
        c(0.975, 0.95, 0.90) < q_thresh
    ]))
    if (is.null(pool)) {
        stopifnot(!is.null(z), !is.null(pathways))
        config <- result$config
        if (!is.null(config) &&
            (length(config$stage1) || length(config$stage2) ||
                !identical(config$jaccard_tau, 0.5)))
            message(
                "calibrate_pvalues: reusing the non-default configuration ",
                "in `result` for the null runs"
            )
        pool <- build_null_pool(z, pathways, B = B, config = config,
            out_dir = null_dir, seed_base = seed_base,
            verbose = verbose)$pool_abs
    }
    if (!is.numeric(pool) || !length(pool) || anyNA(pool) ||
        any(!is.finite(pool)) || any(pool < 0)) {
        stop(
            "`pool` must be a non-empty finite numeric vector ",
            "of absolute scores"
        )
    }
    pool <- sort(pool)
    if (!any(!is.na(result$M))) stop("`result$M` contains no testable entries")
    pe <- pooled_pvalues(result$M, pool)
    max_obs <- max(abs(result$M), na.rm = TRUE)
    fit <- NULL
    adequate <- FALSE
    for (q in q_ladder) {
        candidate <- tryCatch(
            gpd_tail_fit(pool, q_thresh = q),
            error = function(e) NULL
        )
        if (is.null(candidate)) next
        if (!is.finite(candidate$support_end) ||
            candidate$support_end > max_obs) {
            fit <- candidate
            adequate <- TRUE
            break
        }
        message(
            "calibrate_pvalues: GPD fit at q=", sprintf("%.3f", q),
            " has endpoint ", sprintf("%.2f", candidate$support_end),
            " <= max|M| ", sprintf("%.2f", max_obs),
            "; lowering threshold"
        )
    }
    if (!adequate) {
        # pool too small for a stable tail fit at any ladder rung: degrade
        # gracefully to the pooled empirical p-values (their resolution floor
        # is 1/(L+1); increase B for tail extrapolation)
        message("calibrate_pvalues: no adequate GPD fit; returning pooled ",
            "empirical p-values only (increase B for tail extrapolation)")
        pg <- pe$pval_pool
    } else {
        fit$q_used <- q
        pg <- fit$pval(abs(result$M))
        dim(pg) <- dim(result$M)
        dimnames(pg) <- dimnames(result$M)
        pg[is.na(result$M)] <- NA
    }
    list(pval_pool = pe$pval_pool, pval_gpd = pg,
        padj_gpd = bh_populated(pg),
        gpd = if (is.null(fit)) NULL else
            fit[c("u", "xi", "beta", "n_exc", "n_tot", "support_end",
                "q_used")],
        pool_abs = pool)
}
