#' Stage 2: per-pathway-column GSEA in the perturbation universe.
#'
#' For each column k of the Stage 1 NES matrix, ranks the perturbations by
#' NES(., k) and asks, for every pathway j, whether the perturbations of
#' S_j members cluster at either end. Pathways are restricted to the
#' perturbation universe (rownames of `nes`) and must retain at least
#' `min_pu` members.
#'
#' @param nes numeric matrix, M x K Stage 1 NES values (rownames =
#'   perturbed gene symbols, colnames = pathway names).
#' @param pathways named list of character vectors: the full pathway
#'   collection (gene symbols); restricted internally.
#' @param method character(1): "fgsea" (fgseaMultilevel, p-values) or
#'   "pooled".
#' @param min_pu integer(1): minimum members of a set in the perturbation
#'   universe
#'   (rownames of `nes`) for the set to be tested at all.
#' @param min_col_size integer(1): minimum members of a set actually
#'   present in one
#'   column's non-NA ranking for that entry to be scored. The universe
#'   filter `min_pu` is global; within a single column, NA entries can
#'   reduce a set's effective size below `min_pu`, and `min_col_size`
#'   bounds it from below (default 2, ruling out degenerate
#'   single-member scores). Raise it (e.g. to `min_pu`) for a stricter
#'   per-entry guarantee.
#' @param b_prime integer(1), pooled method only: random null draws per
#'   column pooled
#'   into the per-pathway-row null.
#' @param score_type character(1): "std" (default GSEA two-sided score),
#'   or "both" to run
#'   scoreType "pos" and "neg" separately and return both extremes
#'   (M_pos / M_neg), which resolves bidirectionally-regulated pairs that
#'   cancel under "std".
#' @param max_size integer(1): maximum effective set size within a column's
#'   ranking.
#' @param min_valid integer(1): minimum non-NA entries a column needs to be
#'   ranked at all.
#' @param eps numeric(1): \code{fgsea::fgseaMultilevel} accuracy floor
#'   (0 = exact
#'   estimation of arbitrarily small p-values).
#' @param BPPARAM a \code{BiocParallel} parameter object controlling
#'   parallelisation over readout columns; default
#'   \code{BiocParallel::bpparam()}.
#' @details When ranked statistics contain ties, the function emits one
#'   summary warning for the stage rather than one warning per readout column.
#' @return list with K x K matrices: M, pval (std) or
#'   M_pos, M_neg, pval_pos, pval_neg (both). Rows cover all pathways;
#'   pathways failing the min_pu filter are all-NA rows.
#' @examples
#' data(demo_z)
#' data(demo_pathways)
#' s1 <- stage1_nes(demo_z, demo_pathways, nperm = 100, min_size = 5,
#'     BPPARAM = BiocParallel::SerialParam())
#' s2 <- stage2_matrix(s1$nes, demo_pathways, min_pu = 5,
#'     BPPARAM = BiocParallel::SerialParam())
#' s2$M[1:3, 1:3]
#' @export
stage2_matrix <- function(nes, pathways,
                          method     = c("fgsea", "pooled"),
                          score_type = c("std", "both"),
                          min_pu     = 10L,
                          min_col_size = 2L,
                          max_size   = 500L,
                          min_valid  = 20L,
                          eps        = 0,
                          b_prime    = 1L,
                          BPPARAM    = BiocParallel::bpparam()) {
    method     <- match.arg(method)
    score_type <- match.arg(score_type)
    .check_matrix_input(nes, "nes")
    .check_pathways(pathways)
    .check_positive_integer(min_pu, "min_pu")
    .check_positive_integer(min_col_size, "min_col_size")
    .check_positive_integer(max_size, "max_size")
    .check_positive_integer(min_valid, "min_valid")
    .check_positive_integer(b_prime, "b_prime")
    if (method == "pooled" && score_type != "std")
        stop("score_type='both' is only implemented for method='fgsea'")
    if (min_col_size > max_size)
        stop("`min_col_size` cannot exceed `max_size`")
    if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps < 0)
        stop("`eps` must be one non-negative number")
    perturbed <- rownames(nes)
    # Rows of the output = the gene sets being tested (restricted to the
    # perturbation universe); columns = the columns of `nes`. In the standard
    # pipeline these coincide, but they are deliberately decoupled so that
    # arbitrary set families (e.g. negative-control pseudo-pathways) can be
    # scored against an existing NES matrix.
    set_names <- names(pathways)
    col_names <- colnames(nes)
    pathways_pu <- lapply(pathways, intersect, perturbed)
    pathways_pu <- pathways_pu[lengths(pathways_pu) >= min_pu]
    bp <- BiocParallel::SerialParam(progressbar = FALSE)

    empty_KK <- function() matrix(NA_real_, length(pathways), length(col_names),
        dimnames = list(set_names, col_names))

    if (!length(pathways_pu)) {
        empty <- empty_KK()
        if (score_type == "std") return(list(M = empty, pval = empty))
        return(list(M_pos = empty, pval_pos = empty,
            M_neg = empty, pval_neg = empty))
    }

    col_stats <- function(k) {
        stats_k <- nes[, k]
        stats_k <- stats_k[!is.na(stats_k)]
        if (length(stats_k) < min_valid) return(NULL)
        sort(stats_k, decreasing = TRUE)
    }

    if (method == "fgsea") {
        run_col <- function(k_idx, st) {
            stats_k <- col_stats(col_names[k_idx])
            if (is.null(stats_k)) return(NULL)
            res <- .run_fgsea_ties_quiet(.run_fgsea_multilevel(
                fgsea::fgseaMultilevel(
                    pathways = pathways_pu, stats = stats_k, eps = eps,
                    minSize = min_col_size, maxSize = max_size,
                    scoreType = st, BPPARAM = bp
                )
            ))
            list(nes = stats::setNames(res$NES, res$pathway),
                pval = stats::setNames(res$pval, res$pathway),
                tie_fraction = .tie_fraction(stats_k))
        }
        fill <- function(st) {
            out <- BiocParallel::bplapply(seq_along(col_names), run_col,
                st = st, BPPARAM = BPPARAM)
            M <- empty_KK()
            P <- empty_KK()
            for (k_idx in seq_along(col_names)) {
                v <- out[[k_idx]]
                if (!is.null(v)) {
                    M[names(v$nes),  col_names[k_idx]] <- v$nes
                    P[names(v$pval), col_names[k_idx]] <- v$pval
                }
            }
            tie_fraction <- vapply(out, function(v) {
                if (is.null(v)) 0 else v$tie_fraction
            }, numeric(1))
            list(M = M, pval = P, tie_fraction = tie_fraction)
        }
        if (score_type == "std") {
            r <- fill("std")
            .warn_ties_summary(r$tie_fraction, "stage2_matrix")
            return(list(M = r$M, pval = r$pval))
        } else {
            rp <- fill("pos")
            rn <- fill("neg")
            .warn_ties_summary(rp$tie_fraction, "stage2_matrix")
            return(list(M_pos = rp$M, pval_pos = rp$pval,
                M_neg = rn$M, pval_neg = rn$pval))
        }
    }

    ## ---- pooled ------------------------------------------------------------
    per_col <- function(k_idx) {
        stats_k <- col_stats(col_names[k_idx])
        K_paths <- length(pathways_pu)
        if (is.null(stats_k))
            return(list(es_obs = rep(NA_real_, K_paths),
                es_null = matrix(NA_real_, b_prime, K_paths),
                tie_fraction = 0))
        N_k <- length(stats_k)
        pos <- stats::setNames(seq_len(N_k), names(stats_k))
        es_obs <- rep(NA_real_, K_paths)
        sizes <- integer(K_paths)
        for (j in seq_len(K_paths)) {
            sel <- pos[pathways_pu[[j]]]
            sel <- sort(unique(sel[!is.na(sel)]))
            sizes[j] <- length(sel)
            if (length(sel) >= min_col_size && length(sel) <= max_size)
                es_obs[j] <- fgsea::calcGseaStat(
                    stats_k, selectedStats = sel, gseaParam = 1
                )
        }
        es_null <- matrix(NA_real_, b_prime, K_paths)
        for (b in seq_len(b_prime)) for (j in seq_len(K_paths)) {
            if (sizes[j] >= min_col_size && sizes[j] <= max_size) {
                sel <- sort(sample.int(N_k, sizes[j]))
                es_null[b, j] <- fgsea::calcGseaStat(
                    stats_k, selectedStats = sel, gseaParam = 1
                )
            }
        }
        list(es_obs = es_obs, es_null = es_null,
            tie_fraction = .tie_fraction(stats_k))
    }
    res <- BiocParallel::bplapply(
        seq_along(col_names), per_col, BPPARAM = BPPARAM
    )
    .warn_ties_summary(
        vapply(res, `[[`, numeric(1), "tie_fraction"),
        "stage2_matrix[pooled]")
    es_obs <- matrix(NA_real_, length(pathways_pu), length(col_names),
        dimnames = list(names(pathways_pu), col_names))
    for (k in seq_along(col_names)) es_obs[, k] <- res[[k]]$es_obs
    # transpose orientation: pool nulls per pathway ROW j across columns
    null_by_row <- lapply(seq_along(pathways_pu), function(j)
        do.call(cbind, lapply(res, function(r) r$es_null[, j, drop = FALSE])))
    M_partial <- es_obs
    n_dropped <- 0L
    for (j in seq_along(pathways_pu)) {
        all_null <- as.vector(null_by_row[[j]])
        all_null <- all_null[!is.na(all_null)]
        pos <- all_null[all_null > 0]
        neg <- all_null[all_null < 0]
        idx_pos <- which(!is.na(es_obs[j, ]) & es_obs[j, ] > 0)
        idx_neg <- which(!is.na(es_obs[j, ]) & es_obs[j, ] < 0)
        if (length(pos) >= 10L) {
            M_partial[j, idx_pos] <- es_obs[j, idx_pos] / mean(pos)
        } else if (length(idx_pos)) {
            M_partial[j, idx_pos] <- NA_real_
            n_dropped <- n_dropped + length(idx_pos)
        }
        if (length(neg) >= 10L) {
            M_partial[j, idx_neg] <- es_obs[j, idx_neg] / abs(mean(neg))
        } else if (length(idx_neg)) {
            M_partial[j, idx_neg] <- NA_real_
            n_dropped <- n_dropped + length(idx_neg)
        }
    }
    if (n_dropped)
        warning("stage2_matrix[pooled]: ", n_dropped,
            " entries had too few same-sign null draws and were set to NA")
    M <- empty_KK()
    M[rownames(M_partial), colnames(M_partial)] <- M_partial
    list(M = M, pval = NULL)
}
