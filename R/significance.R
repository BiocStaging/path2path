#' Internal: the single definition of "significant" used across the API.
#'
#' Cutoff semantics (shared by plot_map dots, select_active, and
#' order_by_activity): each non-NULL cutoff is applied; when both are given
#' they must BOTH hold (conjunction). NULL disables a cutoff; at least one
#' must be enabled. With `calibrated` supplied, cutoffs are applied to
#' calibrated$pval_gpd / calibrated$padj_gpd (Route B); otherwise to
#' result$pval / result$padj (Route A).
#' @noRd
.sig_matrix <- function(result, calibrated = NULL,
                        padj_cutoff = 0.05, pval_cutoff = NULL) {
    if (is.null(padj_cutoff) && is.null(pval_cutoff))
        stop("enable at least one of `padj_cutoff` / `pval_cutoff`")
    if (!is.null(padj_cutoff))
        .check_probability(padj_cutoff, "padj_cutoff")
    if (!is.null(pval_cutoff))
        .check_probability(pval_cutoff, "pval_cutoff")
    if (is.null(calibrated)) {
        P <- result$pval
        Q <- result$padj
        if (is.null(P))
            stop(
                "result has no $pval ",
                "(pooled-engine results carry no p-values)"
            )
    } else {
        if (is.null(calibrated$padj_gpd))
            stop(
                "`calibrated` has no $padj_gpd; recompute it with the current ",
                "calibrate_pvalues()"
            )
        P <- calibrated$pval_gpd
        Q <- calibrated$padj_gpd
    }
    S <- !is.na(P)
    if (!is.null(padj_cutoff)) S <- S & !is.na(Q) & Q < padj_cutoff
    if (!is.null(pval_cutoff)) S <- S & P < pval_cutoff
    S
}

#' Pathways with at least `min_sig` significant entries, as regulators and
#' as readouts.
#'
#' Counts significant entries per row (outgoing, pathway as regulator) and
#' per column (incoming, pathway as readout) separately -- in practice
#' out-degrees and in-degrees differ by an order of magnitude, so a joint
#' count would drown the readout side. The result plugs directly into
#' plot_map(): `plot_map(res, rows = sel$rows, cols = sel$cols)`.
#'
#' @param min_sig integer(1): minimum number of significant entries.
#' @param include_diag logical(1): count the self-interaction diagonal
#'   (default FALSE,
#'   matching the paper's off-diagonal focus; only applies to square maps).
#' @param result a list as returned by \code{path2path()}.
#' @param calibrated `NULL` or a list as returned by
#'   \code{calibrate_pvalues()} (switches
#'   the cutoffs to Route B values).
#' @param padj_cutoff,pval_cutoff `NULL` or numeric(1) each: significance
#'   cutoffs; non-NULL cutoffs
#'   conjoin, NULL disables one (at least one must be enabled).
#' @return list(rows = character, cols = character)
#' @examples
#' data(demo_z)
#' data(demo_pathways)
#' res <- path2path(demo_z, demo_pathways,
#'     stage1 = list(nperm = 100, min_size = 5),
#'     stage2 = list(min_pu = 5),
#'     BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
#' select_active(res, padj_cutoff = 0.2)
#' @export
select_active <- function(result, calibrated = NULL,
                          padj_cutoff = 0.05, pval_cutoff = NULL,
                          min_sig = 1L, include_diag = FALSE) {
    .check_positive_integer(min_sig, "min_sig")
    S <- .sig_matrix(result, calibrated, padj_cutoff, pval_cutoff)
    if (!include_diag && nrow(S) == ncol(S) &&
        identical(rownames(S), colnames(S))) diag(S) <- FALSE
    list(rows = rownames(S)[rowSums(S) >= min_sig],
        cols = colnames(S)[colSums(S) >= min_sig])
}

#' Pathway names ordered by number of significant entries (descending).
#'
#' `margin = "rows"` orders by outgoing (regulator) counts, `"cols"` by
#' incoming (readout) counts. Returns ALL pathway names of that margin
#' (zero-count ones last), so the result is a complete ordering usable as
#' plot_map(order = ...).
#' @inheritParams select_active
#' @param margin character(1): order by outgoing (\code{"rows"}) or incoming
#'   (\code{"cols"}) significant counts.
#' @param include_diag logical(1): count the self-interaction diagonal.
#' @return character vector: a complete ordering of that margin's names.
#' @examples
#' data(demo_z)
#' data(demo_pathways)
#' res <- path2path(demo_z, demo_pathways,
#'     stage1 = list(nperm = 100, min_size = 5),
#'     stage2 = list(min_pu = 5),
#'     BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
#' order_by_activity(res, padj_cutoff = 0.2)[1:3]
#' @export
order_by_activity <- function(result, calibrated = NULL,
                              padj_cutoff = 0.05, pval_cutoff = NULL,
                              margin = c("rows", "cols"),
                              include_diag = FALSE) {
    margin <- match.arg(margin)
    S <- .sig_matrix(result, calibrated, padj_cutoff, pval_cutoff)
    if (!include_diag && nrow(S) == ncol(S) &&
        identical(rownames(S), colnames(S))) diag(S) <- FALSE
    counts <- if (margin == "rows") rowSums(S) else colSums(S)
    names(sort(counts, decreasing = TRUE))
}
