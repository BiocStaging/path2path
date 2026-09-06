#' Stage 1: per-perturbation preranked GSEA.
#'
#' For each row m of `z` (a perturbation), ranks the measured genes by the
#' CRISPRi-negated statistic s = -z[m, ] and scores every pathway with
#' preranked GSEA. Returns the M x K NES matrix (and p-values for
#' method = "fgsea").
#'
#' @param z numeric matrix, perturbations x genes; rownames = perturbed gene
#'   symbols, colnames = measured gene symbols. NAs are dropped per row.
#' @param pathways named list of gene-symbol vectors.
#' @param method "fgsea" (fgseaSimple, per-perturbation permutations) or
#'   "pooled" (one pooled null per pathway across perturbations; no p-values).
#' @param exclude_target drop the perturbed gene from its own ranked list
#'   before scoring. This removes the tautological self-hit: the target is at
#'   the top of L_m by construction, giving every pathway that contains it a
#'   free hit. Default TRUE.
#' @param negate apply the CRISPRi loss-of-function sign flip (s = -z).
#'   Default TRUE; set FALSE only for gain-of-function perturbations.
#' @param b_prime pooled method only: random null draws per perturbation
#'   pooled into the per-pathway null (1 draw x thousands of perturbations
#'   is already a large pool).
#' @param nperm fgsea method only: permutations per ranked list for
#'   \code{fgsea::fgseaSimple}. The default 1000 estimates the NES
#'   accurately; downstream inference does not consume Stage 1 p-values.
#' @param min_size,max_size effective gene-set size window, applied to the
#'   intersection with each ranked list.
#' @param BPPARAM a \code{BiocParallel} parameter object controlling
#'   parallelisation over perturbations; default
#'   \code{BiocParallel::bpparam()}.
#' @details When ranked statistics contain ties, the function emits one
#'   summary warning for the stage rather than one warning per perturbation.
#' @return list(nes = M x K matrix, pval = M x K matrix or NULL)
#' @references
#' Korotkevich G, Sukhov V, Sergushichev A (2019). Fast gene set enrichment
#' analysis. bioRxiv. doi:10.1101/060012.
#' @examples
#' data(demo_z)
#' data(demo_pathways)
#' s1 <- stage1_nes(demo_z, demo_pathways, nperm = 100,
#'     min_size = 5, BPPARAM = BiocParallel::SerialParam())
#' dim(s1$nes)
#' @export
stage1_nes <- function(z, pathways,
                       method         = c("fgsea", "pooled"),
                       exclude_target = TRUE,
                       negate         = TRUE,
                       nperm          = 1000L,
                       min_size       = 15L,
                       max_size       = 500L,
                       b_prime        = 1L,
                       BPPARAM        = BiocParallel::bpparam()) {
    method <- match.arg(method)
    .check_matrix_input(z, "z")
    .check_pathways(pathways)
    .check_positive_integer(nperm, "nperm")
    .check_positive_integer(min_size, "min_size")
    .check_positive_integer(max_size, "max_size")
    .check_positive_integer(b_prime, "b_prime")
    if (min_size > max_size) stop("`min_size` cannot exceed `max_size`")
    if (!is.logical(exclude_target) || length(exclude_target) != 1L ||
        is.na(exclude_target)) stop("`exclude_target` must be TRUE or FALSE")
    if (!is.logical(negate) || length(negate) != 1L || is.na(negate))
        stop("`negate` must be TRUE or FALSE")
    pert_names    <- rownames(z)
    gene_names    <- colnames(z)
    pathway_names <- names(pathways)
    bp <- BiocParallel::SerialParam(progressbar = FALSE)

    row_stats <- function(i) {
        stats <- if (negate) -z[i, ] else z[i, ]
        names(stats) <- gene_names
        if (exclude_target) stats <- stats[names(stats) != pert_names[i]]
        sort(stats, decreasing = TRUE)   # sort() drops NAs
    }

    if (method == "fgsea") {
        run_one <- function(i) {
            stats_i <- row_stats(i)
            res <- .run_fgsea_ties_quiet(fgsea::fgseaSimple(
                pathways = pathways, stats = stats_i, nperm = nperm,
                minSize = min_size, maxSize = max_size, BPPARAM = bp
            ))
            list(nes = stats::setNames(res$NES, res$pathway),
                pval = stats::setNames(res$pval, res$pathway),
                tie_fraction = .tie_fraction(stats_i))
        }
        out <- BiocParallel::bplapply(
            seq_len(nrow(z)), run_one, BPPARAM = BPPARAM
        )
        .warn_ties_summary(
            vapply(out, `[[`, numeric(1), "tie_fraction"), "stage1_nes")
        nes  <- matrix(NA_real_, nrow(z), length(pathways),
            dimnames = list(pert_names, pathway_names))
        pval <- nes
        for (i in seq_along(out)) {
            nes[i,  names(out[[i]]$nes)]  <- out[[i]]$nes
            pval[i, names(out[[i]]$pval)] <- out[[i]]$pval
        }
        return(list(nes = nes, pval = pval))
    }

    ## ---- pooled ------------------------------------------------------------
    # One expected null |ES| per pathway, estimated from b_prime random draws
    # per perturbation pooled across all perturbations. Assumes the |s|^p
    # profile is approximately exchangeable across perturbations.
    per_pert <- function(m) {
        s_sorted <- row_stats(m)
        N <- length(s_sorted)
        gene_pos <- stats::setNames(seq_len(N), names(s_sorted))
        K <- length(pathways)
        es_obs <- rep(NA_real_, K)
        sizes <- integer(K)
        for (k in seq_len(K)) {
            sel <- gene_pos[pathways[[k]]]
            sel <- sort(unique(sel[!is.na(sel)]))
            sizes[k] <- length(sel)
            if (length(sel) >= min_size && length(sel) <= max_size)
                es_obs[k] <- fgsea::calcGseaStat(
                    s_sorted, selectedStats = sel, gseaParam = 1
                )
        }
        es_null <- matrix(NA_real_, b_prime, K)
        for (b in seq_len(b_prime)) for (k in seq_len(K)) {
            if (sizes[k] >= min_size && sizes[k] <= max_size) {
                sel <- sort(sample.int(N, sizes[k]))
                es_null[b, k] <- fgsea::calcGseaStat(
                    s_sorted, selectedStats = sel, gseaParam = 1
                )
            }
        }
        list(es_obs = es_obs, es_null = es_null,
            tie_fraction = .tie_fraction(s_sorted))
    }
    res <- BiocParallel::bplapply(seq_len(nrow(z)), per_pert, BPPARAM = BPPARAM)
    .warn_ties_summary(
        vapply(res, `[[`, numeric(1), "tie_fraction"),
        "stage1_nes[pooled]")
    es_obs <- do.call(rbind, lapply(res, `[[`, "es_obs"))
    dimnames(es_obs) <- list(pert_names, pathway_names)
    nes <- normalize_by_pooled_null(es_obs, lapply(res, `[[`, "es_null"))
    list(nes = nes, pval = NULL)
}

#' Normalize observed ES by sign-matched pooled null means (internal).
#'
#' If one sign side has fewer than `min_side` null draws, entries on that
#' side cannot be normalized; they are set to NA (with a warning) rather
#' than silently left on the raw-ES scale.
#' @noRd
normalize_by_pooled_null <- function(es_obs, es_null_list, min_side = 10L) {
    K <- ncol(es_obs)
    nes <- es_obs
    n_dropped <- 0L
    for (k in seq_len(K)) {
        all_null <- unlist(lapply(es_null_list, function(m) m[, k]))
        all_null <- all_null[!is.na(all_null)]
        pos <- all_null[all_null > 0]
        neg <- all_null[all_null < 0]
        idx_pos <- which(!is.na(es_obs[, k]) & es_obs[, k] > 0)
        idx_neg <- which(!is.na(es_obs[, k]) & es_obs[, k] < 0)
        if (length(pos) >= min_side) {
            nes[idx_pos, k] <- es_obs[idx_pos, k] / mean(pos)
        } else if (length(idx_pos)) {
            nes[idx_pos, k] <- NA_real_
            n_dropped <- n_dropped + length(idx_pos)
        }
        if (length(neg) >= min_side) {
            nes[idx_neg, k] <- es_obs[idx_neg, k] / abs(mean(neg))
        } else if (length(idx_neg)) {
            nes[idx_neg, k] <- NA_real_
            n_dropped <- n_dropped + length(idx_neg)
        }
    }
    if (n_dropped)
        warning("normalize_by_pooled_null: ", n_dropped,
            " entries had too few same-sign null draws (<", min_side,
            ") and were set to NA")
    nes
}
