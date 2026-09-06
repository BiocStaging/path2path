#' Full two-stage GSEA pipeline: z matrix -> pathway x pathway matrix.
#'
#' Orchestrates stage1_nes() -> stage2_matrix() -> jaccard_mask() and global
#' BH adjustment over populated, unmasked cells. This is the single canonical
#' code path: production runs, calibration permutations, and sensitivity
#' analyses must all go through this function.
#'
#' @param z perturbations x genes matrix (rownames = perturbed gene symbols).
#' @param pathways named list of gene-symbol vectors.
#' @param stage1,stage2 named lists of overrides forwarded to stage1_nes()
#'   and stage2_matrix().
#' @param jaccard_tau mask threshold; NA_real_ disables masking.
#' @param BPPARAM a \code{BiocParallel} parameter object forwarded to both
#'   stages (a \code{BPPARAM} entry inside \code{stage1}/\code{stage2}
#'   takes precedence for that stage).
#' @param verbose emit per-stage timing messages.
#' @return list(nes, pval_stage1, M_raw, M, pval, padj, ..., config).
#'   `config` records the stage1/stage2 overrides and jaccard_tau of this
#'   call, so that calibrate_pvalues() can rerun the identical configuration
#'   on permuted data without the user re-passing it.
#' @section Reproducibility:
#' \code{BiocParallel} draws worker seeds from entropy unless the backend
#' carries an explicit seed, so \code{set.seed()} alone does not make a
#' run reproducible. For exact reproducibility pass a seeded backend,
#' e.g. \code{BPPARAM = BiocParallel::MulticoreParam(RNGseed = 1)} (or
#' \code{SerialParam(RNGseed = 1)}).
#' @examples
#' data(demo_z)
#' data(demo_pathways)
#' res <- path2path(demo_z, demo_pathways,
#'     stage1 = list(nperm = 100, min_size = 5),
#'     stage2 = list(min_pu = 5),
#'     BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
#' res$M[1:3, 1:3]
#' sum(res$padj < 0.05, na.rm = TRUE)
#' @export
path2path <- function(z, pathways,
                      stage1     = list(),
                      stage2     = list(),
                      jaccard_tau = 0.5,
                      BPPARAM    = BiocParallel::bpparam(),
                      verbose    = TRUE) {
    .check_matrix_input(z, "z")
    .check_pathways(pathways)
    if (!is.numeric(jaccard_tau) || length(jaccard_tau) != 1L ||
        (!is.na(jaccard_tau) &&
            (!is.finite(jaccard_tau) || jaccard_tau < 0 || jaccard_tau > 1))) {
        stop("`jaccard_tau` must be NA or one number between 0 and 1")
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
        stop("`verbose` must be TRUE or FALSE")
    if (is.null(stage1$BPPARAM)) stage1$BPPARAM <- BPPARAM
    if (is.null(stage2$BPPARAM)) stage2$BPPARAM <- BPPARAM
    t0 <- proc.time()
    s1 <- do.call(stage1_nes, c(list(z = z, pathways = pathways), stage1))
    if (verbose)
        message(sprintf("Stage 1 done in %.0f s (%.1f%% NA)",
            (proc.time() - t0)["elapsed"], 100 * mean(is.na(s1$nes))))

    t0 <- proc.time()
    s2 <- do.call(
        stage2_matrix, c(list(nes = s1$nes, pathways = pathways), stage2)
    )
    if (verbose)
        message(sprintf(
            "Stage 2 done in %.0f s", (proc.time() - t0)["elapsed"]
        ))

    J <- if (!is.na(jaccard_tau)) pathway_jaccard(pathways) else NULL
    mask <- function(M) {
        if (is.null(J) || is.null(M)) M else
            jaccard_mask(M, pathways, jaccard_tau, J)
    }

    out <- list(nes = s1$nes, pval_stage1 = s1$pval)
    if (!is.null(s2$M)) {                       # score_type = "std"
        out$M_raw <- s2$M
        out$M     <- mask(s2$M)
        out$pval  <- mask(s2$pval)
        out$padj  <- bh_populated(out$pval)
    } else {                                    # score_type = "both"
        out$M_pos <- mask(s2$M_pos)
        out$pval_pos <- mask(s2$pval_pos)
        out$M_neg <- mask(s2$M_neg)
        out$pval_neg <- mask(s2$pval_neg)
        out$padj_pos <- bh_populated(out$pval_pos)
        out$padj_neg <- bh_populated(out$pval_neg)
    }
    out$config <- list(stage1 = stage1, stage2 = stage2,
        jaccard_tau = jaccard_tau)
    out
}

#' Benjamini-Hochberg over the populated (non-NA) cells of a matrix.
#'
#' Exported utility: encodes the package's multiple-testing family rule (one
#' BH family over all populated entries), for use on any p-value matrix
#' produced by the package (e.g. negative-control stage2_matrix() output).
#' @param P a numeric matrix of p-values (NA = unpopulated/masked entries).
#' @return matrix of the same shape with BH-adjusted values on populated
#'   entries.
#' @examples
#' P <- matrix(c(0.001, 0.04, NA, 0.5), 2, 2)
#' bh_populated(P)
#' @export
bh_populated <- function(P) {
    if (is.null(P)) return(NULL)
    if (!is.matrix(P) || !is.numeric(P)) stop("`P` must be a numeric matrix")
    valid <- P[!is.na(P)]
    if (any(!is.finite(valid)) || any(valid < 0 | valid > 1))
        stop("non-missing entries of `P` must be finite values between 0 and 1")
    idx <- which(!is.na(P))
    P[idx] <- stats::p.adjust(P[idx], method = "BH")
    P
}
