#' Pairwise Jaccard similarity of pathway gene sets (dense crossproduct).
#' @param pathways named list of gene-symbol vectors.
#' @return K x K similarity matrix (diagonal set to 0; it is never masked).
#' @examples
#' data(demo_pathways)
#' J <- pathway_jaccard(demo_pathways)
#' max(J)
#' @export
pathway_jaccard <- function(pathways) {
    .check_pathways(pathways)
    all_genes <- sort(unique(unlist(pathways)))
    A <- matrix(FALSE, length(all_genes), length(pathways),
        dimnames = list(all_genes, names(pathways)))
    for (pw in names(pathways))
        A[intersect(all_genes, pathways[[pw]]), pw] <- TRUE
    inter <- crossprod(A)
    sizes <- colSums(A)
    J <- inter / (outer(sizes, sizes, "+") - inter)
    diag(J) <- 0   # diagonal never masked downstream
    J
}

#' Mask off-diagonal entries whose pathway pair overlaps above `tau`.
#'
#' High-overlap pairs (Reactome parent/child etc.) score mechanically —
#' j's members ride k's signal because they ARE k's members — which is
#' gene-set structure, not causal regulation. Diagonal is never masked.
#' @param M K x K matrix to mask.
#' @param pathways named list of gene-symbol vectors.
#' @param tau overlap threshold; entries with Jaccard > tau become NA.
#' @param jaccard optional precomputed \code{pathway_jaccard()} matrix.
#' @return M with high-overlap off-diagonal entries set to NA.
#' @examples
#' data(demo_pathways)
#' M <- matrix(1, length(demo_pathways), length(demo_pathways),
#'     dimnames = list(names(demo_pathways), names(demo_pathways)))
#' sum(is.na(jaccard_mask(M, demo_pathways, tau = 0.1)))
#' @export
jaccard_mask <- function(M, pathways, tau = 0.5, jaccard = NULL) {
    if (!is.matrix(M) || !is.numeric(M)) stop("`M` must be a numeric matrix")
    .check_pathways(pathways)
    .check_probability(tau, "tau")
    if (is.null(jaccard)) jaccard <- pathway_jaccard(pathways)
    if (!is.null(dimnames(M)) && !is.null(dimnames(jaccard))) {
        if (!all(rownames(M) %in% rownames(jaccard)) ||
            !all(colnames(M) %in% colnames(jaccard)))
            stop(
                "jaccard_mask: M dimnames are not covered by the Jaccard matrix"
            )
        jaccard <- jaccard[rownames(M), colnames(M)]
    } else {
        stopifnot(identical(dim(M), dim(jaccard)))
    }
    M[jaccard > tau] <- NA
    M
}
