#' Load an MSigDB collection as a named list of gene-symbol vectors.
#'
#' Gene symbols are de-duplicated within each set. When `size_filter` is TRUE
#' the sets are restricted to `min_size <= |S| <= max_size` measured on the
#' *global* (genome) membership; the effective per-universe filters are
#' applied later inside stage1_nes()/stage2_matrix(), which is where they
#' matter statistically.
#'
#' @param collection one of "hallmark", "reactome", "kegg", "kegg_medicus",
#'   "kegg_legacy", "pid", "wikipathways". Note that "kegg" resolves to
#'   KEGG_MEDICUS (the current KEGG subcollection in MSigDB); use
#'   "kegg_legacy" for the classic 186-set collection.
#' @param min_size,max_size global size window applied when
#'   `size_filter` is active for the collection.
#' @param species passed to \code{msigdbr::msigdbr}.
#' @return named list of character vectors
#' @references
#' Dolgalev I (2026). msigdbr: MSigDB Gene Sets for Multiple Organisms in a
#' Tidy Data Format. R package.
#' @examples
#' hm <- load_pathways("hallmark")
#' length(hm)
#' @export
load_pathways <- function(collection,
                          min_size = 15L, max_size = 500L,
                          species = "Homo sapiens") {
    if (!is.character(collection) || length(collection) != 1L ||
        is.na(collection) || !nzchar(collection))
        stop("`collection` must be one non-empty string")
    .check_positive_integer(min_size, "min_size")
    .check_positive_integer(max_size, "max_size")
    if (min_size > max_size) stop("`min_size` cannot exceed `max_size`")
    if (!is.character(species) || length(species) != 1L ||
        is.na(species) || !nzchar(species))
        stop("`species` must be one non-empty string")
    if (tolower(collection) == "kegg")
        message(
            "load_pathways: 'kegg' resolves to KEGG_MEDICUS; use ",
            "'kegg_legacy' for the classic 186-set collection"
        )
    spec <- switch(tolower(collection),
        hallmark = list(
            collection = "H", subcollection = NA, size_filter = FALSE
        ),
        reactome = list(
            collection = "C2", subcollection = "CP:REACTOME", size_filter = TRUE
        ),
        kegg = list(
            collection = "C2", subcollection = "CP:KEGG_MEDICUS",
            size_filter = TRUE
        ),
        kegg_medicus = list(
            collection = "C2", subcollection = "CP:KEGG_MEDICUS",
            size_filter = TRUE
        ),
        kegg_legacy = list(
            collection = "C2", subcollection = "CP:KEGG_LEGACY",
            size_filter = TRUE
        ),
        pid = list(
            collection = "C2", subcollection = "CP:PID", size_filter = TRUE
        ),
        wikipathways = list(
            collection = "C2", subcollection = "CP:WIKIPATHWAYS",
            size_filter = TRUE
        ),
        stop("Unknown collection '", collection, "'")
    )
    args <- list(species = species, collection = spec$collection)
    if (!is.na(spec$subcollection)) args$subcollection <- spec$subcollection
    db <- do.call(msigdbr::msigdbr, args)
    pw <- lapply(split(db$gene_symbol, db$gs_name), unique)
    if (spec$size_filter)
        pw <- pw[lengths(pw) >= min_size & lengths(pw) <= max_size]
    pw
}
