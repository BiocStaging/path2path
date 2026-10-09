#' Read a z-score CSV (rows = perturbed gene symbols in first column).
#' @param path character(1): CSV or CSV.GZ file; first column = perturbed
#'   gene symbols,
#'   remaining columns = measured genes.
#' @return numeric matrix with dimnames.
#' @examples
#' data(demo_z)
#' f <- tempfile(fileext = ".csv")
#' write.csv(data.frame(perturbation = rownames(demo_z), demo_z,
#'     check.names = FALSE), f, row.names = FALSE)
#' z <- read_input_matrix(f)
#' identical(dim(z), dim(demo_z))
#' @export
read_input_matrix <- function(path) {
    dt <- data.table::fread(path)
    z <- as.matrix(dt[, -1, with = FALSE])
    rownames(z) <- dt[[1]]
    .check_matrix_input(z, "z")
    z
}

#' Convert a SummarizedExperiment assay to a path2path input matrix.
#'
#' A SummarizedExperiment stores features in rows and samples or perturbations
#' in columns. path2path uses the transpose: perturbations in rows and measured
#' genes in columns. This helper performs that conversion while retaining the
#' assay dimnames.
#'
#' @param x a `SummarizedExperiment` with genes in rows and perturbations
#'   in columns.
#' @param assay character(1) or integer(1): assay name or index passed to
#'   `SummarizedExperiment::assay()`.
#' @return a numeric perturbation-by-gene matrix.
#' @examples
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     data(demo_z)
#'     se <- SummarizedExperiment::SummarizedExperiment(
#'         assays = list(z = t(demo_z))
#'     )
#'     z <- as_path2path_matrix(se, "z")
#'     stopifnot(identical(z, demo_z))
#' }
#' @export
as_path2path_matrix <- function(x, assay = 1L) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
        stop("Install `SummarizedExperiment` to use this function")
    if (!methods::is(x, "SummarizedExperiment"))
        stop("`x` must be a SummarizedExperiment")
    z <- t(as.matrix(SummarizedExperiment::assay(x, assay)))
    .check_matrix_input(z, "converted assay")
    z
}

#' Convert a path2path map to a SummarizedExperiment.
#'
#' The returned object stores the map and its p-value matrices as assays,
#' allowing the pathway-by-pathway result to be used with standard
#' Bioconductor infrastructure. The Stage 1 matrix remains available in the
#' original result object.
#'
#' @param result a list, a standard `path2path()` result containing `M`.
#' @return a `SummarizedExperiment` with pathways as rows and readout
#'   pathways as columns.
#' @examples
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     data(demo_z)
#'     data(demo_pathways)
#'     res <- path2path(
#'         demo_z, demo_pathways,
#'         stage1 = list(nperm = 100, min_size = 5),
#'         stage2 = list(min_pu = 5),
#'         BPPARAM = BiocParallel::SerialParam(), verbose = FALSE
#'     )
#'     se <- as_path2path_se(res)
#'     SummarizedExperiment::assayNames(se)
#' }
#' @export
as_path2path_se <- function(result) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
        stop("Install `SummarizedExperiment` to use this function")
    if (is.null(result$M) || !is.matrix(result$M))
        stop("`result` must be a standard path2path result containing `$M`")
    assays <- Filter(Negate(is.null), list(
        M = result$M,
        M_raw = result$M_raw,
        pval = result$pval,
        padj = result$padj
    ))
    SummarizedExperiment::SummarizedExperiment(
        assays = assays,
        metadata = list(path2path_config = result$config)
    )
}

#' Write pipeline outputs with dataset/collection namespacing.
#'
#' Files: nes_matrix_<ns>.csv.gz, pval_matrix_stage1_<ns>.csv.gz,
#'   M_matrix_<ns>_raw.csv, M_matrix_<ns>.csv, pval/padj_matrix_<ns>.csv
#'   (or *_pos/_neg variants for score_type = "both").
#'
#' If `calibrated` (a calibrate_pvalues() result) is supplied, Route B
#' outputs are written alongside: pval_gpd_matrix_<ns>.csv,
#'   pval_pool_matrix_<ns>.csv, gpd_params_<ns>.csv, and the null pool
#'   as null_pool_<ns>.rds (for reuse via calibrate_pvalues(pool =)).
#' @param result a list as returned by \code{path2path()}.
#' @param dir character(1): output directory.
#' @param tag,collection,method character(1) each: components of the file
#'   namespace.
#' @param calibrated `NULL` or a list as returned by
#'   \code{calibrate_pvalues()}.
#' @return (invisible) the namespace string used in filenames.
#' @examples
#' data(demo_z)
#' data(demo_pathways)
#' res <- path2path(demo_z, demo_pathways,
#'     stage1 = list(nperm = 100, min_size = 5),
#'     stage2 = list(min_pu = 5),
#'     BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
#' d <- tempdir()
#' write_result(res, d, "demo", "demo6", method = "fgsea")
#' list.files(d, pattern = "demo")
#' @export
write_result <- function(result, dir, tag, collection, method = "fgsea",
                         calibrated = NULL) {
    if (!is.character(dir) || length(dir) != 1L || is.na(dir))
        stop("`dir` must be one non-missing path")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    ns <- sprintf("%s_%s_%s", tag, collection, method)
    wm <- function(M, name, rn = "pathway_j") {
        if (is.null(M)) return(invisible(NULL))
        data.table::fwrite(
            data.table::as.data.table(M, keep.rownames = rn),
            file.path(dir, name))
    }
    data.table::fwrite(
        data.table::as.data.table(result$nes, keep.rownames = "perturbation"),
        file.path(dir, sprintf("nes_matrix_%s.csv.gz", ns)))
    if (!is.null(result$pval_stage1))
        data.table::fwrite(
            data.table::as.data.table(
                result$pval_stage1, keep.rownames = "perturbation"
            ),
            file.path(dir, sprintf("pval_matrix_stage1_%s.csv.gz", ns)))
    if (!is.null(result$M)) {
        wm(result$M_raw, sprintf("M_matrix_%s_raw.csv", ns))
        wm(result$M,     sprintf("M_matrix_%s.csv", ns))
        wm(result$pval,  sprintf("pval_matrix_%s.csv", ns))
        wm(result$padj,  sprintf("padj_matrix_%s.csv", ns))
    } else {
        wm(result$M_pos,    sprintf("M_matrix_%s_pos.csv", ns))
        wm(result$M_neg,    sprintf("M_matrix_%s_neg.csv", ns))
        wm(result$pval_pos, sprintf("pval_matrix_%s_pos.csv", ns))
        wm(result$pval_neg, sprintf("pval_matrix_%s_neg.csv", ns))
        wm(result$padj_pos, sprintf("padj_matrix_%s_pos.csv", ns))
        wm(result$padj_neg, sprintf("padj_matrix_%s_neg.csv", ns))
    }
    if (!is.null(calibrated)) {
        wm(calibrated$pval_gpd,  sprintf("pval_gpd_matrix_%s.csv", ns))
        wm(calibrated$padj_gpd,  sprintf("padj_gpd_matrix_%s.csv", ns))
        wm(calibrated$pval_pool, sprintf("pval_pool_matrix_%s.csv", ns))
        if (!is.null(calibrated$gpd))
            data.table::fwrite(data.table::as.data.table(calibrated$gpd),
                file.path(dir, sprintf("gpd_params_%s.csv", ns)))
        saveRDS(calibrated$pool_abs,
            file.path(dir, sprintf("null_pool_%s.rds", ns)))
    }
    invisible(ns)
}

#' Dataset tag from a zscore_<tag>.csv[.gz] filename.
#' @param path character(1): a file path such as \code{zscore_myexp.csv.gz}.
#' @return the tag string (\code{"myexp"}).
#' @examples
#' dataset_tag("zscore_k562_gwps.csv.gz")
#' @export
dataset_tag <- function(path)
    sub("\\.csv(\\.gz)?$", "", sub("^zscore_", "", basename(path)))
