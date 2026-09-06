#' Internal input checks shared by the pipeline entry points.
#' @noRd
.check_matrix_input <- function(z, what = "z") {
    if (!is.matrix(z) || !is.numeric(z))
        stop("`", what, "` must be a numeric matrix (rows x named columns)")
    if (!nrow(z) || !ncol(z))
        stop("`", what, "` must have at least one row and one column")
    if (is.null(rownames(z)) || is.null(colnames(z)))
        stop("`", what, "` must have both rownames and colnames (gene symbols)")
    if (anyNA(rownames(z)) || anyNA(colnames(z)) ||
        any(!nzchar(rownames(z))) || any(!nzchar(colnames(z))))
        stop("`", what, "` must have non-empty, non-missing dimnames")
    if (any(is.infinite(z)))
        stop("`", what, "` cannot contain infinite values")
    if (anyDuplicated(rownames(z))) {
        duplicates <- unique(rownames(z)[duplicated(rownames(z))])
        stop(
            "duplicated rownames in `", what,
            "` (fgsea requires unique names); e.g. ",
            paste(utils::head(duplicates, 3), collapse = ", ")
        )
    }
    if (anyDuplicated(colnames(z))) {
        duplicates <- unique(colnames(z)[duplicated(colnames(z))])
        stop(
            "duplicated colnames in `", what,
            "` (fgsea requires unique names); e.g. ",
            paste(utils::head(duplicates, 3), collapse = ", ")
        )
    }
    invisible(TRUE)
}

#' Internal check for a pathway collection.
#' @noRd
.check_pathways <- function(pathways) {
    if (!is.list(pathways) || length(pathways) == 0L ||
        is.null(names(pathways)) || any(!nzchar(names(pathways))))
        stop("`pathways` must be a non-empty named list of gene-symbol vectors")
    if (anyDuplicated(names(pathways)))
        stop("duplicated pathway names in `pathways`")
    valid <- vapply(pathways, function(x) {
        is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x))
    }, logical(1))
    if (!all(valid))
        stop(
            "every element of `pathways` must be a character vector ",
            "without missing or empty gene symbols")
    invisible(TRUE)
}

#' Check that an argument is a scalar integer value.
#' @noRd
.check_integer <- function(x, what) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        !is.finite(x) || x != floor(x)) {
        stop("`", what, "` must be one finite integer")
    }
    invisible(TRUE)
}

#' Check that an argument is a positive scalar integer.
#' @noRd
.check_positive_integer <- function(x, what) {
    .check_integer(x, what)
    if (x < 1L) stop("`", what, "` must be positive")
    invisible(TRUE)
}

#' Check that an argument is a scalar probability.
#' @noRd
.check_probability <- function(x, what, open = FALSE) {
    valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
    valid <- valid && if (open) x > 0 && x < 1 else x >= 0 && x <= 1
    if (!valid) {
        interval <- if (open) "strictly between 0 and 1" else "between 0 and 1"
        stop("`", what, "` must be one finite number ", interval)
    }
    invisible(TRUE)
}

#' Run fgseaMultilevel while silencing its non-actionable log2err warning.
#'
#' fgsea returns the conservative p-value when its log2 error estimate cannot
#' be calculated. In a pathway-by-pathway map this warning can otherwise be
#' repeated once per readout pathway. All other fgsea warnings remain visible.
#' @noRd
.run_fgsea_multilevel <- function(expr) {
    withCallingHandlers(force(expr), warning = function(w) {
        if (grepl("P-values were likely overestimated", conditionMessage(w),
            fixed = TRUE)) {
            invokeRestart("muffleWarning")
        }
    })
}

#' Muffle only fgsea's per-ranking tied-statistics warning.
#'
#' The repeated warnings are replaced by one stage-level summary through
#' `.warn_ties_summary()`. All other warnings remain visible.
#' @noRd
.run_fgsea_ties_quiet <- function(expr) {
    withCallingHandlers(force(expr), warning = function(w) {
        if (grepl("ties in the preranked stats", conditionMessage(w),
            fixed = TRUE)) {
            invokeRestart("muffleWarning")
        }
    })
}

#' Emit one summary warning for tied statistics in ranked lists.
#'
#' Ties can make enrichment scores depend on the relative input order of
#' equally ranked genes. Reporting the condition once per stage preserves this
#' information without producing one warning for every ranked list.
#' @noRd
.warn_ties_summary <- function(tie_fraction, where) {
    affected <- which(tie_fraction > 0)
    if (!length(affected)) return(invisible(FALSE))
    median_percent <- formatC(
        100 * stats::median(tie_fraction[affected]),
        format = "f", digits = 1)
    maximum_percent <- formatC(
        100 * max(tie_fraction[affected]),
        format = "f", digits = 1)
    warning(
        where, ": ", length(affected), " of ", length(tie_fraction),
        " ranked lists contain tied statistics (median tie fraction ",
        median_percent, "%, maximum ", maximum_percent,
        "%). The relative order of tied genes depends on input order and ",
        "can affect enrichment scores.",
        call. = FALSE)
    invisible(TRUE)
}

#' Fraction of entries in a numeric vector that repeat an earlier value.
#' @noRd
.tie_fraction <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(0)
    sum(duplicated(x)) / length(x)
}
