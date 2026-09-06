#' Internal: prettify MSigDB pathway names for axis labels.
#' @noRd
.pretty_names <- function(x, max_chars = 24L) {
    x <- sub("^(HALLMARK|REACTOME|KEGG|WP|PID)_", "", x)
    x <- tolower(gsub("_", " ", x))
    ifelse(
        nchar(x) > max_chars,
        paste0(substr(x, 1, max_chars - 1), "\u2026"), x
    )
}

#' Internal: two-directional clustering order (rows of [A | t(A)]).
#' @noRd
.cluster_order <- function(A) {
    A0 <- A
    A0[is.na(A0)] <- 0
    if (nrow(A) == ncol(A) && identical(rownames(A), colnames(A))) {
        ord <- stats::hclust(
            stats::dist(cbind(A0, t(A0))), method = "average"
        )$order
        list(rows = rownames(A)[ord], cols = colnames(A)[ord])
    } else {
        ro <- stats::hclust(stats::dist(A0), method = "average")$order
        co <- stats::hclust(stats::dist(t(A0)), method = "average")$order
        list(rows = rownames(A)[ro], cols = colnames(A)[co])
    }
}

#' Internal: apply a user-supplied order to surviving names.
#' @noRd
.resolve_order <- function(order, current, axis) {
    ord <- if (is.list(order)) order[[axis]] else order
    if (is.null(ord)) return(current)
    keep  <- intersect(ord, current)
    extra <- setdiff(current, keep)
    if (length(extra))
        message("plot_map: ", length(extra), " ", axis,
            " not covered by `order` appended at the end")
    c(keep, extra)
}

#' Heatmap of a path2path map with significance dots.
#'
#' Processing pipeline (in this fixed priority): (1) `rows`/`cols` subset;
#' (2) all-NA rows/columns dropped if `drop_empty_*`; (3) ordering --
#' `order` if given (a character vector, or a list(rows=, cols=) as
#' returned invisibly by a previous plot_map call, so a second map can be
#' drawn under the first one's ordering), else two-directional
#' hierarchical clustering if `cluster = TRUE`, else matrix order.
#' Dots mark entries passing the significance cutoffs (see .sig_matrix
#' semantics: non-NULL cutoffs conjoin; `calibrated` switches to Route B
#' p-values). Masked/filtered entries are gray.
#'
#' @inheritParams select_active
#' @param rows,cols optional name subsets applied before anything else.
#' @param drop_empty_rows,drop_empty_cols drop all-NA rows/columns.
#' @param order a character vector or list(rows=, cols=) fixing the order.
#' @param cluster two-directional hierarchical clustering when no `order`.
#' @param labels force axis labels on/off (default: on up to 60 pathways).
#' @param zlim color range in NES units.
#' @param main plot title.
#' @param colorbar draw the NES colorbar.
#' @param dot_cex significance-dot size (default scales with map size).
#' @return (invisible) list(rows, cols): the plotted ordering, reusable as
#'   the `order` argument for another map.
#' @examples
#' data(demo_z)
#' data(demo_pathways)
#' res <- path2path(demo_z, demo_pathways,
#'     stage1 = list(nperm = 100, min_size = 5),
#'     stage2 = list(min_pu = 5),
#'     BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
#' plot_map(res, padj_cutoff = 0.2)
#' @export
plot_map <- function(result, calibrated = NULL,
                     padj_cutoff = 0.05, pval_cutoff = NULL,
                     rows = NULL, cols = NULL,
                     drop_empty_rows = TRUE, drop_empty_cols = TRUE,
                     order = NULL, cluster = TRUE,
                     labels = NULL, zlim = c(-3, 3),
                     main = NULL, colorbar = TRUE, dot_cex = NULL) {
    M <- result$M
    if (is.null(M))
        stop(
            "result has no $M ",
            "(score_type='both' results are not supported by plot_map)"
        )
    S <- .sig_matrix(result, calibrated, padj_cutoff, pval_cutoff)

    if (!is.null(rows)) M <- M[rows, , drop = FALSE]
    if (!is.null(cols)) M <- M[, cols, drop = FALSE]
    if (drop_empty_rows) M <- M[rowSums(!is.na(M)) > 0, , drop = FALSE]
    if (drop_empty_cols) M <- M[, colSums(!is.na(M)) > 0, drop = FALSE]
    if (!nrow(M) || !ncol(M)) stop("nothing left to plot after subsetting")

    if (!is.null(order)) {
        ro <- .resolve_order(order, rownames(M), "rows")
        co <- .resolve_order(order, colnames(M), "cols")
    } else if (cluster && nrow(M) > 2 && ncol(M) > 2) {
        ord <- .cluster_order(M)
        ro <- ord$rows
        co <- ord$cols
    } else {
        ro <- rownames(M)
        co <- colnames(M)
    }
    M <- M[ro, co, drop = FALSE]
    S <- S[ro, co, drop = FALSE]
    nr <- nrow(M)
    nc <- ncol(M)

    pal  <- grDevices::colorRampPalette(c("#3d6fb6", "#f0f0ee", "#c0392b"))(255)
    gray <- "#dddcd8"
    Mc <- pmin(pmax(M, zlim[1]), zlim[2])

    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    if (colorbar)
        graphics::layout(matrix(seq_len(2), 1, 2), widths = c(1, 0.16))

    show_lab <- if (is.null(labels)) max(nr, nc) <= 60 else isTRUE(labels)
    # label margins scale with the longest label so small maps fit small devices
    lab_lines <- function(v) min(10, max(3.5, 0.42 * max(nchar(v)) + 1))
    mar_left   <- if (show_lab) lab_lines(.pretty_names(ro)) else 2
    mar_bottom <- if (show_lab) lab_lines(.pretty_names(co)) else 3
    graphics::par(mar = c(mar_bottom, mar_left, if (is.null(main)) 1 else 3, 1))

    graphics::plot(NA, xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5),
        xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "",
        main = main)
    graphics::rect(0.5, 0.5, nc + 0.5, nr + 0.5, col = gray, border = NA)
    graphics::image(x = seq_len(nc), y = seq_len(nr),
        z = t(Mc[nr:1, , drop = FALSE]),
        col = pal, zlim = zlim, add = TRUE)
    w <- which(S, arr.ind = TRUE)
    if (nrow(w)) {
        cex <- if (is.null(dot_cex))
            max(0.15, min(0.7, 35 / max(nr, nc))) else dot_cex
        graphics::points(
            w[, 2], nr + 1 - w[, 1], pch = 16, cex = cex, col = "black"
        )
    }
    if (show_lab) {
        cexl <- max(0.25, min(0.8, 45 / max(nr, nc)))
        graphics::axis(1, at = seq_len(nc), labels = .pretty_names(co),
            las = 2, cex.axis = cexl, tick = FALSE, line = -0.5)
        graphics::axis(2, at = nr:1, labels = .pretty_names(ro),
            las = 2, cex.axis = cexl, tick = FALSE, line = -0.5)
    }
    graphics::mtext(
        "readout pathway k", side = 1, line = mar_bottom - 1, cex = 0.85
    )
    graphics::mtext(
        "perturbed pathway j", side = 2, line = mar_left - 1, cex = 0.85
    )
    graphics::box(lwd = 0.5)

    if (colorbar) {
        graphics::par(
            mar = c(mar_bottom, 0.6, if (is.null(main)) 1 else 3, 2.0)
        )
        zseq <- seq(zlim[1], zlim[2], length.out = 255)
        graphics::image(1, zseq, matrix(zseq, 1), col = pal, axes = FALSE,
            xlab = "", ylab = "")
        graphics::axis(4, cex.axis = 0.7, las = 1)
        graphics::mtext("NES", side = 4, line = 1.4, cex = 0.8)
        graphics::box(lwd = 0.5)
    }
    invisible(list(rows = ro, cols = co))
}

#' GSEA-style enrichment plot for one map entry: why is M[j, k] strong?
#'
#' Shows the Stage 2 ranking (perturbations ranked by their Stage 1 NES on
#' readout pathway k) with tick marks at the perturbations targeting
#' members of pathway j, via fgsea::plotEnrichment(). The title reports the
#' entry's NES and its p-value (Route B if `calibrated` is given).
#'
#' @param result a \code{path2path()} result.
#' @param pathways the collection used to build `result`.
#' @param j,k pathway names: the perturbed (row) and readout (column) side.
#' @param calibrated optional \code{calibrate_pvalues()} result; switches
#'   the title to Route B values.
#' @return a ggplot object.
#' @examples
#' data(demo_z)
#' data(demo_pathways)
#' res <- path2path(demo_z, demo_pathways,
#'     stage1 = list(nperm = 100, min_size = 5),
#'     stage2 = list(min_pu = 5),
#'     BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
#' if (requireNamespace("ggplot2", quietly = TRUE))
#'     plot_entry(res, demo_pathways, names(demo_pathways)[1],
#'         names(demo_pathways)[2])
#' @export
plot_entry <- function(result, pathways, j, k, calibrated = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("plot_entry needs ggplot2 (installed with fgsea)")
    if (!k %in% colnames(result$nes)) stop("unknown readout pathway `k`: ", k)
    if (!j %in% names(pathways)) stop("`j` not found in `pathways`: ", j)
    stats_k <- result$nes[, k]
    stats_k <- sort(stats_k[!is.na(stats_k)], decreasing = TRUE)
    members <- intersect(pathways[[j]], names(stats_k))
    if (length(members) < 2)
        warning("only ", length(members), " member(s) of ", j,
            " present in this column's ranking")
    g <- fgsea::plotEnrichment(members, stats_k)
    pv <- if (is.null(calibrated)) {
        sprintf("padj = %.2g", result$padj[j, k])
    } else {
        sprintf("calibrated p = %.2g (padj_gpd = %.2g)",
            calibrated$pval_gpd[j, k], calibrated$padj_gpd[j, k])
    }
    g + ggplot2::labs(
        title = sprintf(
            "%s \u2192 %s", .pretty_names(j, 40), .pretty_names(k, 40)
        ),
        subtitle = sprintf("NES = %+.2f, %s; %d perturbed members in ranking",
            result$M[j, k], pv, length(members)))
}

#' Route B diagnostics: null-pool tail versus the fitted GPD.
#'
#' Left: empirical survival curve of the null pool (log10 scale) with the
#' fitted GPD tail overlaid beyond the threshold u. Right: QQ plot of the
#' exceedances against fitted GPD quantiles. Both panels support the
#' judgment the user makes before quoting a pval_gpd: does the tail fit?
#' @param calibrated a \code{calibrate_pvalues()} result carrying
#'   \code{pool_abs} and \code{gpd}.
#' @return invisible(NULL); draws two base-graphics panels.
#' @examples
#' pool <- sort(-log(seq(0.9999, 0.0001, length.out = 5000)))
#' M <- matrix(c(0.5, 1, 1.5, 2), 2,
#'     dimnames = list(c("A", "B"), c("A", "B")))
#' cal <- calibrate_pvalues(list(M = M), pool = pool)
#' plot_calibration(cal)
#' @export
plot_calibration <- function(calibrated) {
    pool <- calibrated$pool_abs
    g <- calibrated$gpd
    if (is.null(pool) || is.null(g))
        stop("`calibrated` must carry pool_abs and gpd")
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))

    L <- length(pool)
    keep <- unique(round(seq(1, L, length.out = 4000)))   # thin for speed
    x <- pool[keep]
    surv <- (L - keep + 1) / L
    graphics::plot(x, surv, log = "y", type = "s", col = "gray40", lwd = 1.4,
        xlab = "|M'| (null pool)", ylab = "P(|M'| >= x)",
        main = "null tail vs GPD fit",
        xlim = c(0, max(pool) * 1.15))
    xs <- seq(g$u, min(max(pool) * 1.15,
        if (is.finite(g$support_end)) g$support_end else Inf),
    length.out = 200)
    ps <- .gpd_survival(xs, g$u, g$xi, g$beta, g$n_exc / g$n_tot)
    graphics::lines(xs[ps > 0], ps[ps > 0], col = "#c0392b", lwd = 2)
    graphics::abline(v = g$u, lty = 2, col = "gray55")
    graphics::mtext(sprintf("u = %.2f, xi = %.3f, beta = %.3f",
        g$u, g$xi, g$beta), cex = 0.75, line = 0.2)
    graphics::legend("topright", bty = "n", cex = 0.8,
        legend = c("empirical", "GPD tail"),
        col = c("gray40", "#c0392b"), lwd = c(1.4, 2))

    exc <- sort(pool[pool > g$u] - g$u)
    n <- length(exc)
    q <- (seq_len(n) - 0.5) / n
    theo <- if (abs(g$xi) > 1e-8) g$beta / g$xi * ((1 - q)^(-g$xi) - 1)
    else -g$beta * log(1 - q)
    graphics::plot(theo + g$u, exc + g$u, pch = 16, cex = 0.45, col = "gray30",
        xlab = "fitted GPD quantile", ylab = "observed exceedance",
        main = "exceedance QQ")
    graphics::abline(0, 1, col = "#c0392b", lwd = 1.5)
    invisible(NULL)
}
