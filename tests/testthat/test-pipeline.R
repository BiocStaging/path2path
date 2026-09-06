data(demo_z)
data(demo_pathways)
SER <- BiocParallel::SerialParam(RNGseed = 100)
CFG <- list(stage1 = list(nperm = 200, min_size = 5),
    stage2 = list(min_pu = 5))

run_demo <- function(seed = 1) {
    set.seed(seed)
    path2path(demo_z, demo_pathways, stage1 = CFG$stage1, stage2 = CFG$stage2,
        BPPARAM = BiocParallel::SerialParam(RNGseed = seed),
        verbose = FALSE)
}

test_that("pipeline recovers the planted directional relations", {
    res <- run_demo()
    expect_identical(
        dim(res$M), c(length(demo_pathways), length(demo_pathways))
    )
    # A activates B; C suppresses D (planted); direction matters
    expect_gt(res$M["PATHWAY_A", "PATHWAY_B"], 1.5)
    expect_lt(res$M["PATHWAY_C", "PATHWAY_D"], -1.5)
    expect_lt(res$padj["PATHWAY_A", "PATHWAY_B"], 0.05)
    expect_lt(res$padj["PATHWAY_C", "PATHWAY_D"], 0.05)
    # the reverse entries are untestable rows (B, D members unperturbed): NA
    expect_true(is.na(res$M["PATHWAY_B", "PATHWAY_A"]))
    # E has no planted targets: its off-diagonal entries are clearly weaker
    # than the planted ones (an absolute no-significance assertion would be
    # flaky by construction: BH thresholds relax around the two planted
    # ~1e-29 entries in this tiny 18-entry family)
    e_row <- res$padj["PATHWAY_E", setdiff(colnames(res$padj), "PATHWAY_E")]
    expect_gt(min(e_row, na.rm = TRUE),
        1e6 * res$padj["PATHWAY_A", "PATHWAY_B"])
})

test_that("results are reproducible under a fixed seed (serial)", {
    r1 <- run_demo(7)
    r2 <- run_demo(7)
    expect_identical(r1$M, r2$M)
    expect_identical(r1$pval, r2$pval)
})

test_that("target exclusion removes the self-hit from Stage 1", {
    set.seed(1)
    s_excl <- stage1_nes(demo_z, demo_pathways, nperm = 100, min_size = 5,
        exclude_target = TRUE, BPPARAM = SER)
    set.seed(1)
    s_incl <- stage1_nes(demo_z, demo_pathways, nperm = 100, min_size = 5,
        exclude_target = FALSE, BPPARAM = SER)
    # with the target retained, each A-perturbation's own pathway rides the
    # -(-8) self-hit; excluded, that boost disappears
    a <- demo_pathways$PATHWAY_A
    expect_gt(mean(s_incl$nes[a, "PATHWAY_A"], na.rm = TRUE),
        mean(s_excl$nes[a, "PATHWAY_A"], na.rm = TRUE))
})

test_that("jaccard_mask masks overlapping pairs and never the diagonal", {
    pw <- list(X = letters[1:10], Y = letters[1:9], Z = letters[15:24])
    M <- matrix(1, 3, 3, dimnames = list(names(pw), names(pw)))
    Mm <- jaccard_mask(M, pw, tau = 0.5)
    expect_true(is.na(Mm["X", "Y"]) && is.na(Mm["Y", "X"]))
    expect_false(any(is.na(diag(Mm))))
    expect_false(is.na(Mm["X", "Z"]))
})

test_that("bh_populated adjusts only populated cells", {
    P <- matrix(c(0.01, NA, 0.02, 0.5), 2, 2)
    Q <- bh_populated(P)
    expect_identical(is.na(P), is.na(Q))
    expect_identical(sort(Q[!is.na(Q)]),
        sort(stats::p.adjust(P[!is.na(P)], "BH")))
})

test_that("input validation rejects malformed matrices", {
    bad <- demo_z
    rownames(bad) <- NULL
    expect_error(path2path(bad, demo_pathways), "rownames")
    expect_error(path2path(demo_z, list(1, 2)), "named list")
    expect_error(path2path(demo_z, list(A = 1:3)), "character vector")
    expect_error(path2path(demo_z, demo_pathways, jaccard_tau = 2),
        "between 0 and 1")
})

test_that("SummarizedExperiment conversion preserves matrix orientation", {
    skip_if_not_installed("SummarizedExperiment")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(z = t(demo_z))
    )
    expect_identical(as_path2path_matrix(se, "z"), demo_z)
    res <- run_demo(13)
    map_se <- as_path2path_se(res)
    expect_identical(SummarizedExperiment::assay(map_se, "M"), res$M)
    expect_identical(SummarizedExperiment::assayNames(map_se),
        c("M", "M_raw", "pval", "padj"))
})

test_that("write_result creates its output directory", {
    res <- run_demo(17)
    d <- tempfile("path2path-output-")
    write_result(res, d, "demo", "hallmark")
    expect_true(dir.exists(d))
    expect_true(file.exists(file.path(
        d, "M_matrix_demo_hallmark_fgsea.csv"
    )))
})

capture_warning_messages <- function(expr) {
    messages <- character()
    withCallingHandlers(
        invisible(force(expr)),
        warning = function(w) {
            messages <<- c(messages, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    messages
}

test_that("tied Stage 1 rankings produce one summary warning", {
    z_tied <- round(demo_z * 2) / 2
    serial <- BiocParallel::SerialParam(RNGseed = 1)

    messages <- capture_warning_messages(stage1_nes(
        z_tied, demo_pathways, nperm = 50, min_size = 5,
        BPPARAM = serial))
    expect_equal(sum(grepl("tied statistics", messages)), 1L)

    messages <- capture_warning_messages(stage1_nes(
        z_tied, demo_pathways, method = "pooled", min_size = 5,
        b_prime = 2, BPPARAM = serial))
    expect_equal(sum(grepl("tied statistics", messages)), 1L)

    expect_no_warning(stage1_nes(
        demo_z, demo_pathways, nperm = 50, min_size = 5,
        BPPARAM = serial))
})

test_that("tied Stage 2 rankings produce one summary in every mode", {
    n <- nrow(demo_z) * length(demo_pathways)
    nes_tied <- matrix(
        rep(c(-2, -1, 0, 1, 2), length.out = n),
        nrow = nrow(demo_z), ncol = length(demo_pathways),
        dimnames = list(rownames(demo_z), names(demo_pathways)))
    serial <- BiocParallel::SerialParam(RNGseed = 2)
    run_stage2 <- function(...) {
        stage2_matrix(
            nes_tied, demo_pathways, min_pu = 5, min_valid = 5,
            BPPARAM = serial, ...)
    }

    messages <- capture_warning_messages(run_stage2())
    expect_equal(sum(grepl("tied statistics", messages)), 1L)

    messages <- capture_warning_messages(run_stage2(score_type = "both"))
    expect_equal(sum(grepl("tied statistics", messages)), 1L)

    messages <- capture_warning_messages(run_stage2(
        method = "pooled", b_prime = 2))
    expect_equal(sum(grepl("tied statistics", messages)), 1L)
})

test_that("the ties handler passes unrelated warnings through", {
    expect_warning(
        path2path:::.run_fgsea_ties_quiet({
            warning("sentinel warning")
            1
        }),
        "sentinel warning")
})
