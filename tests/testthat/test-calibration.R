data(demo_z)
data(demo_pathways)
SER <- BiocParallel::SerialParam(RNGseed = 5)

test_that("calibrate_pvalues returns coherent Route B values", {
    set.seed(2)
    res <- path2path(demo_z, demo_pathways,
        stage1 = list(nperm = 100, min_size = 5),
        stage2 = list(min_pu = 5),
        BPPARAM = SER, verbose = FALSE)
    cal <- calibrate_pvalues(res, demo_z, demo_pathways, B = 2, seed_base = 11)
    expect_identical(dim(cal$pval_gpd), dim(res$M))
    expect_identical(is.na(res$M), is.na(cal$pval_pool))
    ok <- !is.na(res$M) & !is.na(cal$pval_pool)
    # pooled p decreases with |M| (same pool for every entry)
    o <- order(abs(res$M[ok]))
    expect_true(all(diff(cal$pval_pool[ok][o]) <= 1e-12))
    # permutations are reproducible via seed_base
    cal2 <- calibrate_pvalues(res, demo_z, demo_pathways, B = 2, seed_base = 11)
    expect_identical(cal$pool_abs, cal2$pool_abs)
})

test_that("permute_input preserves values and labels", {
    set.seed(3)
    zp <- permute_input(demo_z)
    expect_identical(dimnames(zp), dimnames(demo_z))
    expect_identical(sort(as.vector(zp)), sort(as.vector(demo_z)))
    expect_false(identical(zp, demo_z))
})

test_that("seed_base does not alter the caller's random-number stream", {
    set.seed(91)
    before <- .Random.seed
    build_null_pool(
        demo_z, demo_pathways, B = 1, seed_base = 7,
        stage1 = list(nperm = 100, min_size = 5),
        stage2 = list(min_pu = 5),
        BPPARAM = SER
    )
    expect_identical(.Random.seed, before)
})

test_that("calibration rejects malformed pools and thresholds", {
    res <- path2path(
        demo_z, demo_pathways,
        stage1 = list(nperm = 100, min_size = 5),
        stage2 = list(min_pu = 5),
        BPPARAM = SER, verbose = FALSE
    )
    expect_error(calibrate_pvalues(res, pool = numeric()), "non-empty")
    expect_error(calibrate_pvalues(res, pool = c(0.1, NA)), "finite")
    expect_error(calibrate_pvalues(res, pool = 1:100, q_thresh = 1),
        "strictly between")
})
