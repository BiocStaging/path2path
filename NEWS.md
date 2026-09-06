# path2path 0.99.2

* Replaced repeated per-ranking `fgsea` tied-statistics warnings with one
  summary per stage. The summary reports the number of affected rankings and
  their median and maximum tie fractions; other warnings remain visible.
* Applied the same tie reporting to the pooled implementations and added
  regression tests across Stage 1 and Stage 2 modes.

# path2path 0.99.1

* Added conversion to and from `SummarizedExperiment`.
* Strengthened input validation, reproducibility, and documentation in
  preparation for Bioconductor review.
* When no adequate GPD tail fit can be obtained, `calibrate_pvalues()` returns
  pooled empirical p-values instead of extrapolating from an inadequate fit.

# path2path 0.99.0

* Initial Bioconductor submission candidate, implementing two-stage pathway
  mapping, two significance routes, overlap masking, plotting, and file
  input/output.
