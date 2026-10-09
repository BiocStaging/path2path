# path2path 0.99.3

Changes responding to the first round of Bioconductor review:

* `path2path()` now defaults to `verbose = FALSE`; `build_null_pool()` and
  `calibrate_pvalues()` gained a `verbose` argument (default `FALSE`) gating
  the per-permutation progress messages, which were previously unconditional.
  Messages that signal a change of behavior (GPD fit-adequacy fallback,
  collection-name aliasing, non-default configuration reuse, names
  appended beyond a user-supplied plot order) remain.
* All imported functions are now declared with `importFrom()` in NAMESPACE.
* Function documentation states the class of every parameter.
* The vignette gained a section on how the package integrates with the
  Bioconductor ecosystem and how it differs from per-contrast enrichment
  packages; the quick-tour chunks explain each parameter and interpret the
  printed output; and all code chunks are now evaluated, including the
  Route B calibration (run at `B = 50` on the demo data) and the
  file-based input example, which reads the demo matrix shipped in
  `inst/extdata/zscore_demo.csv.gz`.
* README rewritten around the user-facing question the package answers,
  with a copy-paste quick start that runs on the bundled demo data.
* Added the WholeGenome biocView.

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
