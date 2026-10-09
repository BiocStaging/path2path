# path2path

**What did my screen's perturbations do to each other's pathways?**

A genome-scale Perturb-seq experiment knocks down thousands of genes and
reads out the whole transcriptome for each one. Standard enrichment tools
summarize each perturbation separately; the systems-level question — which
*pathway*, when perturbed, moves which *other pathway*, and in which
direction — is usually left to manual inspection.

`path2path` answers that question in one call. It turns the screen's
perturbation-by-gene matrix of differential-expression statistics into a
**directional pathway-by-pathway interaction map**: entry (j, k) is a
signed enrichment score with associated significance estimates,
measuring whether perturbing members of pathway *j* coherently shifts
the transcriptional program of pathway *k*, in the direction that, under
the loss-of-function sign convention, suggests activation (> 0) or
suppression (< 0). Because rows are perturbed pathways and columns are
readouts, the map is asymmetric — the candidate relationship "j acts on
k" is scored separately from "k acts on j", a distinction no
correlation- or co-expression-based pathway similarity can make.

Under the hood this is two stacked runs of preranked GSEA (fgsea): once
per perturbation over genes, then once per readout pathway over
perturbations, with three safeguards (on-target exclusion,
effective-size filters, pathway-overlap masking) and two significance
routes — fast nominal p-values with BH screening for exploration, and a
whole-pipeline permutation calibration with generalized-Pareto tail
extrapolation for individually reported entries. In the accompanying
manuscript, maps built this way replicate across independent screens,
platforms, and laboratories (for example, 98% strong-entry sign agreement
between K562 and iPSC screens analyzed with different statistics), and
condition contrasts of the maps highlight pathway rewiring that
single-map readings miss.

Use it when you have (or can export) a perturbation × gene matrix of
signed DE statistics from any Perturb-seq-style screen — CRISPRi or
knockout, or CRISPRa with `stage1 = list(negate = FALSE)` — and want a
quantitative, significance-annotated directional pathway map rather
than per-pathway lists.

## Installation

```r
# once accepted on Bioconductor:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("path2path")

# development version:
remotes::install_github("jli-stat/path2path")
```

## Quick start

Runs as-is on the bundled synthetic demo screen (two planted
relationships; both are recovered with the right sign and direction):

```r
library(path2path)
data(demo_z)          # 120 perturbations x 400 genes
data(demo_pathways)   # 6 pathways

res <- path2path(demo_z, demo_pathways,
    stage1 = list(nperm = 200, min_size = 5),
    stage2 = list(min_pu = 5),
    BPPARAM = BiocParallel::SerialParam(RNGseed = 1))

res$M["PATHWAY_A", "PATHWAY_B"]      # planted activation:  +3.5
res$M["PATHWAY_C", "PATHWAY_D"]      # planted suppression: -3.5
plot_map(res, padj_cutoff = 0.05)

# whole-pipeline permutation calibration for individually examined entries:
cal <- calibrate_pvalues(res, demo_z, demo_pathways, B = 50, seed_base = 42)
cal$pval_gpd["PATHWAY_A", "PATHWAY_B"]
```

## Using your own data

Replace the demo objects with your own matrix and an MSigDB collection.
`zscore_myscreen.csv.gz` below is a placeholder for your screen's
perturbation x gene CSV of signed DE statistics (first column =
perturbed gene symbols):

```r
z <- read_input_matrix("zscore_myscreen.csv.gz")  # <- your own file
pathways <- load_pathways("hallmark")             # or "reactome", ...
res <- path2path(z, pathways)
```

See the vignette (`vignette("path2path")`) for the full tour, and the
accompanying manuscript for the method's calibration and validation
program.

## License

MIT (c) 2026 Jun Li
