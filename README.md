# path2path

Directional pathway-by-pathway interaction maps from genome-scale
Perturb-seq, by two stacked runs of preranked GSEA: once per perturbation
over genes, then once per readout pathway over perturbations. Entry (j, k)
of the resulting map measures whether perturbing pathway j's members
produces a coordinated transcriptional response in pathway k.

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

```r
library(path2path)
z <- read_input_matrix("zscore_myscreen.csv.gz")  # perturbations x genes
pathways <- load_pathways("hallmark")
res <- path2path(z, pathways)
plot_map(res)

# whole-pipeline calibration for individually examined entries:
cal <- calibrate_pvalues(res, z, pathways, B = 20, seed_base = 1)
```

See the vignette (`vignette("path2path")`) for the full tour, and the
accompanying manuscript for the method's calibration and validation
program.

## License

MIT (c) 2026 Jun Li
