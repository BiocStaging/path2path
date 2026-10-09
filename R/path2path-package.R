#' path2path: directional pathway-by-pathway maps from Perturb-seq
#'
#' The package converts a perturbation-by-gene matrix of signed
#' differential-expression statistics into a directional pathway-by-pathway
#' map using two successive preranked gene set enrichment analyses. It also
#' provides pathway-overlap safeguards, two routes for significance
#' assessment, plotting functions, and conversion to and from
#' `SummarizedExperiment` objects.
#'
#' See `vignette("path2path")` for an introduction and worked example.
#'
#' @references
#' Subramanian A, et al. (2005). Gene set enrichment analysis: a
#' knowledge-based approach for interpreting genome-wide expression profiles.
#' Proceedings of the National Academy of Sciences 102:15545-15550.
#' doi:10.1073/pnas.0506580102.
#'
#' Korotkevich G, Sukhov V, Sergushichev A (2019). Fast gene set enrichment
#' analysis. bioRxiv. doi:10.1101/060012.
#' @author Jun Li
#' @keywords internal
#'
#' @importFrom BiocParallel MulticoreParam SerialParam bplapply bpparam
#' @importFrom data.table as.data.table data.table fread fwrite
#' @importFrom fgsea calcGseaStat fgseaMultilevel fgseaSimple plotEnrichment
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline axis box image layout legend lines mtext par
#'   plot points rect
#' @importFrom methods is
#' @importFrom msigdbr msigdbr
#' @importFrom stats dist hclust median optim p.adjust quantile sd setNames
#' @importFrom utils head modifyList write.table
#' @importFrom withr with_seed
"_PACKAGE"
