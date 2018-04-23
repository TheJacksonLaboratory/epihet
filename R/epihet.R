#' @importFrom foreach foreach %dopar%
#' @importFrom GenomicRanges seqnames  GRanges findOverlaps values
#' @importFrom doMC registerDoMC
#' @importFrom stats p.adjust t.test
#' @importFrom ggplot2 element_text geom_boxplot scale_fill_manual
#' scale_x_discrete scale_y_continuous theme autoplot geom_text geom_bar ggsave
#' @importFrom graphics plot points
#' @importFrom stats complete.cases quantile sd prcomp
#' @importFrom grDevices colorRampPalette dev.off pdf rgb
#' @importFrom stats sd
#' @importFrom Rtsne Rtsne
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors Rle queryHits
#' @importFrom data.table :=
#' @importFrom Rcpp sourceCpp
#' @importFrom  WGCNA pickSoftThreshold blockwiseModules moduleEigengenes
#' labeledHeatmap exportNetworkToCytoscape
#' @importFrom  ReactomePA enrichPathway
#' @importFrom  pheatmap pheatmap
#' @importFrom igraph layout_on_sphere
#' @useDynLib epihet
NULL
