#' @title Summarize Data
#' @description
#' Summarizes pdr, epi, and shannon values over the annotation regions
#' @param  gr1 A GenomicRaes object to be compared
#' @param  gr2 A GenomicRaes object to be compared
#' @param  value1 The value of gr1 to be compared
#' @param  value2 The value of gr2 be compared
#' @param  cutoff1 The first cutoff value for the number of reads
#' @param  cutoff2 The second cutoff value for the number of reads
#' @return A data frame containing a summary of the GenomicRaes object
#' @examples
#' p1.GR<-GRanges(seqnames = Rle(c("chr22"), c(5)),
#' ranges = IRanges(c(327,821,838,755,761), end = c(364,849,858,773,781)),
#' strand = Rle(strand(c("-", "+", "+", "+", "-"))),
#' values.loci = c("327:350:361:364","821:837:844:849",
#' "838:845:850:858","755:761:771:773","761:771:773:781"),
#' values.read1 = c(92,72,68,176,176),values.meth1=c(84,93,94,96,95),
#' values.shannon=c(0.4,0.5,0.5,0.2,0.5),values.pdr=c(0.6,0.25,0.23,0.15,0.17),
#' values.epipoly=c(0.48,0.42,0.38,0.27,0.3))
#'
#' p2.GR<-GRanges(seqnames = Rle(c("chr22"), c(5)),
#' ranges = IRanges(c(327,821,838,755,761), end = c(364,849,858,773,781)),
#' strand = Rle(strand(c("-", "+", "+", "+", "-"))),
#' values.loci = c("327:350:361:364","821:837:844:849",
#' "838:845:850:858","755:761:771:773","761:771:773:781"),
#' values.read1 = c(107,102,102,76,76),values.meth1=c(88,66,69,71,94),
#' values.shannon=c(0.12,0.25,0.54,0.23,0.25),
#' values.pdr=c(0.38,1,0.97,1,0.13),
#' values.epipoly=c(0.57,0.42,0.28,0.18,0.23))
#'
#' GR.List=list(p1=p1.GR,p2=p2.GR)
#' summary = epihet::summarize(gr1 = GR.List[[1]], gr2 = GR.List[[2]],
#' value1 = 'pdr', value2 = 'epipoly',
#' cutoff1 = 10, cutoff2 = 60)
#' @importFrom stats cor
#' @importFrom GenomicRaes findOverlaps values
#' @importFrom S4Vectors queryHits
#' @export
summarize = function(gr1, gr2, value1, value2, cutoff1 = 10,
    cutoff2 = 60) {
    values = c("pdr", "shannon", "epipoly")
    if (!(value1 %in% values)) {
        stop("Invalid value '", value1, "': Possible values are 'pdr',
             'epipoly', or 'shannon'")
    }
    if (!(value2 %in% values)) {
        stop("Invalid value '", value2, "': Possible values are 'pdr',
             'epipoly', or 'shannon'")
    }
    o = findOverlaps(gr1, gr2)
    x.anno = values(gr1[unique(queryHits(o))])
    sub1 = x.anno[x.anno$values.read1 >= cutoff1, ]
    sub2 = x.anno[x.anno$values.read1 >= cutoff2, ]
    mean.value1.c1 = mean(sub1$values.pdr)
    mean.value2.c1 = mean(sub1$values.epi)
    mean.value1.c2 = mean(sub2$values.pdr)
    mean.value2.c2 = mean(sub2$values.epi)
    name = c(paste0("mean.", value1, ".", cutoff1),
        paste0("mean.", value2, ".", cutoff1), paste0("mean.",
            value1, ".", cutoff2), paste0("mean.",
            value2, ".", cutoff2), paste0("corr.",
            cutoff1), paste0("corr.", cutoff2), paste0("loci.",
            cutoff1), paste0("loci.", cutoff2))
    df = data.frame(mean.value1.cutoff1 = mean.value1.c1,
        mean.value2.cutoff1 = mean.value2.c1,
        mean.value1.cutoff2 = mean.value1.c2,
        mean.value2.cutoff2 = mean.value2.c2,
        corr.cutoff1 = cor(sub1$values.pdr,
            sub1$values.epi), corr.cutoff2 = cor(sub2$values.pdr,
            sub2$values.epi), loci.cutoff1 = nrow(sub1),
        loci.cutoff2 = nrow(sub2))
    colnames(df) = name
    df
}
