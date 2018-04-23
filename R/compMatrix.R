#' @title Make Comparison Matrix
#'
#' @description
#' A matrix is created for pca/hclust/tsne which contains
#' read number, average methylation levels, pdr, epipoly,
#' and Shannon entropy values across multiple samples at
#' the same loci using read number in a GenomicRanges object
#'
#' @param epi.gr An input file containing the read number, locus,
#' pdr, epipoly, and Shannon entropy values stored in a list of
#' GenomicRanges objects
#' @param outprefix The prefix name of the outputted matrix file.
#' 'sve' must be set to TRUE (default: NULL)
#' @param readNumber The lowest number of reads required for each loci
#' (default: 60)
#' @param p Percentage (as decimal) of matching samples required
#' to determine a match at a given locus, e.g. a value of 0.75
#' requires 75\% of the samples to have an epiallele at a common
#' loci in order to add the loci to the matrix (default: 1)
#' @param cores The number of cores to be used for parallel execution
#' (default: 5)
#' @param sve A boolean to save the comparison matrix (default: FALSE)
#' @return A large matrix containing values (epi, pdr, etc.) at the same loci
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
#' N1.GR<-GRanges(seqnames = Rle(c("chr22"), c(5)),
#' ranges = IRanges(c(327,821,838,755,761), end = c(364,849,858,773,781)),
#' strand = Rle(strand(c("-", "+", "+", "+", "-"))),
#' values.loci = c("327:350:361:364","821:837:844:849",
#' "838:845:850:858","755:761:771:773","761:771:773:781"),
#' values.read1 = c(112,112,112,68,76),values.meth1=c(82,60,91,71,90),
#' values.shannon=c(0.15,0.26,0.34,0.24,0.15),
#' values.pdr=c(0.32,0.57,0.37,0.37,0.13),
#' values.epipoly=c(0.57,0.42,0.28,0.38,0.23))
#'
#' N2.GR<-GRanges(seqnames = Rle(c("chr22"), c(5)),
#' ranges = IRanges(c(327,821,838,755,761), end = c(364,849,858,773,781)),
#' strand = Rle(strand(c("-", "+", "+", "+", "-"))),
#' values.loci = c("327:350:361:364","821:837:844:849",
#' "838:845:850:858","755:761:771:773","761:771:773:781"),
#' values.read1 = c(385,78,70,96,96),values.meth1=c(96,81,87,87,93),
#' values.pdr=c(0.15,0.52,0.48,0.25,0.25),
#' values.epipoly=c(0.26,0.58,0.58,0.37,0.37),
#' values.shannon=c(0.12,0.25,0.54,0.23,0.25))
#'
#' GR.List=list(p1=p1.GR,p2=p2.GR,N1=N1.GR,N2=N2.GR)
#' comp.Matrix = epihet::compMatrix(epi.gr = GR.List, outprefix = NULL,
#' readNumber = 60, p = 1, cores = 1, sve = FALSE)
#' @export
compMatrix = function(epi.gr, outprefix = NULL, readNumber = 60,
    p = 1, cores = 5, sve = FALSE) {
    doMC::registerDoMC(cores = cores)
    i = NULL
    print("Getting all loci")
    sub.ids = foreach(x = epi.gr, .combine = c) %dopar% {
        x = x[values(x)$values.read1 >= readNumber,]
        paste(seqnames(x), values(x)$values.loci, sep = "-")
    }
    print(paste0("Shared ids by ", p * 100, "% samples"))
    shared.ids = names(which(table(sub.ids) >= p * length(epi.gr)))
    print("Shared epimatrix")
    dat.na = data.frame(value = rep(NA, length(shared.ids)),
        row.names = shared.ids)
    epi.shared.matrix = foreach(x = epi.gr, .combine = cbind) %dopar% {
        x = x[values(x)$values.read1 >= readNumber,]
        ids = paste(seqnames(x), values(x)$values.loci, sep = "-")
        x1 = x[ids %in% shared.ids, ]
        x2 = values(x1)[, 2:6]
        rownames(x2) = paste(seqnames(x1), values(x1)[,1], sep = "-")
        use.ids = intersect(ids, shared.ids)
        res = foreach(i = colnames(values(x1))[2:6],.combine = rbind) %dopar% {
            dat = dat.na
            # dat$type=gsub('values.|1','',i)
            dat[use.ids, "value"] = x2[use.ids, i]
            dat
        }
        res
    }
    colnames(epi.shared.matrix) = names(epi.gr)
    epi.shared.matrix$type = rep(gsub("values.|1", "",
        colnames(values(epi.gr[[1]]))[2:6]), each = length(shared.ids))
    epi.shared.matrix$location = rownames(epi.shared.matrix)
    epi.shared.matrix[epi.shared.matrix$type != "read", "location"] =
        gsub(".{1}$", "", epi.shared.matrix[epi.shared.matrix$type !=
        "read", "location"])
    rownames(epi.shared.matrix) = 1:nrow(epi.shared.matrix)
    if (sve) {
        save(epi.shared.matrix, file = paste0(outprefix,
            "_epi_shared.matrix.rda"))
    } else {
        epi.shared.matrix
    }
}
