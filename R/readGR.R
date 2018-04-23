#' @title Make GenomicRanges Object
#'
#' @description
#' Creates a GenomicRanges file for a singular bisulfite sequencing
#' data file
#'
#' @param files A vector of files containing bisulfite sequencing data
#' @param ids A vector of sample id's for the files
#' @param n The index of the file vector to be read
#' @return A GenomicRanges object containing pdr, epipolymorphism, and
#' Shannon entropy values for the nth file
#' @examples
#' path = system.file('extdata', package = 'epihet')
#' files = dir(path = path, pattern = '\\.bam$',
#'             recursive = TRUE, full.names = TRUE)
#' ids = basename(dirname(files))
#' GR.Object = epihet::readGR(files = files, ids = ids, n = 3)
#' @export
readGR = function(files, ids, n) {
    f = files[n]
    sampleid = ids[n]

    parts = unlist(strsplit(f, "[.]"))

    extension = tail(parts, n = 1)
    if (extension != "bam") {
        message = paste(f,
            "file type is not supported, only supports BAM files", sep = " ")
        stop(message)
    }
    index_file = paste(f, "bai", sep = ".")
    if (!file.exists(index_file)) {
        message = paste("BAM index file for", f, "is not provided", sep = " ")
        stop(message)
    }
    f_out = paste(f, "methClone_out.gz", sep = ".")
    MethClone_one_sample(f, f_out, sampleid, 72, 60)
    x = data.table::fread(paste("zcat", f_out), sep = "\t")  #[,-27]
    x=x[,-27]
    x$pdr = 1 - rowSums(x[, c(11, 26), with = FALSE])/100
    x$epipoly = 1 - rowSums((x[, c(11:26), with = FALSE]/100)^2)
    x$shannon = apply(x[, c(11:26), with = FALSE], 1, shannon)
    x = as.data.frame(x)
    x.gr = GRanges(Rle(x$chr), IRanges(start = x$start,
        end = x$end), strand = x$strand, values = x[, c(7, 8, 9, 27:29)])
    x.gr
}
