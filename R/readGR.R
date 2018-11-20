#' @title Make GenomicRanges Object
#'
#' @description
#' Creates a GenomicRanges file for a singular methclone ouput file
#'
#' @param files A vector of files containing methcolone output
#' @param ids A vector of sample ids for the files
#' @param n The index of the file vector to be read
#' @return A GenomicRanges object containing pdr, epipolymorphism, and
#' Shannon entropy values for the nth file
#' @examples
#' files = c(system.file("extdata","D-2238.chr22.region.methClone_out.gz",package = "epihet"),
#' system.file("extdata","D-2668.chr22.region.methClone_out.gz",package = "epihet"),
#' system.file("extdata","N-1.chr22.region.methClone_out.gz",package = "epihet"),
#' system.file("extdata","N-2.chr22.region.methClone_out.gz",package = "epihet"))
#' ids = epihet::splitn(basename(files),"[.]",1)
#' GR.Object = epihet::readGR(files = files, ids = ids, n = 3)
#' @export
readGR = function(files, ids, n) {
    f = files[n]
    sampleid = ids[n]

    parts = unlist(strsplit(f, "[.]"))
    parts = tail(parts,n=2)
    extension = paste(parts[1],parts[2],sep=".")
    if (extension != "methClone_out.gz") {
        message = paste(f,
            "file type is not supported, only supports files generated from methclone", sep = " ")
        stop(message)
    }

    x = data.table::fread(paste("gzip -dc", f), sep = "\t")  #[,-27]
    x=x[,-27]
    x$pdr = rowSums(x[, c(12:25), with = FALSE])/100
    x$epipoly = 1 - rowSums((x[, c(11:26), with = FALSE]/100)^2)
    x$shannon = apply(x[, c(11:26), with = FALSE], 1, shannon)
    x = as.data.frame(x)
    x.gr = GRanges(Rle(x$chr), IRanges(start = x$start,
        end = x$end), strand = x$strand, values = x[, c(7, 8, 9, 27:29)])
    x.gr
}
