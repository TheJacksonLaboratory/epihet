#' @title Make List of GenomicRanges Object
#'
#' @description
#' Creates a GenomicRanges object for each bisulfite sequencing data
#' file
#'
#' @param files A vector of input files containing bisulfite sequencing data
#' @param ids A vector of sample id's for the files
#' @param cores The number of cores to be used for parallel execution
#' (default: 5)
#' @param sve A boolean to save the GenomicRanges object (default: FALSE)
#' @return A data frame of GenomicRanges objects containing pdr,
#' epipolymorphism, and Shannon entropy values for all input files.
#' Saves as an epi.gr.rda extension
#' @examples
#' path = system.file('extdata', package = 'epihet')
#' files = dir(path = path, pattern = '\\.bam$',
#'             recursive = TRUE, full.names = TRUE)[c(1,2)]
#' ids = basename(dirname(files))[c(1,2)]
#' GR.List = epihet::makeGR(files = files, ids = ids,
#' cores = 1, sve = FALSE)
#' @export
makeGR = function(files, ids, cores = 5, sve = FALSE) {

    # validate Methclone files

    ## check if the file is in bam format
    for (f in files) {

        parts = unlist(strsplit(f, "[.]"))

        extension = tail(parts, n = 1)
        if (extension != "bam") {
            message = paste(f,
                "file type is not supported, only supports BAM files",
                sep = " ")
            stop(message)
        }


        index_file = paste(f, "bai", sep = ".")
        if (!file.exists(index_file)) {
            message = paste("BAM index file for", f,
                "is not provided", sep = " ")
            stop(message)
        }
    }



    doMC::registerDoMC(cores = cores)
    n = NULL
    epi.gr = foreach(n = 1:length(files)) %dopar% {
        f = files[n]
        sampleid = ids[n]
        f_out = paste(f, "methClone_out.gz", sep = ".")
        MethClone_one_sample(f, f_out, sampleid, 72, 60)
        x = data.table::fread(paste("zcat", f_out), sep = "\t")  #[,-27]
        x=x[,-27]
        x$pdr = 1 - rowSums(x[, c(11, 26), with = FALSE])/100
        x$epipoly = 1 - rowSums((x[, c(11:26), with = FALSE]/100)^2)
        x$shannon = apply(x[, c(11:26), with = FALSE], 1, shannon)
        x = as.data.frame(x)
        x.gr = GRanges(Rle(x$chr), IRanges(start = x$start, end = x$end),
            strand = x$strand, values = x[, c(7, 8, 9, 27:29)])
        x.gr
    }
    names(epi.gr) = ids
    if (sve) {
        save(epi.gr, file = "epi.gr.rda")
    } else {
        epi.gr
    }
}
