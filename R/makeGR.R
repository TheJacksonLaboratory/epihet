#' @title Make List of GenomicRanges Object
#'
#' @description
#' Creates a GenomicRanges object for each methclone output file
#'
#' @param files A vector of input files containing methclone output files,
#' the suffix of files should be methClone_out.gz
#' @param ids A vector of sample ids for the files
#' @param cores The number of cores to be used for parallel execution
#' (default: 5)
#' @param sve A boolean to save the GenomicRanges object (default: FALSE)
#' @return A list, each element is a data frame of GenomicRanges objects 
#' containing pdr, epipolymorphism, and Shannon entropy values for each 
#' input file. Saves as an epi.gr.rda extension
#' @examples
#' path = system.file('extdata', package = 'epihet')
#' files = dir(path = path, pattern = 'methClone_out.gz',
#'             recursive = TRUE, full.names = TRUE)
#' ids = basename(dirname(files))
#' GR.List = epihet::makeGR(files = files, ids = ids,
#' cores = 1, sve = FALSE)
#' @export
makeGR = function(files, ids, cores = 5, sve = FALSE) {


  # validate Methclone files

  for (f in files) {
    parts = unlist(strsplit(f, "[.]"))
    parts = tail(parts,n=2)
    extension = paste(parts[1],parts[2],sep=".")
    if (extension != "methClone_out.gz") {
      message = paste(f,
          "file type is not supported, only supports files generated from methclone",
          sep = " ")
      stop(message)
    }
  }
  doParallel::registerDoParallel(cores = cores)
  n = NULL
  epi.gr = foreach(n = 1:length(files)) %dopar% {
    f = files[n]
    x = data.table::fread(paste("gzip -dc", f), sep = "\t")  #[,-27]
    x=x[,-27]
    x$pdr = rowSums(x[, c(12:25), with = FALSE])/100
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
    return(epi.gr)
  }
}
