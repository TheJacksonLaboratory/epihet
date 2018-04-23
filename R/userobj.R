#' @title GenomicRanges object generation
#' @description generate GenomicRanges object for DEH loci
#' @param  data a data frame containing the chromosome number,
#' loci and strand information of DEH loci generated from diffHet() function.
#' @return A GenomicRanges object
#' @examples
#' data=data.frame(chromosome=c('chr1','chr1','chr1'),
#' loci=c('6480531:6480554:6480561:6480565','6647655:6647696:6647701:6647705',
#' '7130155:7130172:7130179:7130188'),
#' strand=c('+','-','+'),stringsAsFactors = FALSE)
#' userobj(data)
#' @export
userobj = function(data) {
    data$start = as.numeric(splitn(data$loci, ":", 
        1))
    data$end = as.numeric(splitn(data$loci, ":", 4))
    userset = GenomicRanges::GRanges(seqnames = Rle(data$chromosome), 
        ranges = IRanges(start = data$start, end = data$end), 
        strand = data$strand, loci = paste(data$chromosome, 
            data$loci, sep = "-"))
    userset
    return(userset)
}
