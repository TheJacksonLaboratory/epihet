#' @title  The subtring extraction of a character vector
#' @description Extract the subtrings of a character vector according to
#' the matches to substring split within them.
#' @param  strings A GenomicRanges object to be compared
#' @param  field A GenomicRanges object to be compared
#' @param  n The value of gr1 to be compared
#' @return A data frame containing a summary of the GenomicRanges object
#' @examples
#' x='chr1:10000-10005'
#' splitn(x,':',1)
#' @export
splitn = function(strings, field, n) {
    sapply(strsplit(strings, field), "[", n)
}
