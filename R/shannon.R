#' @title Shannon Entropy
#'
#' @description
#' Calculates the Shannon entropy value
#'
#' @param  p A vector of epiallele probabilities

#' @return The Shannon entropy value
#' @examples
#' a=c(rep(0,13),60.86960,0.00000,39.1304)
#' shannon(a)
#' @export
shannon = function(p) {
    p1 = p[p > 0]
    sum(-p1/100 * log(p1/100))
}
