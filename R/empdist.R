#' calculate the emperical cumulative density
#'
#' @param x1 A vector of signal values. For example, the read counts in genomic bins.
#' @param x2 Another vector of same length signal values that to be test reproducibility against x1.
#' @return A data frame with two columns and the same length with x1 and x2.
#' @importFrom stats ecdf
#' @export

empdist <- function(x1, x2){

    x1r <- rank(x1, ties.method = "random")
    x2r <- rank(x2, ties.method = "random")

    x1.cdf.func <- ecdf(x1r)
    x2.cdf.func <- ecdf(x2r)

    afactor <- length(x1r)/(length(x1r) + 1)

    x1.cdf <- x1.cdf.func(x1r) * afactor
    x2.cdf <- x2.cdf.func(x2r) * afactor

    u <- cbind(x1.cdf, x2.cdf)
    return(u)
}
