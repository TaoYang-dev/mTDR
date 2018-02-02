#' function to caluculate the density of Clayton copula
#'
#' @param u The cumulative density of a random variable. It can be calculated using `empdist` function.
#' @param v The cumulative density of another random variable.
#' @param beta The shape parameter for Clayton copula.
#' @return The density of Clayton Copula.
#' @export


d.clayton <- function(u, v, beta){

    m <- u^(-beta)+v^(-beta)-1

    dclay <- (1+beta)*((u*v)^(-1-beta))*(u^(-beta)+v^(-beta)-1)^(-1/beta-2)

    return(dclay)
}
