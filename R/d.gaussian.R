#' function to caluculate the density of Guassian copula
#'
#' @param u The cumulative density of a random variable. It can be calculated using `empdist` function.
#' @param v The cumulative density of another random variable.
#' @param rho The correlation paramter of Guassian copula.
#' @return The density of Clayton Copula.
#' @importFrom stats qnorm
#' @export


d.gaussian <- function(u, v, rho){

    dgau <- (1/sqrt(1 - rho^2))*exp((2*rho*qnorm(u)*qnorm(v)
                                   - rho^2*(qnorm(u)^2 + qnorm(v)^2))/(2*(1-rho^2)))
    return(dgau)

}
