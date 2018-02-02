#' EM algorithm to estimate the parameters of mixture model of survival Clayton
#' copula and Gaussian copula.
#'
#' @param x1 A vector of signal values. For example, the read counts in genomic bins.
#' @param x2 Another vector of same length signal values that to be test reproducibility
#' against x1.
#' @param p The starting value of proportion of Guassion copula for the EM algorithm.
#' The default is 0.5.
#' @param rho The starting value of Guassian correlation parameter for the EM algorithm.
#' The default is 0.2. It can be any value betwen 0 and 1.
#' @param beta The starting value of Clayton shape parameter for the EM algorithm.
#' The default is 2. It can be any value greater than 0.
#' @param max.it The maximum number of iterations for EM algorithm. The default is 100.
#' @param eps The threhold of convergence. Iteration of EM stops if all the changes of
#' parameters are smaller than this threshold.
#' @return A list of results including parameter traces,log likelihood trace and final
#' estimated parameters.
#' \itemize{
#'   \item{para.trace    }{The trace of parameters along the estimation iterations.}
#'    \item{likelihood.trace    }{The trace of log likelihood. The trace of log likelihood
#'    can use to monitor the estimation. It should be monotonically increasing afte
#'    first few iterations.}
#'    \item{para    }{The estimated parameters, including pi, rho, lambda.}
#' }
#' @details This function estimates the parameters for the mixture copula of survival
#' Clayton and Gaussian.
#' @references To be added.
#' @importFrom stats optimize qnorm ecdf
#' @export
#' @examples
#' p <- 0.25
#' rho <- 0.2
#' beta <- 2
#' max.it <- 100
#' eps <- 0.001
#' data(chipseq)
#' TDR.out = est.TDR(chipseq$R1, chipseq$R2, p, rho, beta, max.it, eps)
#' TDR.out$para

est.TDR <- function(x1, x2, p=0.25, rho=0.2, beta=2, max.it=100, eps=0.001) {

    U <- empdist(x1, x2)
    u <- U[,1]
    v <- U[,2]

    if (p <=0 && p>=1 && !is.numeric(p))
        stop("invalid argument: intialize p with numeric number between 0 and 1 (not include)\n")
        else{
            p.trace <- array()
            rho.trace <- array()
            beta.trace <- array()
            loglik.gau <- array()
            loglik.cla <- array()
            loglik.T <- array()

            i <- 0
            repeat {
                i <- i+1

                d.gau <- d.gaussian(u, v, rho)
                d.clay <- d.clayton(1-u, 1-v, beta)

                ki.p <- p*d.clay/((1-p)*d.gau+p*d.clay)
                p.p <- sum(ki.p)/length(ki.p)

                fgau <- function(rho){
                    sum((1-ki.p)*(log(1-p.p)+log(d.gaussian(u, v, rho))))
                }

                optima.gau <- optimize(f = fgau, c(0 + sqrt(.Machine$double.eps), 1), maximum = T)
                rho.p <- optima.gau$maximum
                loglik.gau[i] <- optima.gau$objective


                fclay <- function(beta){
                    sum(ki.p*(log(p.p)+log(d.clayton(1-u, 1-v, beta))))
                }
                optima.cla <- optimize(f = fclay, c(0 + sqrt(.Machine$double.eps), 20), maximum = T)
                beta.p <- optima.cla$maximum
                loglik.cla[i] <- optima.cla$objective

                dif.p <- abs(p.p-p)
                dif.rho <- abs(rho.p-rho)
                dif.be <- abs(beta.p-beta)

                if(dif.p > eps){p = p.p} else {p = p}
                if(dif.rho > eps){rho = rho.p} else {rho = rho}
                if(dif.be > eps){beta = beta.p} else {beta = beta}

                p.trace[i] <- p
                rho.trace[i] <- rho
                beta.trace[i] <- beta

                loglik.T[i] <- loglik.gau[i] + loglik.cla[i]

                if (dif.p < eps & dif.rho < eps & dif.be < eps | i == max.it){
                    if (i == max.it) warning("You estimation may not converge, try a greater 'max.it'")
                    break
                }
            }


            p <- rev(p.trace)[1]
            rho <- rev(rho.trace)[1]
            beta <- rev(beta.trace)[1]
            tailD <-  2^(-1/beta)
            a.rho <- (1-p)*rho
            a.tailD <- p*tailD
            para <- c(p = 1-p, rho = rho, lambda = a.tailD)

            para.trace <- cbind(p.trace, rho.trace, beta.trace)

            likelihood.trace <- cbind(loglik.gau, loglik.cla, loglik.T)
            est <- list(para.trace = para.trace, likelihood.trace = likelihood.trace, para = para)
        }

    return(est)
}
