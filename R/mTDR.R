#' mTDR fit a mixture copula model to comprehensively evaluate the reproducibile
#' of chromatin profiling sequencing data, including ChIP-seq (narrow peak and
#' broad peak), ATAC-seq, DNAse-seq.
#'
#' @details
#' \itemize{
#'   \item{Package:    }{mTDR}
#'   \item{Type:    }{Package}
#'   \item{Version:    }{0.99.1}
#'   \item{Date:    }{2018-2-1}
#'   \item{License:    }{GPL-2}
#'   \item{LazyLoad:    }{Yes}
#' }
#'
#' The main functions are \code{\link{est.TDR}}. The function \code{\link{est.TDR}}
#' implements expectation maximization (EM) algorithm to estimate the parameters of
#' mixture model of survival Clayton and Guassian.
#' @author
#' Tao Yang
#' Maintainer: Tao Yang <xadmyangt@gmail.com>
#' @references To be added.
#' @examples
#' p <- 0.25
#' rho <- 0.2
#' beta <- 2
#' max.it <- 100
#' eps <- 0.001
#' data(chipseq)
#' TDR.out = est.TDR(chipseq$R1, chipseq$R2, p, rho, beta, max.it, eps)
#' TDR.out$para


"_PACKAGE"
