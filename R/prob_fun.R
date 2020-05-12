# dependencies: -

#' @title Event Distribution
#' @description !!!!!!!!!!!!! shall not be exported !!!!!!!!!!!!!!!!!!!
#'
#' @param r vector of recruitment timepoints.
#' @param lambda,sigma scale and shape parameters for the survival process, the latter defaulting to 1,
#' resulting in an exponential distribution.
#' @param gamma,kappa scale and shape parameters fot the censor process, the latter defaulting to 1,
#' resulting in an exponential distribution.
#'
#' @return
#' \itemize{
#' \item Gives the distribution of the events.
#' \item Invalid arguments will result in return value \code{NaN}, with warning.
#' }
#'
#' @export
#' @examples
#' prob_fun(r=1:10,.014,sigma=1,.009,kappa=1)
#' prob_fun(r=c(1:10),.014,sigma=1,.009,kappa=1)
prob_fun <- function(r,lambda,sigma=1,gamma,kappa=1){
  dweibull(r,sigma,lambda^(-1/sigma)) * (1-pweibull(r,kappa,gamma^(-1/kappa)))
}
