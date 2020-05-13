# dependencies: prob_fun()

#' @title Calculates the expected number of events
#' @description !!!!!!!!!!!!! shall not be exported !!!!!!!!!!!!!!!!!!!
#' Calculates the expected number of events of a population given their respective
#' recruiting timepoint, assumed survival time of the first group, the hazard ratio
#' and ce censor rate. If the distributions are ....
#'
#' @param N a two column matrix containing the size of both recruitment groups per
#' row, representing the number of recruited patients per group and timepoint (\code{tn}).
#' @param tn a numeric vector defining the recruitment timepoint. Default is \code{0:(NROW(N)-1)}.
#' @param theta a number greater 0 defining the assumed hazard ratio.
#' @param L a number greater 0 defining the timepoint for administrative censoring.
#' If no administrative censoring is planned \code{L=Inf}. Default is L=Inf.
#' @param lambda a number greater 0 defining the hazard rate for the survival process of
#' group1. For the parametrization see \code{details}.
#' @param gamma a number greater 0 defining the overall hazard rate for the censor process.
#' @param sigma,kappa a number greater 0 defining the shape parameter for the
#' survival- and the censor process. Only needed if distS or distC is set to weibull.
#' Default is 1 resulting in an exponential distribution.
#' @param distS,distC a character string defining the distribution of the survival-
#' and the censor process. Default is 'exponential'.
#' @param col_Sums a boolen value. If set to FALSE, the function  will return the expected
#' number of events per recruitment batch at time T instead of the sums per group.
#' @details  With constant hazard rate \eqn{\lambda} and shape parameter \eqn{\sigma} the survival function for the weibull event process is given by
#' \deqn{S ( t ) = exp{ \lambda t ^ \sigma },}
#' with density function
#' \deqn{f ( t ) = \lambda \sigma t ^ { \sigma - 1 } exp{ \lambda t ^ \sigma }.}
#' @return a vector containing the expected number of events per group at time L.
#' @export
#'
#' @examples
#' # with a recruitment time of 5 month and recruiting 30 Patients per month
#' # assuming a hazard ratio of .7 and a rate of two year survival of .7
#' # as well as a two year censor rate of .2 under exponential distribution
#' # and administrative censoring after 10 month.
#' N <- matrix(rep(50,10),ncol=2)
#' calc_expEvents(N=N ,theta = .7, L=10, lambda=-log(.7)/24, gamma=-log(.8)/24)
calc_expEvents <- function(N, theta, L=Inf, tn,
                           lambda, sigma=1, distS="exponential",
                           gamma, kappa=1, distC="exponential",
                           col_Sums = TRUE){
  if (NCOL(N) !=2)
    stop("N must be a tow-column matrix")
  R = min(L,NROW(N))
  if(missing(tn))
    tn = 0:(R-1)
  if(L<NROW(N)){
    N <- N[1:L,]
    warning("The administrative censoring takes place before recruitment is finished.")
  }


  lambda.vec = matrix(c(lambda,lambda*theta),ncol=2)
  ### ??? might be fastend by slecting and exchanging pdweibull by pdesxp
  if ((distC=="exponential")&(distS=="exponential")){
    if (col_Sums)
      return(colSums((N*(1-exp(-(lambda.vec + gamma) %x% (L-tn)))) *
                     matrix(rep(lambda.vec/(lambda.vec+gamma),each=R),ncol=2)))
    return((N*(1-exp(-(lambda.vec + gamma) %x% (L-tn)))) *
                     matrix(rep(lambda.vec/(lambda.vec+gamma),each=R),ncol=2))
  } else {
    tmpN <- N*0
    for (i in 1:R){
      tmpN[i,] <- N[i,]*c(integrate(prob.fun,0,(L-tn[i]),lambda=lambda.vec[1],gamma=gamma,sigma=sigma,kappa=kappa)$value,
                          integrate(prob.fun,0,(L-tn[i]),lambda=lambda.vec[2],gamma=gamma,sigma=sigma,kappa=kappa)$value)
    }
    if (col_Sums)
      return(colSums(tmpN))
    return(tmpN)
  }
}
