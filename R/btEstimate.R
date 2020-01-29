# dependend on: rEstimate(), bParEstimate()

#' @title Blinded Study Time Estimation
#' @description !!!!!!!!!!!!! shall not be exported !!!!!!!!!!!!!!!!!!!
#'
#' @param dfSurv a data.frame containing the variables:
#' \itemize{
#'   \item status: (0,1) := (no event,event)
#'   \item time: time to status
#'   \item group: grouping variable for a
#' }
#' @param N a two column matrix containing the size of both recruitment groups per
#' row, representing the number of recruited patients per group and timepoint (\code{tn}).
#' @param tn a vecotr of length = \code{NROW(N)} defining the timepoint of recruiting.
#' @param M a two column matrix containing the size of both recruitment groups per
#' row, representing the number of future recruited patients per group and timepoint (\code{tm}).
#' If not defined it will be taken from the last row of \code{N}.
#' @param tm a vecotr of length = \code{NROW(M)} defining the timepoint of future recruiting.
#' If not defined it will be interpolated from the last two entries of \code{tn}.
#' @param theta.star a number greater 0 defining the assumed hazard ratio.
#' @param nEvents a number greater 0 defining the targeted number of events to achieve the
#' planned power.
#' @param L a number defining the timepoint for administrative censoring.
#' If no administrative censoring is planned \code{L=Inf}. Default is L=Inf.
#' @param from a number defining the minimum recruitment time.
#' @param lambda a number greater 0 defining the rate for the survival process of
#' group1. For the parametrization see \code{details}.
#' @param gamma a number greater 0 defining the overall rate for the censor process.
#' @param sigma,kappa a number greater 0 defining the shape parameter for the
#' survival- and the censor process. Only needed if distS or distC is set to weibull.
#' Default is 1 resulting in an exponential distribution.
#' @param distS,distC a character string defining the distribution of the survival-
#' and the censor process. Default is 'exponential'.
#'
#' @return a number defining the minimum recruitmenttime, given the recruitments N, rn and
#' the blinded estimation of the event rates from the supplied data.frame.
#' @export
#'
#' @examples
#' N <- matrix(rep(50,20),ncol=2)
#' tn <- 0:(NROW(N)-1)
#' dfSurv <- simSurvData(N=N, lambda=exp(.7)/12, sigma=1, theta=.7, gamma=exp(.8)/12, kappa=1, distS = "exponential",distC = "exponential",L=10)
#' btEstimate(dfSurv=dfSurv, theta.star=.7, N=N, tn=tn, nEvents=40, L=20,
#'            sigma=1, distS="exponential",
#'            gamma=-log(.3)/24, kappa=1, distC="exponential")
btEstimate <- function(dfSurv, theta.star, N, tn, M, tm, nEvents, L=L,sigma=sigma, distS=distS,
                       gamma=gamma, kappa=kappa, distC=distC){
  lambda <- BSSRed::bParEstimate(dfSurv=dfSurv, theta=theta.star)

  ### fix all of this
  if (missing(M)){
    return(BSSRed::rEstimate(N=N, tn=tn, theta=theta.star, nEvents=nEvents,
                             L=L, from = NROW(N),
                             lambda=lambda$lambdaP, sigma=sigma, distS=distS,
                             gamma=gamma, kappa=kappa, distC=distC))
  }
  return(BSSRed::rEstimate(N=N, tn=tn, M=M, tm=tm, theta=theta.star, nEvents=nEvents,
                           L=L, from = NROW(N),
                           lambda=lambda$lambdaP, sigma=sigma, distS=distS,
                           gamma=gamma, kappa=kappa, distC=distC))
}
