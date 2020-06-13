# depends on: -

# N = recruitment either a numberm vector or 2-col vector
# k = allokation parameter
# lambda = scale parameter in h0()
# sigma = shape parameter in h0()
# theta = fixed effect parameter (log(HR))
# gamma = rate parameter of the exponential distribution of the censoring process
# kappa = shape parameter of censor process
# L administrative censor

#' @title Simulation of a Survival Data.Frame
#' @description Simulates a two group data.frame including time to event, recruiting time and status.
#'
#' @details !!!!!!!!!!!!!!!!!!!!11how the data is simulated!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#'
#' @param N a two column matrix containing the size of both recruitment groups per
#' row, representing the number of recruited patients per group and timepoint (\code{tn}).
#' @param tn a vecotr of length = \code{NROW(N)} defining the timepoint of recruiting. Default is \code{0:(NROW(N)-1)}.
#' @param L a number defining the timepoint for administrative censoring.
#' If no administrative censoring is planned \code{L=Inf}. Default is L=Inf.
#' @param lambda a number greater 0 defining the hazard rate for the survival process of
#' group1. For the parametrization see \code{details}.
#' @param gamma a number greater 0 defining the overall hazard rate for the censor process.
#' @param sigma,kappa a number greater 0 defining the shape parameter for the
#' survival- and the censor process. Only needed if distS or distC is set to weibull.
#' Default is 1 resulting in an exponential distribution.
#' @param distS,distC a character string defining the distribution of the survival-
#' and the censor process. Default is 'exponential'.
#'
#' @return a data.frame including id (\code{id}), time to event (\code{time}), recruiting timepoint (\code{tn}),
#' grouping (\code{group}) and status (\code{status} (1 = event occured)).
#' @export
#'
#' @examples
#' N <- matrix(rep(50,10),ncol=2)
#' simSurvData(N=N, lambda=exp(.7)/12, sigma=1, theta=.7, gamma=exp(.3)/12, kappa=1, distS = "exponential",distC = "exponential",L=10)
simSurvData <- function(N, tn, lambda, sigma, theta, gamma,kappa,distS = "exponential",distC = "exponential",L){
  n1 <- sum(N[,1])
  n2 <- sum(N[,2])

  if (missing(tn)){
    tn <- 0:(NROW(N)-1)
  }
  rec_time <- c(rep(tn,N[,1]),rep(tn,N[,2]))
  group <- c(rep(0,n1),rep(1,n2))

  # Weibull event times
  if (distS == "weibull"){
    Tlat <- rweibull(sum(N),shape=sigma,scale= (theta^(group) * lambda)^(-1/sigma))
  }
  if (distS == "exponential"){
    Tlat <- rexp(sum(N),rate= theta^(group) * lambda)
  }

  # censor event times
  if (distC == "weibull")
    C <- rweibull(sum(N),shape=kappa,scale=gamma^(-1/kappa))
  if (distC == "exponential")
    C <- rexp(n=sum(N), rate=gamma)

  # follow-up times and event indicators
  if (!missing(L)){
    L=L-rec_time
    time <- pmin(Tlat, C, L)
    status <- as.numeric((Tlat < C)&(Tlat < L))
  } else {
    time <- pmin(Tlat, C)
    status <- as.numeric(Tlat < C)
  }
  # data set
  out <- data.frame(id=1:sum(N),
             tn=rec_time,
             time=time,
             status=status,
             group=group)
  class(out) <- c("dfSim","data.frame")
  out
}
