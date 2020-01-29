# depends on: calc_expEvents()

#--- assumed recruits ---#
# r: at start may be a vector
# rn: infinite continuition
# k: allocation parameter
#--- assumed parameters ---#
# theta: HR
# lambda: event rate
# sigma: event rate shape parameter
# distS: distribution of event process
# gamma: censor rate
# kappa: censor rate shape parameter
# distC: distribution of censor process
#--- wished model ---#
# alpha: signifikanz niveau
# beta: beta niveau
# alternative: # one- or two-sided
#--- needed events for Power and alpha ---#
# nEvents <- sse.schoenfeld()
#--- needed time to reach Power by given recruitments and nuicence parameters ---#

#' @title Recruitment Time Estimation
#' @description Calculates the minimum recruitment time and number of patients needed to achieve the targeted
#' number of events befor administrative study end.
#'
#' @param N a two column matrix containing the size of both recruitment groups per
#' row, representing the number of recruited patients per group and timepoint (\code{t}).
#' @param tn a vecotr of length = \code{NROW(N)} defining the timepoint of recruiting.
#' @param M a two column matrix containing the size of both recruitment groups per
#' row, representing the number of future recruited patients per group and timepoint (\code{tm}).
#' If not defined it will be taken from the last row of \code{N}.
#' @param tm a vecotr of length = \code{NROW(M)} defining the timepoint of future recruiting.
#' If not defined it will be interpolated from the last two entries of \code{tn}.
#' @param theta a number greater 0 defining the assumed hazard ratio.
#' @param nEvents a number greater 0 defining the targeted number of events to achieve the
#' planned power.
#' @param L a number defining the timepoint for administrative censoring.
#' If no administrative censoring is planned \code{L=Inf}. Default is L=Inf.
#' @param from a number > 0 defining the minimum number of recruitment steps. Default is 1,
#' starting with the first recruited batch of patients.
#' @param lambda a number greater 0 defining the hazard rate for the survival process of
#' group1. For the parametrization see \code{details}.
#' @param gamma a number greater 0 defining the overall hazard rate for the censor process.
#' @param sigma,kappa a number greater 0 defining the shape parameter for the
#' survival- and the censor process. Only needed if distS or distC is set to weibull.
#' Default is 1 resulting in an exponential distribution.
#' @param distS,distC a character string defining the distribution of the survival-
#' and the censor process. Default is 'exponential'.
#'
#' @return a list with the minimum needed batches of recruits defined by N (\code{rec_batch}),
#' the minimum time for recruitment dependend on tn (\code{rec_time}) and the minimum amount of
#' patients (\code{rec_total_number}).
#' @export
#'
#' @examples
#' N <- matrix(rep(50,20),ncol=2)
#' tn <- 0:(NROW(N)-1)
#' rEstimate(N=N, tn=tn, theta=.7, nEvents=40, L=20, from = 1,
#'           lambda=-log(.7)/24, sigma=1, distS="exponential",
#'           gamma=-log(.3)/24, kappa=1, distC="exponential")
rEstimate <- function(N,tn,M,tm,theta,nEvents,L,from = 1,
                      lambda, sigma=1, distS="exponential",
                      gamma, kappa=1, distC="exponential"){
  # N <- par.list$N
  # M <- par.list$M
  # tn <- par.list$tn
  # tm <- par.list$tm
  if (missing(M)){
    M = N[NROW(N),]
  }
  if (missing(tm)){
    if (length(tn) > 1){
      a <- tn[length(tn)]-tn[length(tn)-1]
    } else {
      a <- 1
    }
    tm <- 1:NROW(M)*a + max(tn)
  }
  recruitments <- rbind(N,M)
  rec.time <- c(tn,tm)
  inf.recruitments <- recruitments[NROW(recruitments),]
  if (length(tm)==1){
    inf.rec.time <- (tm - tn[length(tn)])
  } else {
    inf.rec.time <- tm[length(tm)] - tm[length(tm)-1]
  }
  i=1
  t = min(NROW(recruitments),from)
  expEvents <- numeric(L-t+1)
  while(TRUE){
    expEvents[i] <- sum(calc_expEvents(recruitments[1:t,,drop=FALSE], L=L, tn = rec.time[1:t], theta=theta,
                                       lambda=lambda, sigma=sigma, distS=distS,
                                       gamma=gamma, kappa=kappa, distC=distC))
    if (t==L){
      if (nEvents > expEvents[i]){
        warning(sprintf("With the specified recruitments the targeted number of %.0f events can not be achieved before administrative censoring. At the timepoint L= %.0f only %.2f events will be expected to occur.",nEvents,L,expEvents[i]))
      }
      return(list(rec_batch = t,
                  rec_time = rec.time[t],
                  rec_total_number = colSums(recruitments[1:t,,drop=FALSE])))
    }
    if(nEvents < expEvents[i]){
      return(list(rec_batch = t,
                  rec_time = rec.time[t],
                  rec_total_number = colSums(recruitments[1:t,,drop=FALSE])))
    }
    if (NROW(recruitments)<=t){
      recruitments = rbind(recruitments,inf.recruitments)
      rec.time = c(rec.time,max(rec.time) + inf.rec.time)
    }
    t = t+1
    i=i+1
  }
  return("from was set too high")
}
