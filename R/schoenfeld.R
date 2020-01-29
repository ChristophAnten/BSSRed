#' @title Schnoenfeld Formula
#' @description Calculates the number of events needed to achieve a targeted
#' power (or the power given the number of events) in a two group scenario with
#' the hazard ratio \code{theta}, balance between the groups \code{k} and the alpha niveau \code{alpha}.
#' @name schoenfeld
#'
#' @param theta vecor of hazard ratio theta.
#' @param nEvent vecor of number of events.
#' @param k number >0. Defines the relative group sizes with the relation k:1.
#' If k=1, equal group sizes are assumed (Note: that the formula is symmetric around k=1).
#' @param alpha number between 0 and 1 defining the deserved type I error rate.
#'  Default is 0.05.
#' @param beta number between 0 and 1 defining the deserved type II error rate.
#' Default is 0.1 resulting in a power of 1-0.1=0.9 (90\%).
#' @param alternative a string; either \code{'one-sided'} or \code{'two-sided'}. Default is
#'  \code{'two-sided'}.
#'
#' @return Returns the number of events needed to satisfy the given assumptions.
#' @export
#'
#' @examples
#' # for a assumed hazard ration of theta = .7
#' # in a design with a treatment to control group allocation of 2:1 resulting in k = 2
#' # a power of 1-beta=80%
#' # and a significance niveau of alpha = .025
#' nschoenfeld(theta = .7,k=2,alpha=.05,beta=.1,alternative = "two-sided")
#'
#' @aliases pschoenfeld
nschoenfeld <- function(theta,k=1,alpha=.05,beta=.2,alternative = "one-sided"){
  if (theta <= 0){
    stop("The assumed effect size must be greater than 0!")
  }
  if ((k==Inf)|(k<=0)){
    stop("The allokation parameter k must be greater 0 and smaller than Inf!")
  }
  if ((alpha>=1)|(alpha<=0)){
    stop("The allokation parameter k must be greater 0 and smaller than Inf!")
  }
  if ((beta>=1)|(beta<=0)){
    stop("The allokation parameter k must be greater 0 and smaller than Inf!")
  }
  if (!((alternative =="two-sided")|(alternative == "one-sided"))){
    stop("The alternative must be one of 'two-sided' or 'one-sided' (default)!")
  }
  if (alternative == "two-sided")
    alpha = alpha/2

  return(
    (1+k)^2 * (qnorm(1-alpha) + qnorm(1-beta))^2 /
      (k * log(theta)^2)
  )
}

#' @rdname schoenfeld
#' @export
pschoenfeld <- function(theta,nEvent,k=1,alpha=.05,alternative = "one-sided"){
  #stopifnot(is.numeric())
  if (theta <= 0){
    stop("The assumed effect size must be greater than 0!")
  }
  if ((k==Inf)|(k<=0)){
    stop("The allokation parameter k must be greater 0 and smaller than Inf!")
  }
  if ((alpha>=1)|(alpha<=0)){
    stop("The allokation parameter k must be greater 0 and smaller than Inf!")
  }
  if (sum(nEvent<1)>0){
    stop("The nuber of Events must be greater equal 1.")
  }
  if (!((alternative =="two-sided")|(alternative == "one-sided"))){
    stop("The alternative must be one of 'two-sided' or 'one-sided' (default)!")
  }
  if (alternative == "two-sided")
    alpha = alpha/2

  return(
    1-pnorm((sqrt(nEvent*k)*log(theta)/(1+k) + qnorm(1-alpha)))
  )
}
