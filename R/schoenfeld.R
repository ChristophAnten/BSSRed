# dependencies: -

#' @title Number of events or power calculation with the 'Schnoenfeld Formula'
#' @description Calculates the number of events needed to achieve a targeted
#' power (or the power given the number of events) in a two group scenario for all combinations of provided
#' hazard ratios \code{theta}, balance between the groups \code{k} and the alpha niveaus \code{alpha}.
#'
#' @name schoenfeld
#' @aliases pschoenfeld
#'
#' @param theta vecor of hazard ratios.
#' @param nEvent vecor of number of events.
#' @param k positive vektor. Defines the relative group sizes with the relation k:1.
#' If k=1, equal group sizes are assumed (Note: that the formula is symmetric around k=1).
#' @param alpha number between 0 and 1 defining the deserved type I error rate.
#'  Default is 0.05.
#' @param beta number between 0 and 1 defining the deserved type II error rate.
#' Default is 0.1 resulting in a power of 1-0.1=0.9 (90\%).
#' @param alternative a string; either \code{'one-sided'} or \code{'two-sided'}. Default is
#'  \code{'two-sided'}.
#'
#' @return A data frame containing the input values and calculated number of events needed to satisfy the given assumptions or
#' the power achieved with the given number of events.
#' @export
#'
#' @references Schoenfeld, D. (1983). \emph{Sample-Size Formula for the Proportional-Hazards Regression Model}. Biometrics, 39(2), 499-503. doi:10.2307/2531021
#'
#' @examples
#' # for a assumed hazard ratio of theta = .7
#' # in a design with a treatment to control group allocation of 2:1 resulting in k = 2
#' # power of 1-beta=90%
#' # and a two sided significance niveau of alpha = .05
#' x <- nschoenfeld(theta = .7, k=2, alpha=.05, beta=.1, alternative = "two-sided")
#'
#'
#' # with 20 events less the power would be
#' pschoenfeld(theta=.7, nEvents=x$nEvents-20, k=2, alpha=.05, alternative = "two-sided")
nschoenfeld <- function(theta,k=1,alpha=.05,beta=.2,alternative = "two-sided"){
  if (min(theta) <= 0){
    stop("The assumed effect size 'theta' must be greater than 0!")
  }
  if ((max(k)==Inf)|(min(k)<=0)){
    stop("The allokation parameter k must be greater 0 and smaller than Inf!")
  }
  if ((max(alpha)>=1)|(min(alpha)<=0)){
    stop("The allokation parameter k must be greater 0 and smaller than Inf!")
  }
  if ((max(beta)>=1)|(min(beta)<=0)){
    stop("The allokation parameter k must be greater 0 and smaller than Inf!")
  }
  if (sum(1-(alternative %in% c("two-sided","one-sided")))>0){
    stop("The alternative must be one of 'two-sided' or 'one-sided' (default)!")
  }

  par <- expand.grid(alternative=alternative, beta=beta, alpha=alpha, k=k, theta=theta)[,5:1]

  return(
    cbind(par,
          nEvents = (1+par$k)^2 / (par$k * log(par$theta)^2) *
            (qnorm(1-par$alpha/2^(par$alternative=="two-sided")) + qnorm(1-par$beta))^2)
  )
}

#' @rdname schoenfeld
#' @export
pschoenfeld <- function(theta,nEvents,k=1,alpha=.05,alternative = "two-sided"){
  if (min(theta) <= 0){
    stop("The assumed effect size 'theta' must be greater than 0!")
  }
  if (min(nEvents)<1){
    stop("The nuber of Events must be greater equal 1.")
  }
  if ((max(k)==Inf)|(min(k)<=0)){
    stop("The allokation parameter k must be greater 0 and smaller than Inf!")
  }
  if ((max(alpha)>=1)|(min(alpha)<=0)){
    stop("The allokation parameter k must be greater 0 and smaller than Inf!")
  }
  if (sum(1-(alternative %in% c("two-sided","one-sided")))>0){
    stop("The alternative must be one of 'two-sided' or 'one-sided' (default)!")
  }

  par <- expand.grid(alternative=alternative, alpha=alpha, k=k, nEvents=nEvents, theta=theta)[,5:1]

  return(
    cbind(par,
          power = 1-pnorm((sqrt(par$nEvents*par$k)*log(par$theta)/(1+par$k) +
                               qnorm(1-par$alpha/2^(par$alternative=="two-sided")))))
  )
  return(
    1-pnorm((sqrt(nEvent*k)*log(theta)/(1+k) + qnorm(1-alpha)))
  )
}
