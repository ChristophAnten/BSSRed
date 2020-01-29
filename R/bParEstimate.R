# dependend on: -

#' @title Blinded Parameter Estimation
#' @description Performs a blinded parameter estimation.
#'
#' @param dfSurv a data.frame containing the variables:
#' \itemize{
#'   \item status: (0,1) := (no event, event)
#'   \item time: time to status
#'   \item group: grouping variable with two levels
#' }
#' @param theta a number defining the assumed hazard ratio.
#' @param k a number defining the balance between the two groups.
#' If not supplied it is calculated from the data.
#'
#' @return a data.frame containing the overall event rate as well as the event rates per group.
#' @export
#'
#' @examples
#' N <- matrix(rep(50,20),ncol=2)
#' df.sim <- simSurvData(N=N, lambda=-log(.7)/24, theta=.7, gamma=-log(.3)/24, distS="exponential",distC="exponential", L=20)
#' bParEstimate(dfSurv=df.sim, theta=.7, k=1)
bParEstimate <- function(dfSurv,theta,k){
  if(missing(k))
    k = sum(dfSurv$group==unique(dfSurv$group)[2])/sum(dfSurv$group==unique(dfSurv$group)[1])
  out = data.frame(lambda_ = sum(dfSurv$status)/sum(dfSurv$time))
  out$lambdaT = (1+k) / (k+1/theta) * out$lambda_
  out$lambdaP = out$lambdaT / theta
  return(out)
}
