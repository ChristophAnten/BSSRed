#' Plots a event based data frame
#'
#' @param dfSim a data.frame with time status and recruitmenttime = tn.
#' @param type grafical parameter for the plottype. Default and recommended is "l".
#' @param lty grafical parameter for the linetype. Default and recommended is 1.
#' @param ... additional \code{matplot()} related arguments.
#'
#' @return a plot of the supplied data.frame.
#' @export
#'
#' @examples
#' N <- matrix(rep(50,10),ncol=2)
#' dfSim <- simSurvData(N=N, lambda=-log(1-.3)/24, sigma=1, theta=.7, gamma=-log(1-.2)/24, kappa=1, distS = "exponential",distC = "exponential",L=10)
#' plot_dfSim(dfSim)
plot_dfSim <- function(dfSim,type="l",lty=1,...){
  dfSim$time.end = dfSim$time+dfSim$tn
  index <- with(dfSim,order(group,tn,-status,time.end))
  time <- t(dfSim[index,c("tn","time.end")])
  count <- matrix(rep(1:length(index),each=2),nrow=2)
  status <- dfSim$status[index]
  matplot(x=time,
          y=count,
          type=type,
          col=ifelse(status==0,"darkgreen","blue"),
          lty=lty,...)
}
