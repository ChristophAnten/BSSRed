#' Plots a event based data frame
#'
#' @param dfSim a data.frame with time status and recruitmenttime = tn.
#' @param type grafical parameter for the plottype. Default and recommended is "l".
#' @param lty grafical parameter for the linetype. Default and recommended is 1.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
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
