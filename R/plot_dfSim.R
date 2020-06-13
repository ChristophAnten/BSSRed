#' Plots an event based data frame
#'
#' @param dfSim a data.frame with time status and recruitmenttime = tn.
#' @param col a vector of two colours. The first one defines the color of
#' observations without an event and the second one the color of the observations with an event.
#' Default is c("grey","black").
#' @param ... additional \code{matplot()} related arguments.
#'
#' @details
#' The graph shows the observation time per patient from recruitment until administrative censoring.
#'
#' Patients with an event are colored differently.
#'
#' The Observations are sorted by observation time, status and group.
#'
#' @return Generates a new plot.
#' @export
#'
#' @seealso \code{\link[graphics]{matplot}} \code{\link{simSurvData}}
#'
#' @examples
#' N <- matrix(rep(50,10),ncol=2)
#' dfSim <- simSurvData(N=N, lambda=-log(1-.3)/24, sigma=1, theta=.7, gamma=-log(1-.2)/24, kappa=1, distS = "exponential",distC = "exponential",L=10)
#' plot(dfSim)
#'
#' # to see the corresponding Kaplan-Meyer curve
#' require(survival)
#' fit <- survfit(Surv(time=time-tn, event = status) ~group,data = dfSim)
#' plot(fit)
plot.dfSim <- function(dfSim,col=c("grey","black"),...){
  dfSim$time.end = dfSim$time+dfSim$tn
  index <- with(dfSim,order(group,tn,-status,time.end))
  time <- t(dfSim[index,c("tn","time.end")])
  count <- matrix(rep(1:length(index),each=2),nrow=2)
  status <- dfSim$status[index]
  matplot(x=time,
          y=count,
          type="l",
          col=ifelse(status==0,col[1],col[2]),
          lty=1,...)
}

# cols <- rainbow(3, s = 0.5)
# boxplot(study.end ~ type + lambda, data = res$results,
#         at = rep(seq(0,4*(9-1),4),each=3)+rep(0:2,9), col = cols, xaxt = "n")
# axis(side = 1, at = seq(1,4*9-1,4), labels = seq(propP-.1,propP+.1,length.out = 9))
# legend("topleft", fill = rainbow(3, s = 0.5), legend = unique(), horiz = T)
#
# boxplot(nPop ~ type + lambda, data = res$results,
#         at = rep(seq(0,4*(9-1),4),each=3)+rep(0:2,9), col = cols, xaxt = "n")
# axis(side = 1, at = seq(1,4*9-1,4), labels = seq(propP-.1,propP+.1,length.out = 9))
# legend("topleft", fill = rainbow(3, s = 0.5), legend = unique(), horiz = T)
#
# plot(y=nPop, x=lambda,
#         data = aggregate(res$results[,-(1:3)],list(lambda=res$results$lambda,type=res$results$type),mean),
#      type="l")
# axis(side = 1, at = seq(1,4*9-1,4), labels = seq(propP-.1,propP+.1,length.out = 9))
# legend("topleft", fill = rainbow(3, s = 0.5), legend = unique(), horiz = T)
#
#
# aggregate(res$results[,-(1:3)],list(lambda=res$results$lambda,type=res$results$type),mean)
