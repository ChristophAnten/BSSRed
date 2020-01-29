# dependend on: -

#' Cuts the Data.Frame at a specific timepoint, changing status accordingly.
#'
#' @param df.sim a data.frame containing the variables:
#' \itemize{
#'   \item status: (0,1) := (no event,event)
#'   \item time: time to status
#'   \item tn: recruitment timepoint
#' }
#' @param t a number defining the timepoint at which the dataframe is cut.
#'
#' @details the status of events occuring at a later timepoint time > t are set to zero. Patients with a
#' recruitment timepoint > t are removed.
#'
#' @return a data.frame with the same structure of the handed data.
#' @export
#'
#' @examples
#' N <- matrix(rep(50,10),ncol=2)
#' df.sim <- simSurvData(N=N, lambda=exp(.7)/12, theta=.7, gamma=exp(.3)/12, distS="exponential",distC="exponential", L=10)
#' df.sim_cut3 <- cut.dfSim(df.sim,3)
cut_dfSim <- function(df.sim,t){
  df.sim = df.sim[df.sim$tn<=t,]
  df.sim$status <- ifelse((df.sim$time+df.sim$tn)<t,df.sim$status,0)
  df.sim$time <- ifelse((df.sim$time+df.sim$tn)<t,df.sim$time,t-df.sim$tn)
  return(df.sim)
}
