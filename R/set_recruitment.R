# depends on: -

#' @title generate recruitment table
#'
#' @description A function to build a recruitment dataframe by some general conditions e.g. targeted total number of recruitments per group
#' and recruitment timePoints.
#'
#' @param targetN a 1 or 2 length vector defining the targeted total number of recruitments.
#' If \code{targetN} is a two-length vector the total number recruitments per group is defined,
#' where the treatment group is announced first.
#' @param timePoints a vector defining the timepoint of each recruitment or a positive number defining the
#' amount of recruitments (needs \code{time} to be set).
#' @param fn a string defining the shape of the cummulative recruitment. Default is set to "linear".
#' @param k a positive number defining the k:1 treatment allocation (treatment:control)
#' @param timeLimits a two length vector defining the timepoint of the first and the last recruitment.
#' Optional if \code{timePoints} is a vector defining each timepoint of recruitment.
#'
#' @return a \code{data.frame} containing the number of recruitments of the treament and control group \code{"T"} and \code{"C"} per timepoint \code{time}.
#' @export
#'
#' @examples
#' # Two different approches to define a recruitment of 1071 patients
#' # with an allokation of 2:1 (treatment:control), where 24 batches of
#' # recruits are wanted. The recruitment is assumed to take place at the
#' # beginning of each month over 24 month.
#'
#' set_recruitment(targetN=c(714,357), timePoints=0:23)
#' set_recruitment(targetN=1071, k=2, timePoints=24, timeLimits=c(0,23))
set_recruitment <- function(targetN,timePoints,fn="linear",...){

  k <- list(...)$k
  timeLimits <- list(...)$timeLimits

  if (length(targetN)==1){
    if (is.null(k))
      stop("'targetN' has only one entry and therefore 'k' needs to be defined!")
    if (k<=0)
      stop("'k' must be a positove number!")
    nT <- round(k/(k+1)*targetN)
    nC <- targetN-nT
  } else {
    if (length(targetN)>2)
      warning(sprintf("'targetN' is of length %i and therfore only the first two entrys will be used!"),length(targetN))
    nT <- targetN[1]
    nC <- targetN[2]
  }


  if (length(timePoints)==1){
    if (is.null(timeLimits))
      stop("'timePoints'is of length 1 and therefore 'timeLimits' must be defined!")
    if (length(timeLimits)==1)
      stop("'timeLimits' must be of length 2!")
    if (length(timeLimits)>2)
      warning(sprintf("'timeLimits' is of length %i and therfore only the first two entrys will be used!"),length(timeLimits))
    rec_time <- seq(timeLimits[1],timeLimits[2],length.out = timePoints)
  } else {
    rec_time <- timePoints
  }

  nBreaks <- length(rec_time)

  if(fn=="linear"){
    inc <- t(c(nT,nC)/nBreaks)
    r = floor(rep(1,nBreaks) %x% inc)
    cs = inc %% 1
    for(i in 1:nBreaks){
      r[i,] = r[i,] + (cs > 0.5)
      cs = cs+(inc %% 1) - (cs > 0.5)
    }
  }
  return(data.frame("T" = r[,1],
                    "C" = r[,2],
                    "time" = rec_time))
}

