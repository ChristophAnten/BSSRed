# depends on: -


#' @title Get hazard rates
#' @describeIn Functions to calculate hazard rates by the proportion of confirmed events at given time.
#'
#' @usage get_hazardRate(prop,time,theta=NULL)
#' @usage get_eventProportion(hazard,time,theta=NULL)
#'
#' @param prop a vector of numbers greater 0 and smaller 1 defining the proportion of patients having an event before 'time'.
#' @param hazard a vector of numbers greater 0 defining the hazard rates per timeintervall of length 1.
#' @param time a number defining the timeframe.
#' @param theta optional number defining the hazard rate. Returning(-log(1-prop) /time *  theta)

#' @return a vector defining the hazard rate per timeintevall of length 1.
#' @export
#'
#' @examples
#' # Assuming, that 30% of patients in the control group have an event after 24 month,
#' # the monthly hazard rate for the control group can be calculated by:
#' get_hazardRate(.3,24)
#'
#' # With a assumed reduction of the hazard by 30% in the treatment group we can calculate
#' # the hazard rate as follows:
#' get_hazardRate(.3,24,1-.3)
#' @aliases get_eventProportion
#' @name get_hazardRate
get_hazardRate <- function(prop,time,theta=NULL){
  if (is.null(theta)){
    return(-log(1-prop) / time)
  }
  if (is.numeric(theta)){
    return(-log(1-prop) /time *  theta)
  }
}

#' @rdname get_hazardRate
#' @export
get_eventProportion <- function(hazard,time,theta=NULL){
  if (is.null(theta)){
    return(1-exp(-hazard * time))
  }
  if (is.numeric(theta)){
    return(1-exp(-hazard * time * theta))
  }
}
