# depends on: nschoenfeld, pschoenfeld, calc_expEvents, rEstimate

#' Title
#'
#' @param recruitment
#' @param rec_time
#' @param study_end
#' @param lambdaC a number defining the hazard rate of the control group.
#' @param lambdaT a number defining the hazard rate of the treatment group.  Alternatively \code{theta} can be defined.
#' @param theta a number defining the hazard ratio. Alternatively \code{lambdaT} can be defined.
#' @param gamma a number defining the hazard of the censor process.
#' @param dist_surv a numeric defining the distribution of the survival process.
#' @param dist_censor
#' @param alpha
#' @param beta
#' @param alternative
#'
#' @return a list of class \code{studyPlan}
#' @export
#'
#' @examples
#' set_studyPlan(recruitment=set_recruitment(c(2400,1200),1:20),
#'               rec_time = 10,
#'               study_end = 22,
#'               lambdaC = get_hazardRate(.3,24),
#'               theta = .7,
#'               gamma = get_hazardRate(.2,24),
#'               dist_surv = "exponential",
#'               dist_censor = "exponential",beta=.1)


recruitment=set_recruitment(c(2400,1200),1:20);
              rec_time = 10;
              study_end = 22;
              lambdaC = get_hazardRate(.3,24);
              theta = .7;
              gamma = get_hazardRate(.2,24);
              dist_surv = "exponential";
              dist_censor = "exponential";
              beta=.1;
              alpha=.05;
              alternative="two-sided"
# eine größe weglassen: study_end, rec_time, power, recruitment?
# wenn alle gegeben sind berechnet alle optionen um das Modell einzuhalten.
# recruitment
set_studyPlan <- function(recruitment,rec_time,study_end,
                          lambdaC,theta,gamma,dist_surv,dist_censor,
                          alpha=.05,beta=.2,alternative="two-sided",...){
  if (missing(rec_time))
    rec_time = NULL

  if (!(NCOL(recruitment)==3))
    stop("'recruitmet' must be a three-column matrix or data.frame, where the first two columns consists of the number of recruitmens per recruitment step and the third column defines the timepoint at which the recruitment takes place.")

  k <- sum(recruitment[,1])/sum(recruitment[,2])

  lambdaT <- list(...)$lambdaT
  sigma <- list(...)$sigma
  kappa <- list(...)$kappa

  if (missing(theta))
    theta = lambdaT/lambdaC

  if (is.null(lambdaT)){
    lambdaT = lambdaC*theta
  }
  # print this shit
  nEvent <- nschoenfeld(theta=theta,k=k,alpha=alpha,beta=beta,alternative=alternative)$nEvents

  if (!(sum(missing(rec_time),
            missing(study_end),
            missing(beta)) %in% 1:2))
    stop("At least two of the parameters rec_time, study_end and beta must be defined.")

  # estimate power
  # estimate study length
  # estimate rec_time

  # estimate all three
  if (missing(beta)){
    expEvents_planned <- calc_expEvents(N=recruitment[1:rec_time,1:2], theta, L=study_end, tn=recruitment[1:rec_time,3],
                                        lambda=lambdaC, sigma=sigma, distS=dist_surv,
                                        gamma, kappa=1, distC=dist_censor)
    expPower_planned <- pschoenfeld(theta = theta,nEvent = sum(expEvents_planned), k=k, alpha=alpha, alternative=alternative)

  }
  if (missing(study_end)){
    study_end <- rec_time
      study_end <- rec_time+1
      nEvent <- nschoenfeld(theta=theta,k=k,alpha=alpha,beta=beta,alternative=alternative)$nEvents
      fn <- function(x,N,theta,tn,lambda,sigma,distS,gamma,kappa,distC,nEvent){
        abs(sum(calc_expEvents(N=recruitment[1:rec_time,1:2], theta, L=x, tn=recruitment[1:rec_time,3],
                                              lambda=lambdaC, sigma=sigma, distS=dist_surv,
                                              gamma, kappa=1, distC=dist_censor))-nEvent)
      }
      optimize(fn,interval=c(rec_time,rec_time*10),N=recruitment[1:rec_time,1:2], theta=theta, tn=recruitment[1:rec_time,3],
               lambda=lambdaC, sigma=sigma, distS=dist_surv,
               gamma=gamma, kappa=kappa, distC=dist_censor, nEvent=nEvent)
  }
  if (missing(rec_time)){
    nEvent <- nschoenfeld(theta=theta,k=k,alpha=alpha,beta=beta,alternative=alternative)$nEvents
    r <- rEstimate(N=recruitment[,1:2],tn=recruitment[,3],theta= theta,nEvents=nEvent,L=study_end,from = 1,
                   lambda=lambdaC, sigma=sigma, distS=dist_surv,
                   gamma=gamma, kappa=kappa, distC=dist_censor)
    expEvents_est <- calc_expEvents(N=recruitment[1:r$rec_time,1:2], theta, L=study_end, tn=recruitment[1:r$rec_time,3],
                                    lambda=lambdaC, sigma=sigma, distS=dist_surv,
                                    gamma, kappa=1, distC=dist_censor)
    expPower_est <- pschoenfeld(theta = theta,nEvent = sum(expEvents_est), k=k, alpha=alpha, alternative=alternative)
  }



  studyPlan <- list(alpha = alpha,
                    beta = beta,
                    alternative = alternative,
                    study_end = study_end,
                    k = k,
                    theta=theta,
                    recruitment = recruitment,
                    nEvents = nEvents,
                    par = list(surv = list(lambdaT=lambdaT,
                                           lambdaC=lambdaC,
                                           sigma=sigma,
                                           dist=dist_surv),
                               censor = list(gamma=gamma,
                                             kappa=kappa,
                                             dist=dist_censor)))


  if (!is.null(rec_time)){
    expEvents_planned <- calc_expEvents(N=recruitment[1:rec_time,1:2], theta, L=study_end, tn=recruitment[1:rec_time,3],
                                        lambda=lambdaC, sigma=sigma, distS=dist_surv,
                                        gamma, kappa=1, distC=dist_censor)
    expPower_planned <- pschoenfeld(theta = theta,nEvent = sum(expEvents_planned), k=k, alpha=alpha, alternative=alternative)

    planned <- list(rec_time = rec_time,
                    expEvents = expEvents_planned,
                    expPower = expPower_planned)
    studyPlan <- append(studyPlan,
                        list(planned=planned))
  }

  r <- rEstimate(N=recruitment[,1:2],tn=recruitment[,3],theta= theta,nEvents=nEvent,L=study_end,from = 1,
                 lambda=lambdaC, sigma=sigma, distS=dist_surv,
                 gamma=gamma, kappa=kappa, distC=dist_censor)
  expEvents_est <- calc_expEvents(N=recruitment[1:r$rec_time,1:2], theta, L=study_end, tn=recruitment[1:r$rec_time,3],
                                  lambda=lambdaC, sigma=sigma, distS=dist_surv,
                                  gamma, kappa=1, distC=dist_censor)
  expPower_est <- pschoenfeld(theta = theta,nEvent = sum(expEvents_est), k=k, alpha=alpha, alternative=alternative)
  estimated <- list(rEstimate = r,
                    expEvents = expEvents_est,
                    expPower = expPower_est)
  studyPlan <- append(studyPlan,
                      list(estimated=estimated))

  class(studyPlan) <- c(class(studyPlan),"studyPlan")
  studyPlan
}
# set_studyPlan(recruitment=set_recruitment(c(2400,1200),1:20),
#                              rec_time = 10,
#                              study_end = 22,
#                            lambdaC = get_hazardRate(.3,24),
#                              theta = .7,
#                              gamma = get_hazardRate(.2,24),
#                              dist_surv = "exponential",
#                              dist_censor = "exponential",beta=.1)
print.studyPlan <- function(studyPlan){
  cat("Initial number of event calculation with the schoenfeld formula:\n")
  print(data.frame(alpha = studyPlan$alpha,
                   beta = studyPlan$beta,
                   k = studyPlan$k,
                   theta = studyPlan$theta,
                   alternative = studyPlan$alternative,
                   nEvents = studyPlan$nEvents), row.names = FALSE)

  cat("\nAssumed hazard rates for the survival processes:\n")
  out <- data.frame(1)
  for(i in names(studyPlan$par$surv)){
    if (!is.null(studyPlan$par$surv[[i]])){
      tmp <- data.frame(studyPlan$par$surv[[i]])
      colnames(tmp) <- i
      out <- cbind(out,tmp)

    }
  }
  print(out[,-1], row.names = FALSE)

  cat("\nAssumed hazard rates for the censor processes:\n")
  out <- data.frame(1)
  for(i in names(studyPlan$par$censor)){
    if (!is.null(studyPlan$par$censor[[i]])){
      tmp <- data.frame(studyPlan$par$censor[[i]])
      colnames(tmp) <- i
      out <- cbind(out,tmp)

    }
  }
  print(out[,-1], row.names = FALSE)

  estimated <- data.frame(nbatches = studyPlan$estimated$rEstimate$rec_batch,
                          rec_nT = studyPlan$estimated$rEstimate$rec_total_number[1],
                          rec_nC = studyPlan$estimated$rEstimate$rec_total_number[2],
                          expEvents = sum(studyPlan$estimated$expEvents),
                          expPower = studyPlan$estimated$expPower$power)
  rownames(estimated) <- "estimated"

  if(!is.null(studyPlan$planned)){
    planned <- data.frame(nbatches = studyPlan$planned$rec_time,
                          rec_nT = sum(studyPlan$recruitment[1:studyPlan$planned$rec_time,1]),
                          rec_nC = sum(studyPlan$recruitment[1:studyPlan$planned$rec_time,2]),
                          expEvents = sum(studyPlan$planned$expEvents),
                          expPower = studyPlan$planned$expPower$power)
    rownames(planned) <- "planned"

    cat("\nRecruitments and expected number of events and power:\n")
    print(rbind(planned,estimated))
  } else {
    cat("\nEstimated recruitments and expected number of events and power:\n")
    print(estimated, row.names = FALSE)
  }

  invisible(studyPlan)
}
#print(studyPlan)

# output #
# needed events
# exp calc events
# exp calc needed rec time
# allow weibull
# in par
