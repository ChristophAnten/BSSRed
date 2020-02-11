
#' Simulation
#'
#' @usage
#' sim_BSSRed(N, tn, M, tm, addRecT, sim.lambdaP, sigma=1, theta.star, gamma, kappa=1,
#'            distS=c("exponential","weibull"), distC=c("exponential","weibull"),
#'            nEvents, L, nSim=100, par.est=FALSE, fixed.min.studytime=TRUE, seed=235711)
#'
#' @param N a two-column matrix which represents the size of recruitment per group.
#' @param tn a vector of length \code{NROW(N)} defining the timepoint of each recruitment-batch.
#' @param M a matrix defining the ongoining recruitment after the planned N. Default is the
#' last row of N.
#' @param tm a vector of length \code{NROW(M)} defining the timepoint of each future
#' recruitment-batch. Default is the next equidistant timepoint after the last of \code{tn}.
#' @param addRecT a vector defining the allowed additional recruitment batches.
#' @param sim.lambdaP a vector defining the hazard rates of the control group.
#' @param sigma a number defining the shape parameter of the survival prozess.
#' @param theta.star a number defining the true hazard ratio.
#' @param gamma a number defining the hazard of the censor process.
#' @param kappa a number defining the shape parameter of the censor prozess.
#' @param distS a string defining the distribution of the survival process. Default is \code{"exponential"}.
#' @param distC a string defining the distribution of the censor process. Default is \code{"exponential"}.
#' @param nEvents a number defining the targeted number of events. (alternatively alpha, beta)
#' @param L a number defining the administrative censoring time.
#' @param nSim a number defining the number of simulations per parametrisation.
#' @param par.est a boolean. If \code{TRUE} the parameter-estimates and teststatistik will be calculated. Warning: will greatliy increase runtime.
#' @param fixed.min.studytime a boolean. If \code{TRUE} the simulated studies will not stop before \code{L}. If \code{FALSE} the study will end when the targeted amount of events is reached.
#' @param seed
#'
#' @return a list containing the results of each simulation as a data.frame.
#' @export
#'
#' @examples
sim_BSSRed <- function(N,tn,M,tm,addRecT,sim.lambdaP, sigma,
                       theta.star,gamma,kappa,distS,distC,nEvents,
                       L,nSim=100, par.est=FALSE, fixed.min.studytime=TRUE, seed=235711){

  # lambdaP.vec <- -log(seq(.6,.8,length.out = 9))/24
  res <- list()
  nRep <- length(addRecT)
  nPar <- length(sim.lambdaP)
  cat(sprintf("Start simulation with nLambda = %i, nsim = %i and nLimits = %i.\n",nPar,nSim,nRep))
  cat(sprintf("Total number of iterations = %i.\n",nPar*nSim*nRep))
  res$results <- data.frame(
    lambda = rep(sim.lambdaP,each=nRep*nSim),
    sim = rep(rep(1:nSim,each=nRep),nPar),
    type = rep(addRecT,nSim*nPar),
    recT = numeric(nRep*nSim*nPar),
    loglik = numeric(nRep*nSim*nPar),
    waldtest = numeric(nRep*nSim*nPar),
    sctest = numeric(nRep*nSim*nPar),
    study.end = numeric(nRep*nSim*nPar),
    events = numeric(nRep*nSim*nPar),
    nPop = numeric(nRep*nSim*nPar),
    theta.est = numeric(nRep*nSim*nPar))
  # l=2;i=1;j=1;
  # l=1;i=1;j=1;

  #l=nPar;i=nSim;j=nRep;
  for ( l in 1:nPar){
    set.seed(seed)
    cat(sprintf("\nSim %i of %i parameters [lambdaP = %1.4f]\n",l,nPar,sim.lambdaP[l]))
    MM <- M
    tmm <- tm
    for (i in 1:nSim){
      if ((i %% (nSim/10)) == 0)
        cat(paste("Sim",i,"of",nSim,"\n"))
      if (NROW(M) < max(addRecT)){
        MM <- rbind(M,M[NROW(M),,drop=FALSE] %x% rep(1,max(addRecT)-NROW(M)))
        if (length(tm)==1){
          tmm <- max(tn) + (tm - tn[length(tn)]) * 1:max(addRecT)
        } else {
          tmm <- max(tm) + (tm[length(tm)] - tm[length(tm)-1]) * 1:max(addRecT)
        }
      }
      NN <- rbind(N,MM)
      tnn <- c(tn,tmm)
      df.sim <- BSSRed::simSurvData(N=NN,
                                    tn = tnn,
                                    lambda=sim.lambdaP[l],
                                    sigma=sigma,
                                    theta=theta.star,
                                    gamma=gamma,
                                    kappa=kappa,
                                    distS = distS,
                                    distC = distC,
                                    L=Inf)
      df.sim$total.time = df.sim$time+df.sim$tn
      #plot_simData(df.sim)
      for (j in 1:nRep){

        recT=length(tn) # always checks at end of tn (no 'from' statement)
        while((recT < BSSRed::btEstimate(dfSurv = BSSRed::cut_dfSim(df.sim,tnn[recT]),theta.star=theta.star,
                                         N=N, tn = tn, M=M,
                                         tm=tm, L=L,nEvents = nEvents,
                                         sigma=sigma, distS=distS,
                                         gamma=gamma, kappa=kappa, distC=distC)$rec_batch) &
              (recT<(length(tn)+addRecT[j]))){
          recT = recT+1
        }

        # get the right timepoint
        df.sim.fix <- df.sim[df.sim$tn<=tnn[recT],]
        df.sim.fix <- df.sim.fix[order(df.sim.fix$total.time),]
        res.index <- (l-1) * nSim * nRep + (i-1) * nRep + j
        res$results$recT[res.index] = recT
        res$results$study.end[res.index] = ifelse(
          fixed.min.studytime,
          max(L,ceiling(min(df.sim.fix$total.time[which(cumsum(df.sim.fix$status)==ceiling(nEvents))]))),
          ceiling(min(df.sim.fix$total.time[which(cumsum(df.sim.fix$status)==ceiling(nEvents))])))
        res$results$events[res.index] = sum(BSSRed::cut_dfSim(df.sim.fix,res$results$study.end[res.index])$status)
        res$results$nPop[res.index] = NROW(df.sim.fix)

        if (par.est){
          cat("parameter estimation not yet implemented")
        } else {
          res$results$loglik[res.index] = NA
          res$results$waldtest[res.index] = NA
          res$results$sctest[res.index] = NA
          res$results$theta.est[res.index] = NA
        }
      } # end loop j
    } # end loop i
  } # end loop l
  return(res)
} # end function

