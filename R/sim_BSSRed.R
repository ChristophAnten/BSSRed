
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
#' @param seed a number defining the seed.
#'
#' @return a list containing the results of each simulation as a data.frame.
#' @export
#'
#' @examples

# recruitment <- data.frame("T"=c(seq(3,30,by=3),rep(34,5),rep(35,5)),
#                           "C"=c(seq(6,60,by=6),rep(68,5),rep(70,5)),
#                           "time"=0:19)
# class(recruitment) <- c("recruit_frame","data.frame")


# addRecT = c(0,3,6);
# lambdaC=get_hazardRate(seq(propP-.1,propP,length.out = 6),timeP)
# sigma=1;theta.star = .7;gamma=0.009297648;kappa=1;
# distS = "exponential";distC = "exponential";nEvents=374#371.6752;
# nSim=1000;par.est = F;fixed.min.studytime = FALSE;seed = 1
#
# d <- 1
# L=39;
# addRecruitment <- data.frame("T"=34,
#                              "C"=68,
#                              "time"=20+d)
# class(addRecruitment) <- c("recruit_frame","data.frame")
# recruitment <- data.frame("T"=c(seq(3,30,by=3),rep(34,5),rep(35,5)),
#                           "C"=c(seq(6,60,by=6),rep(68,5),rep(70,5)),
#                           "time"=0:19+d)
# class(recruitment) <- c("recruit_frame","data.frame")
#
# res <- sim_BSSRed(recruitment=recruitment,addRecT = addRecT,addRecruitment=addRecruitment,
#            lambdaC=lambdaC,sigma=sigma,theta.star = theta.star,gamma=gamma,kappa=kappa,
#            distS = distS,distC = distC,nEvents=nEvents,L=L,
#            nSim=nSim,par.est = par.est,fixed.min.studytime = fixed.min.studytime,seed = seed)
# resWithPlus <- res
# resWoPlus
# resWoPlusOnL
# resWithPlus

sim_BSSRed <- function(recruitment,addRecT,lambdaC, sigma,...,
                       theta.star,gamma,kappa,distS,distC,nEvents,
                       L,nSim=100, par.est=FALSE, fixed.min.studytime=FALSE, seed=235711){

  # if(!("recruit_frame" %in% class(recruitment))){
  #
  #   recruitment <- list(...)$tn
  # }

  addRecruitment <- list(...)$addRecruitment

  # res$results <- data.frame(
  #   lambda = rep(lambdaC,each=nRep*nSim),
  #   sim = rep(rep(1:nSim,each=nRep),nPar),
  #   type = rep(addRecT,nSim*nPar),
  #   recT = numeric(nRep*nSim*nPar),
  #   loglik = numeric(nRep*nSim*nPar),
  #   waldtest = numeric(nRep*nSim*nPar),
  #   sctest = numeric(nRep*nSim*nPar),
  #   study.end = numeric(nRep*nSim*nPar),
  #   events = numeric(nRep*nSim*nPar),
  #   nPop = numeric(nRep*nSim*nPar),
  #   theta.est = numeric(nRep*nSim*nPar))

  # simPar <- expand.grid(distS=distS,
  #                       lambdaC=lambdaC,
  #                       sigma=sigma,
  #                       theta.star=theta.star,
  #                       distC=distC,
  #                       gamma=gamma,
  #                       kappa=kappa,
  #                       nEvents=nEvents,
  #                       L=L,
  #                       addRecT=addRecT)
  # lambdaP.vec <- -log(seq(.6,.8,length.out = 9))/24
  res <- list()
  res$results <- expand.grid(addRecT=addRecT,
                             nEvents=nEvents,
                             L=L,
                             sim = 1:nSim,
                             distS=distS,
                             lambdaC=lambdaC,
                             sigma=sigma,
                             theta.star=theta.star,
                             distC=distC,
                             gamma=gamma,
                             kappa=kappa)
  res$results[,"recT"] <- numeric()
  res$results[,"study.end"] <- numeric()
  res$results[,"nPop"] <- numeric()
  res$results["events"] <- numeric()
  # res$results[,"loglik"] <- numeric()
  # res$results[,"waldtest"] <- numeric()
  # res$results[,"sctest"] <- numeric()
  # res$results[,"theta.est"] <- numeric()
  nRep <- length(addRecT)*length(nEvents)*length(L)
  nPar <- NROW(res$results)/nRep/nSim
  cat(sprintf("Start simulation with nSim = %i resulting in a total of %i simulations.\n",nSim,NROW(res$results)))
  cat(sprintf("Number of different Study-frames [L,addRecT,nEvents]: %i.\n",nRep))
  cat(sprintf("Number of different Parametrisations: %i.\n",nPar))

  N <- recruitment[,1:2]
  tn <- recruitment[,3]
  if(is.null(addRecruitment)){
    M <- N[NROW(N),,drop=FALSE]
    tm <- tn[NROW(N)]+1
  } else {
    M <- addRecruitment[,c("T","C")]
    tm <- addRecruitment[,"time"]
  }
  i=1;l=1;j=1;
  for (l in 1:nPar){
    set.seed(seed)
    par.index <- (l-1) * nSim * nRep+1
    cat(sprintf("\nSim %i of %i parametrisations.\n",l,nPar))
    print(res$results[par.index+(0:(nRep-1)),1:11])
    MM <- M
    tmm <- tm
    if (NROW(M) < max(addRecT)){
      MM[(NROW(M)):max(addRecT),] <- M[NROW(M),]
      if (length(tm)==1){
        tmm <- max(tn) + (tm - tn[length(tn)]) * 1:max(addRecT)
        } else {
          tmm <- max(tm) + (tm[length(tm)] - tm[length(tm)-1]) * 1:max(addRecT)
        }
    }
    NN <- rbind(N,MM)
    tnn <- c(tn,tmm)
    # actually generate new data
    for (i in 1:nSim){
      if ((i %% (nSim/10)) == 0)
        cat(paste("Sim",i,"of",nSim,"\n"))
      df.sim <- BSSRed::simSurvData(N=NN,
                                    tn = tnn,
                                    lambda=res$results$lambdaC[par.index],
                                    sigma=res$results$sigma[par.index],
                                    theta=res$results$theta.star[par.index],
                                    gamma=res$results$gamma[par.index],
                                    kappa=res$results$kappa[par.index],
                                    distS = res$results$distS[par.index],
                                    distC = res$results$distC[par.index],
                                    L=Inf)
      df.sim$total.time = df.sim$time+df.sim$tn
      #plot_simData(df.sim)
      for (j in 1:nRep){
        res.index <- (l-1) * nSim * nRep + (i-1) * nRep + j
        # print(res.index)
        # print(sprintf("l=%i, i=%i, j=%i, s=%i,                 nSim * nPar =%i, nPar=%i",l,i,j,s,nSim * nPar ,nPar))
        recT=length(tn) # always checks at end of tn (no 'from' statement)
        suppressWarnings(
        while((recT < BSSRed::btEstimate(dfSurv = BSSRed::cut_dfSim(df.sim,tnn[recT]),
                                         theta.star=res$results$theta.star[res.index],
                                         N=N, tn = tn, M=M,
                                         tm=tm, L=res$results$L[res.index],nEvents = res$results$nEvents[res.index],
                                         sigma=res$results$sigma[res.index], distS=res$results$distS[res.index],
                                         gamma=res$results$gamma[res.index], kappa=res$results$kappa[res.index],
                                         distC=res$results$distC[res.index])$rec_batch) &
              (recT<(length(tn)+res$results$addRecT[res.index]))){
          recT = recT+1
        })
        # get the right timepoint
        df.sim.fix <- df.sim[df.sim$tn<=tnn[recT],]
        df.sim.fix <- df.sim.fix[order(df.sim.fix$total.time),]
        res$results$recT[res.index] = recT

        diff <- cumsum(df.sim.fix$status)-ceiling(res$results$nEvents[res.index])
        if(diff[length(diff)]<0){
          res$results$study.end[res.index] = ceiling(max(df.sim.fix$total.time))
        } else if (fixed.min.studytime){
          res$results$study.end[res.index] = max(res$results$L[res.index],
              ceiling(min(df.sim.fix$total.time[which(diff>=0)])))
        } else {
          res$results$study.end[res.index] = ceiling(min(df.sim.fix$total.time[which(diff>=0)]))
        }
        # res$results$study.end[res.index] = ifelse(
        #   fixed.min.studytime,
        #   max(res$results$L[res.index],
        #       ceiling(min(df.sim.fix$total.time[which.max(cumsum(df.sim.fix$status)-ceiling(res$results$nEvents[res.index]))]))),
        #   ceiling(min(df.sim.fix$total.time[which.max(cumsum(df.sim.fix$status)-ceiling(res$results$nEvents[res.index]))])))
        res$results$events[res.index] = sum(BSSRed::cut_dfSim(df.sim.fix,res$results$study.end[res.index])$status)
        res$results$nPop[res.index] = NROW(df.sim.fix)
        # if (par.est){
        #   cat("parameter estimation not yet implemented")
        # } else {
        #   res$results$loglik[res.index] = NA
        #   res$results$waldtest[res.index] = NA
        #   res$results$sctest[res.index] = NA
        #   res$results$theta.est[res.index] = NA
        # }
      } # end loop j
    } # end loop i
  } # end loop l
  return(res)
} # end function

# a <- sim_BSSRed(recruitment=recruitment,addRecT = addRecT,
#                 lambdaC=lambdaC,sigma=sigma,theta.star = theta.star,gamma=gamma,kappa=kappa,
#                 distS = distS,distC = distC,nEvents=nEvents,L=L,
#                 nSim=nSim,par.est = par.est,fixed.min.studytime = F,seed = seed)
#
# a
# plot(NULL,xlim = c(-10,10),ylim = c(-10,10))
# text(0, 0, expression(paste("Enjoy the following equation!: ",hat(rho) == (sigma[t] * bar(lambda))[t[i]] +
#                               sum(x[i]^{-1}, i==1, n)+
#                               frac(x[i], n)^{-1})))
# ?text
