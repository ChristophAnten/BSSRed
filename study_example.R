# install globally:
library("devtools")
install_git(url="https://github.com/ChristophAnten/BSSRed.git",
            upgrade=FALSE)
library("BSSRed")


get_hazardRate <- function(prop,time,theta=NULL){
  if (is.null(theta)){
    return(-log(1-prop) / time)
  }
  if (is.numeric(theta)){
    return(-log(1-prop) /time *  theta)
  }
}
get_eventProportion <- function(hazard,time,theta=NULL){
  if (is.null(theta)){
    return(1-exp(-hazard * time))
  }
  if (is.numeric(theta)){
    return(1-exp(-hazard * time * theta))
  }
}

### The idea of an approach.

## simple approach (all crucial information)
# We have:
# - recruitment time 10 month
R <- 10
# - administrative censoring after 20 month (same time relation as R and t)
L <- 20
# - planned recruitments as a matrix
N <- matrix(round(c(dnorm(seq(-2,2,length.out = 20))*200,dnorm(seq(-2,2,length.out = 20))*400)),ncol=2)[1:R,]
# N <- cbind(c(seq(3,30,by=3),rep(34,5),rep(35,5)),
#            c(seq(6,60,by=6),rep(68,5),rep(70,5)))/2;R=NROW(N);L=39 # example from friedePaper
N;colSums(N)
# - regular recruitment (each month) with all recruitments at start of the month
tn <- 1:R - 1
# - assumed continuation of recruitments
# M <- matrix(round(c(dnorm(seq(-2,2,length.out = 20))*200,dnorm(seq(-2,2,length.out = 20))*400)),ncol=2)[(R+1):20,]
# - regular continuation of recruitment (each month)
# tm <- max(tn) + 1:NROW(M)
# - assumed reduction of hazard rate in placebo group by 30%
theta <- 1-.3
# - assuming 30% of patients have an event after 24 months
propP <- .3
timeP <- 24
lambdaP <- get_hazardRate(propP,timeP);lambdaP #BSSRed::
lambdaT <- get_hazardRate(propP,timeP,theta);lambdaT #BSSRed::
propT <- get_eventProportion(lambdaP,timeP,theta);propT #BSSRed::
# - a discontinuation rate of 20% after two years (24 months)
propC <- .2
timeC <- 24
gamma <- get_hazardRate(propC,timeC);gamma # BSSRed::
# - benötigte Ereigniszahl berechnen
nEvents <- BSSRed::nschoenfeld(theta = theta,k=2,alpha = .05,beta = .1, alternative = "two-sided");nEvents
# - expected Events with this setUp
BSSRed::calc_expEvents(N,theta,L,tn,lambdaP,gamma = gamma)
# - time needed with theese assumptions
BSSRed::rEstimate(N=N, tn=tn, theta=theta, nEvents=nEvents, L=L,
                  lambda=lambdaP, gamma = gamma, distC="exponential",distS="exponential")
# - adjusting the recruitment numbers
BSSRed::rEstimate(N=2*N, tn=tn, theta=theta, nEvents=nEvents, L=L,
                  lambda=lambdaP, gamma = gamma, distC="exponential",distS="exponential")# - building a single sample dataset
df.sim <- BSSRed::simSurvData(N=2*N, tn=tn, lambda=lambdaP, theta=theta,
                              gamma=gamma, distS = "exponential",distC = "exponential",L=L);head(df.sim)

plot_dfSim(df.sim)
abline(v=c(max(tn),L),lty=2,col=c("darkgrey","darkred"))
barplot(names.arg =tn,height=(t(cbind(cumsum(N[,1]),cumsum(N[,2])))),
        main = "Cummulative number of recruitments")



# - simulating study
res <- sim_BSSRed(N=2*N,tn=tn,M=2*N[R,,drop=F],tm=tn[R]+1,addRecT=c(0,3,6), # add number of patients or time
                  sim.lambdaP = get_hazardRate(seq(propP-.1,propP+.1,length.out = 9),timeP), sigma=1,
                  theta=theta,gamma=gamma,kappa=1,distS="exponential",distC="exponential",nEvents=nEvents,
                  L=L,nSim=100, par.est=FALSE, fixed.min.studytime=TRUE, seed=235711)
# res <- sim.BSSRed(N=par.list$N,tn=par.list$tn,M=par.list$M,tm=par.list$tm,addRecT=c(0,3,6), # add number of patients or time
#                   sim.lambdaP = par.list$sim.lambdaP, sigma=1,
#                   theta=theta,gamma=gamma,kappa=1,distS="exponential",distC="exponential",nEvents=nEvents,
#                   L=par.list$L,nSim=100, par.est=FALSE, fixed.min.studytime=TRUE, seed=235711)
res <- resWoPlus
res <- resWoPlusOnL
library(plyr)
library(dplyr)
library(ggplot2)
res$results %>%
  dplyr::select(-distS,-distC) %>%
  group_by(lambdaC,addRecT) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  mutate(lambdaC = 1-exp(-lambdaC*24)) %>%
  ggplot(aes(x=lambdaC,y=study.end,col=factor(addRecT))) +
  geom_point() +
  geom_hline(yintercept = L,lty=2,col="darkgrey") +
  geom_line() +
  ggtitle(sprintf("Mean Time to End of Study (nEvent=%i)",ceiling(nEvents))) +
  theme_bw() +
  ylab("Mean study duration") +
  xlab("24 month conifrmed progression probability") +
  scale_x_continuous(breaks=1-exp(-lambdaC*24))# +
res$results %>%
  dplyr::select(-distS,-distC) %>%
  group_by(lambdaC,addRecT) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  mutate(lambdaC = 1-exp(-lambdaC*24)) %>%
  ggplot(aes(x=lambdaC,y=nPop-sum(recruitment[,1:2]),col=factor(addRecT))) +
  geom_point() +
  geom_line() +
  ggtitle(sprintf("Mean number of Additional Patients (nEvent=%i)",ceiling(nEvents))) +
  ylab("Mean number of additional patients") +
  xlab("24 month conifrmed progression probability") +
  theme_bw() +
  scale_x_continuous(breaks=1-exp(-lambdaC*24))# +
# scale_y_continuous(breaks = seq(0,700,by=100),limits = c(0,700))
res$results %>%
  dplyr::select(-distS,-distC) %>%
  group_by(lambdaC,addRecT) %>%
  summarise_all(mean) %>%
  mutate(power = BSSRed::pschoenfeld(theta = .7,#theta.star,
                                     nEvent = events,
                                     k = 2,
                                     alpha = .05,
                                     alternative = "two-sided")$power) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  mutate(lambdaC = 1-exp(-lambdaC*24)) %>%
  ggplot(aes(x=lambdaC,y=power,col=factor(addRecT))) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1-.1,lty=2,col="darkgrey") +
  ggtitle(sprintf("Power at End of Study (nEvent=%i)",ceiling(nEvents))) +
  ylab("Mean power") +
  xlab("24 month conifrmed progression probability") +
  theme_bw()
