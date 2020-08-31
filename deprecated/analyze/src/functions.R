#!/usr/bin/env Rscript --vanilla
# set expandtab ts=4 sw=4 ai fileencoding=utf-8
#
# Author: KL
# Maintainer(s): PB, MG, JJ
# License: (c) HRDAG 2020, GPL v2 or newer
#
# Covid19-undercount/analyze/src/functions.R
# -----------------------------------------------------------
#

if (!require('pacman')) {install.packages('pacman')}
pacman::p_load(mvtnorm, EpiDynamics, extraDistr)

simulate.sir.synthetic.data <- function(N, beta, gammar, I, T, max.time, theta, p) {
  if (length(beta) == 1) {
    beta <- rep(beta, T - 1)
  }
  
  S <- N - I
  R <- 0
  Nu <- I

  ## add in deaths later
  Nus <- c(Nu, rep(0, T - 1)) #newly infected at time t
  Ss <- c(S, rep(0, T - 1)) #susceptible at time t
  Is <- c(I, rep(0, T - 1)) #infectious at time t
  Rs <- c(R, rep(0, T - 1)) #removed at time t = dead at time t + recovered at time t
  Ds <- rep(0, T + max.time + 1) #newly dead at time t
  save.death.times <- NULL
  
  for(t in 2:T){
    nu <- Ss[t - 1] * Is[t - 1] * beta[t - 1] / N
    
    #use Poisson process
    future.deaths <- rpois(1, p * nu)
    #print(future.deaths)

    death.times <- sample(0:max.time, future.deaths, prob=theta, replace=TRUE)
    
    save.death.times <- c(death.times, save.death.times)
    death.time.counts <- table(c(death.times, 0:max.time)) - 1
    Ds[t:(t+max.time)] <- Ds[t:(t + max.time)] + death.time.counts
    
    Rs[t] <- Rs[t - 1] + gammar * (Is[t - 1] - Ds[t - 1])  + Ds[t - 1]
    Ss[t] <- Ss[t - 1] - nu
    Is[t] <- Is[t - 1] + nu - gammar * (Is[t - 1] - Ds[t - 1]) - Ds[t - 1]
    Nus[t] <- nu
  }
  return(list(Nus=Nus, Ds=Ds, Ss=Ss, Is=Is, Rs=Rs, save.death.times=save.death.times))
}


make.future.deaths <- function(size, n, prob){
  death.counts <- sample.int(n=n, size=size, prob=prob, replace=TRUE)
  death.counts <- tabulate(death.counts, nbins=max.time)
  return(death.counts)
}


sample.Ds <- function(Nus1, theta, p){
  #sample the deaths given the time series of new infections
  #this could be done better with a weighted moving average
  T <- length(Nus1)
  max.time <- length(theta)

  new.deaths <- rpois(T, Nus1 * p) #number of future deaths generated from new infections at each time point
  new.deaths.times <- sapply(new.deaths, make.future.deaths, n=max.time, prob=theta) #number of days in the future for each of the deaths

  #re-align into the right locations
  Ds <- rep(0, T + max.time)
  for (i in 1:T){
    Ds[1:max.time + i] <- Ds[1:max.time + i] + new.deaths.times[,i]
  }

  #hack off the tail because the entries on the end have not yet accumulated all of the deaths there;
  #i.e. in the last position, those are only the number of deaths that were generated from new infections
  #on that day. We'd need to go forward max.time more days to see how many deaths would have accumulated in that day.

  Ds <- Ds[1:T]

  return(Ds)
}


calculate.death.means <- function(Nus1, theta, p){
  max.time <- length(theta)
  T <- length(Nus1)
  #weighted moving average of Nus, weighted by probability of an infection k days ago resulting in a death day
  mus <- stats::filter(c(rep(0, max.time - 1), Nus1 * p), filter=theta, method="convolution", sides=1)
  mus <- mus[max.time:(T + max.time - 1)]

  return(mus)
}


sample.Ds.ma <- function(Nus1, theta, p){
  #sample the deaths given the time series of new infections using a moving average
  max.time <- length(theta)
  T <- length(Nus1)
  mus <- calculate.death.means(Nus1, theta, p)

  Ds <- rpois(T, mus)

  return(Ds)
}


log.prior <- function(x, xi) {
  #the log of the uniform density (used for gammar | beta)
  #x : value at which to compute log of density
  #xi : vector of (a,b) for a Uniform(a,b) distribution

  return(dtnorm(x, mean=xi["mean"], sd=xi["sd"], a=xi["a"], b=xi["b"], log=T))
  
}


log.gprior <- function(x, zeta) {
  #the log of the uniform distribution (used for gammar)
  #x : value at which to compute log of density
  #zeta : vector of (a,b) for a Uniform(a,b) distribution
  #note : redundant with log.prior but want to leave these separate to easily allow possibility to change priors on gamma and beta in the future

  return(dtnorm(x, mean=zeta["mean"], sd=zeta["sd"], a=zeta["a"], b=zeta["b"], log=T))

}


log.T0prior <- function(x, zeta) {
  #the log of the uniform distribution (used for T0)
  #x : value at which to compute log of density
  #zeta : vector of (a,b) for a Uniform(a,b) distribution
  #note : redundant with log.prior but want to leave these separate to easily allow possibility to change priors on gamma and beta in the future

  if (x >= zeta[1] & x <= zeta[2]) {
    return(0)
  } else {
    return(-Inf)
  }
}


simulate.sir <- function(beta, gammar, I0, N, T0, T, S0=NA) {
  if (is.na(S0)) {S0 <- 1 - I0}
  
  integer.start <- T0 == ceiling(T0)
  
  tms <- seq(from=ceiling(T0), to=T, by=1)
  
  if(! integer.start){
    tms <- c(T0, tms)
  }
  
  sim1 <- SIR(pars=list(beta=beta, gamma=gammar),
              init=c(S=S0, I=I0, R=1 - S0 - I0),
              time=tms)
  
  SI.sim <- as.matrix(sim1$results[,2:3])
  if(!integer.start){
    SI.sim <- SI.sim[-1,]
  }
  
  SI <- t(matrix(c(1, 0), 2, T))
  SI[ceiling(T0):T,] <- SI.sim
  
  S <- SI[,1] * N
  Nus.sim <- S[1:(length(S) - 1)] - S[2:length(S)]
  Nus <- rep(0, T)
  Nus[2:T] <- Nus.sim
  
  if(!integer.start){
    tms <- tms[-1]
  }
  
  return(list(SI=SI, Nus=Nus, times=tms))
}


# eliminate max.time, theta

propose.sir <- function(beta, gammar, I, T, p, Nus, Ds, xi, zeta, ll.curr, SI,
                        k.adapt=0, c0=1, BG=NA, k=NA, mn.bg=NA, SXXt=NA, c1=NA,
                        T0=1, zeta.T0) {
# propose new parameters of SIR model
  #beta : beta from last iteration
  #gammar : gammar from last iteration
  #I : intial number of infections
  #T : length of time series
  #theta.t : an expanded version of theta that makes computation easier
  #p : IFR
  #Nus : Nus at last iteration
  #Ds : time series of deaths
  #xi : parameters of prior on beta | gammar
  #zeta : parameters of prior on gamma
  #k.adapt : when do we want to start adapting
  #c0 : before adaptation, metropolis proposal variance
  #BG : matrix of mcmc samples
  #k :  current iteration number
  #mn.bg : running mean of parameters
  #SXXt :  sum_{t=1}^k x_t x_t' where x_t is the t'th sample of the parameters
  #c1 : scaling factor for covariance of metropolis proposals after adaptation begins
  #T0 : day of first infection
  #zeta.T0 : parameters of prior on T0

     if (k < k.adapt) {
       prop.bg <- rnorm(3, c(beta, gammar, T0), sqrt(c0))
     }
    if (k == k.adapt) {
      mn.bg <- matrix(apply(BG[1:(k - 1),], 2, mean), 3, 1)

      SXXt = t(BG[1:(k - 1),]) %*% BG[1:(k - 1),]
      S.bg <- 1 / (k - 2) * SXXt - (k - 1) / (k - 2) * mn.bg %*% t(mn.bg)
      S.bg <- c1 * S.bg #+ c*Eps*diag(d)

      prop.bg <- matrix(rmvnorm(1, c(beta, gammar, T0), S.bg), 3, 1)
    }
    if (k > k.adapt) {
      bg <- matrix(c(beta, gammar, T0), 3, 1)
      mn.bg <- matrix(((k - 2) * mn.bg + bg) / (k - 1), 3, 1)

      SXXt <- SXXt + bg %*% t(bg)
      S.bg <- 1 / (k - 2) * SXXt - (k - 1) / (k - 2) * mn.bg %*% t(mn.bg)
      S.bg <- c1 * S.bg #+ c*Eps*diag(d)

      prop.bg <- matrix(rmvnorm(1, bg, S.bg), 3, 1)
    }
    beta.prop <- prop.bg[1]
    gammar.prop <- prop.bg[2]
    T0.prop <- prop.bg[3]

  lp.prop <- (log.prior(beta.prop / gammar.prop, xi) 
              + log.gprior(1 / gammar.prop, zeta) 
              + log.T0prior(T0.prop, zeta.T0))

  if (!is.infinite(lp.prop)) {

    sim1 <- simulate.sir(beta.prop, gammar.prop, I, N, T0.prop, T)
    Nus.prop <- sim1$Nus
    SI.prop <- sim1$SI

    mus.prop <- calculate.death.means(Nus.prop, theta, p)
    ll.prop <- sum(dpois(Ds, mus.prop, log=T)) + lp.prop

    log.prob <- ll.prop - ll.curr
    acc <- exp(log.prob) > runif(1)
  } else {
    acc <- FALSE
  }

  if (is.na(acc) | is.infinite(acc)) {
    print('check')
  }

  if (acc) {
    beta <- beta.prop
    gammar <- gammar.prop
    T0 <- T0.prop
    Nus <- Nus.prop
    ll.curr <- ll.prop
    SI <- SI.prop
  }

  return(list(beta=beta, Nus=Nus, gammar=gammar, acc=acc, mn.bg=mn.bg, SXXt=SXXt,
              T0=T0, SI=SI, ll.curr=ll.curr))
}


run.mcmc <- function(nmc, theta, beta, gammar, T, I, Ds, Nus, xi, zeta, disp.int,
                     k.adapt=0, c0=1, c1=1, T0=1, zeta.T0) {

# the big function that runs mcmc t oestimate beta, gamma, and T0 from data
  #nmc : the number of iterations to run mcmc
  #theta : probability of dying on each day post infection, conditional on death
  #beta : initial value of beta
  #gammar : initial value of gammar
  #i.max : same
  #T: length of time series
  #I : number of initial infections
  #max.time : maximum days to death allowable
  #Ds: the data, the time series of deaths
  #Nus : initial value for the daily number of new infections
  #xi : parameters of prior on beta | gammar
  #zeta : parameters of prior on gammar
  #disp.int : how often to display trace plots
  #k.adapt : when to start adapting
  #c0: ?
  #c1: ?
  #T0 : time of first infection
  #zeta.T0 : parameters of prior on T0

  #initialize
  
  ACC <- matrix(0, nmc, 1)
  BETA <- matrix(0, nmc, 1)
  GAMMAR <- matrix(0, nmc, 1)
  BG <- matrix(0, nmc, 3)
  NUS <- matrix(0, nmc, length(Nus))
  T0S <- matrix(0, nmc, 1)
  Nus.True <- Nus
  mn.bg <- 0
  SXXt <- 0
  SIS <- array(0, dim=c(T, 2, nmc))

  sim1 <- simulate.sir(beta, gammar, I, N, T0, T)
  Nus <- sim1$Nus
  SI <- sim1$SI
  
  lp.curr <- (log.prior(beta / gammar, xi)
              + log.gprior(1 / gammar, zeta)
              + log.T0prior(T0, zeta.T0))
  mus <- calculate.death.means(Nus, theta, p)
  ll.curr <- sum(dpois(Ds, mus, log=T)) + lp.curr
  
  t1 <- proc.time()
  for (t in 1:nmc) {

    # sample parameters
    eta <- propose.sir(beta, gammar, I, T, p, Nus, Ds, xi, zeta, ll.curr, SI,
                       k.adapt=k.adapt, c0=c0, BG=BG, k=t, mn.bg=mn.bg,
                       SXXt=SXXt, c1=c1, T0=T0, zeta.T0=zeta.T0)
    mn.bg <- eta$mn.bg
    SXXt <- eta$SXXt
    SI <- eta$SI
    ll.curr <- eta$ll.curr

    # save parameters
    Nus <- eta$Nus
    ACC[t] <- eta$acc
    beta <- eta$beta
    gammar <- eta$gammar
    T0 <- eta$T0
    T0S[t] <- eta$T0
    BETA[t] <- beta
    GAMMAR[t] <- gammar
    NUS[t,] <- Nus
    SIS[,,t] <- SI

    BG[t,] <- c(beta, gammar, T0)

    if (t %% disp.int == 0) {
      t2 <- proc.time()
      print(paste('elapsed time', t2[3] - t1[3]))
      t1 <- proc.time()
      
      par(mfrow=c(3, 1))
      plot(GAMMAR[1:t], main=expression(gamma))
      plot(BETA[1:t], main=expression(beta))
      plot(T0S[1:t], main=expression(T[0]))

      print(paste('acceptance rate in the last', disp.int, 'iterations was',round(mean(ACC[(t - disp.int + 1):t]), 3)))
    }
  }

  return(list(BETA=BETA, NUS=NUS, ACC=ACC, GAMMAR=GAMMAR, T0S=T0S, SIS=SIS))
}


log.phiprior <- function(x, zeta) {
  if (x >= zeta[1] & x <= zeta[2]) {
    return(0)
  } else {
    return(-Inf)
  }
}


simulate.sir.timebreak <- function(beta, gammar, I0, N, T0, T, T1, phi){
  #beta : parameter of SIR
  #gammar : parameter of SIR
  #I0 : initial infections at time T0
  #N : population size
  #T : total time in time series
  #T1 : time at which social distancing starts
  #phi : percent of original contacts under social distancing
  
  integer.start <- T0 == ceiling(T0)
  
  #simulate epidemic up to intervention time
  sim1 <- simulate.sir(beta, gammar, I0, N, T0, T1 - 1)
  
  #give simulation at end of sim1 as initial conditions to sim 2, cut beta by phi
  sim2 <- simulate.sir(beta * phi, gammar, sim1$SI[T1 - 1,2], N, T1 - 1, T, sim1$SI[T1 - 1,1])
  
  #append Nus from both sims
  Nus <- rep(0, T)
  Nus[sim2$times] <- sim2$Nus[sim2$times] 
  Nus[sim1$times] <- sim1$Nus[sim1$times]
  
  #append SI from both sims, remove first element from sim 2 because it was 
  # initial conditions of sim 2 and is same as last element of sim1
  SI <- rbind(sim1$SI, sim2$SI[sim2$times[-1],])
  
  return(list(Nus=Nus, SI=SI))
}


# eliminate max.time, theta
propose.sir.state <- function(beta, gammar, I, T, p, Nus, Ds, xi, zeta, ll.curr, SI,
                        k.adapt=0, c0=1, X=NA, k=NA, X.bar=NA, SXXt=NA, c1=NA,
                        T0=1, zeta.T0, phi=NA, T1=NA, zeta.phi=NA) {
# propose new parameters of SIR model
  #beta : beta from last iteration
  #gammar : gammar from last iteration
  #I : intial number of infections
  #T : length of time series
  #theta.t : an expanded version of theta that makes computation easier
  #p : IFR
  #Nus : Nus at last iteration
  #Ds : time series of deaths
  #xi : parameters of prior on beta | gammar
  #zeta : parameters of prior on gamma
  #k.adapt : when do we want to start adapting
  #c0 : before adaptation, metropolis proposal variance
  #BG : matrix of mcmc samples
  #k :  current iteration number
  #mn.bg : running mean of parameters
  #SXXt :  sum_{t=1}^k x_t x_t' where x_t is the t'th sample of the parameters
  #c1 : scaling factor for covariance of metropolis proposals after adaptation begins
  #T0 : day of first infection
  #zeta.T0 : parameters of prior on T0

     d <- 4
     x <- c(beta, gammar, T0, phi)
  
     if (k < k.adapt) {
       prop.x <- rnorm(d, x, sqrt(c0))
     }
    if (k == k.adapt) {
      X.bar <- matrix(apply(X[1:(k - 1),], 2, mean), d, 1)

      SXXt = t(X[1:(k - 1),]) %*% X[1:(k - 1),]
      S.X <- 1 / (k - 2) * SXXt - (k - 1) / (k - 2) * X.bar %*% t(X.bar)
      S.X <- c1 * S.X #+ c*Eps*diag(d)

      prop.x <- matrix(rmvnorm(1, x, S.X), d, 1)
    }
    if (k > k.adapt) {
      X.bar <- matrix(((k - 2) * X.bar + x) / (k - 1), d, 1)

      SXXt <- SXXt + x %*% t(x)
      S.X <- 1 / (k-2) * SXXt - (k - 1) / (k - 2) * X.bar %*% t(X.bar)
      S.X <- c1 * S.X #+ c*Eps*diag(d)

      prop.x <- matrix(rmvnorm(1, x, S.X), d, 1)
    }
    beta.prop <- prop.x[1]
    gammar.prop <- prop.x[2]
    T0.prop <- prop.x[3]
    phi.prop <- prop.x[4]

  lp.prop <- (log.prior(beta.prop / gammar.prop, xi) 
              + log.gprior(1 / gammar.prop, zeta) 
              + log.T0prior(T0.prop, zeta.T0)
              + log.phiprior(phi.prop, zeta.phi))

  if (!is.infinite(lp.prop)) {

    sim1 <- simulate.sir.timebreak(beta.prop, gammar.prop, I, N, T0.prop, T, T1, phi.prop)
    
    Nus.prop <- sim1$Nus
    SI.prop <- sim1$SI

    mus.prop <- calculate.death.means(Nus.prop, theta, p)
    ll.prop <- sum(dpois(Ds, mus.prop, log=T)) + lp.prop

    log.prob <- ll.prop - ll.curr
    acc <- exp(log.prob) > runif(1)
  } else {
    acc <- FALSE
  }

  if (is.na(acc) | is.infinite(acc)) {
    print('check')
  }

  if (acc) {
    beta <- beta.prop
    gammar <- gammar.prop
    T0 <- T0.prop
    Nus <- Nus.prop
    phi <- phi.prop
    ll.curr <- ll.prop
    SI <- SI.prop
  }

  return(list(beta=beta, Nus=Nus, gammar=gammar, acc=acc, X.bar=X.bar, SXXt=SXXt,
              T0=T0, SI=SI, ll.curr=ll.curr, phi=phi))
}



run.mcmc.state <- function(nmc, theta, beta, gammar, T, I, Ds, Nus, xi, zeta, disp.int,
                     k.adapt=0, c0=1, c1=1, T0=1, zeta.T0, phi, T1, zeta.phi, plotting=T) {

# the big function that runs mcmc to estimate beta, gamma, and T0 from data
  #nmc : the number of iterations to run mcmc
  #theta : probability of dying on each day post infection, conditional on death
  #beta : initial value of beta
  #gammar : initial value of gammar
  #i.max : same
  #T: length of time series
  #I : number of initial infections
  #max.time : maximum days to death allowable
  #Ds: the data, the time series of deaths
  #Nus : initial value for the daily number of new infections
  #xi : parameters of prior on beta | gammar
  #zeta : parameters of prior on gammar
  #disp.int : how often to display trace plots
  #k.adapt : when to start adapting
  #c0: ?
  #c1: ?
  #T0 : time of first infection
  #zeta.T0 : parameters of prior on T0

  #initialize
  d <- 4
  
  ACC <- matrix(0, nmc, 1)
  BETA <- matrix(0, nmc, 1)
  GAMMAR <- matrix(0, nmc, 1)
  X <- matrix(0, nmc, d)
  NUS <- matrix(0, nmc, length(Nus))
  T0S <- matrix(0, nmc, 1)
  PHIS <- matrix(0, nmc, 1)
  Nus.True <- Nus
  X.bar <- 0
  SXXt <- 0
  SIS <- array(0, dim=c(T, 2, nmc))

  sim1 <- simulate.sir.timebreak(beta, gammar, I, N, T0, T, T1, phi)
  Nus <- sim1$Nus
  SI <- sim1$SI

  lp.curr <- (log.prior(beta / gammar, xi)
              + log.gprior(1 / gammar, zeta)
              + log.T0prior(T0, zeta.T0)
              + log.phiprior(phi, zeta.phi))
  mus <- calculate.death.means(Nus, theta, p)
  ll.curr <- sum(dpois(Ds, mus, log=T)) + lp.curr
  
  t1 <- proc.time()
  for (t in 1:nmc) {
    # sample parameters
    eta <- propose.sir.state(beta, gammar, I, T, p, Nus, Ds, xi, zeta, ll.curr,
                             SI, k.adapt=k.adapt, c0=c0, X=X, k=t, X.bar=X.bar,
                             SXXt=SXXt, c1=c1, T0=T0, zeta.T0=zeta.T0, phi, T1, zeta.phi)
    X.bar <- eta$X.bar
    SXXt <- eta$SXXt
    SI <- eta$SI
    ll.curr <- eta$ll.curr
    phi <- eta$phi

    # save parameters
    Nus <- eta$Nus
    ACC[t] <- eta$acc
    beta <- eta$beta
    gammar <- eta$gammar
    T0 <- eta$T0
    T0S[t] <- eta$T0
    BETA[t] <- beta
    GAMMAR[t] <- gammar
    NUS[t,] <- Nus
    SIS[,,t] <- SI
    PHIS[t] <- phi

    if(t > k.adapt){
      last.acc <- mean(ACC[t:(t - 99)])
      if(last.acc < 0.001){
        print("check")
      }
    }
    X[t,] <- c(beta, gammar, T0, phi)

    if (t %% disp.int == 0) {
      t2 <- proc.time()
      print(paste('elapsed time', t2[3] - t1[3]))
      t1 <- proc.time()
      
      if (plotting) {
        par(mfrow=c(4, 1))
        plot(GAMMAR[1:t], main=expression(gamma))
        plot(BETA[1:t], main=expression(beta))
        plot(T0S[1:t], main=expression(T[0]))
        plot(PHIS[1:t], main=expression(phi))
      }

      print(paste('acceptance rate in the last', disp.int, 'iterations was', round(mean(ACC[(t - disp.int + 1):t]), 3)))
    }
  }

  return(list(BETA=BETA, NUS=NUS, ACC=ACC, GAMMAR=GAMMAR, T0S=T0S, SIS=SIS, PHIS=PHIS))
}

# done.
