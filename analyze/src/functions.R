if (!require('pacman')) {install.packages('pacman')}
pacman::p_load(mvtnorm)

#sample Xs as in original manuscript, we actually don't do this anymore, but I think we should keep it in case for some reason we decide we want to use it in the future
sample.Xs <- function(Ds,Nus,theta,t.max,i.max) {

  X <- matrix(0,t.max,t.max)


  for (t in 1:t.max) {
      idx <- (max(t-i.max + 1,1)):t
      nu.t <- Nus[idx]
      if (any(nu.t>0)) {
        wts <-rev(theta[1:length(idx)])*nu.t
        X.t <- rmultinom(1,Ds[t],wts)
        X[idx,t] <- X.t
      } else {
        X[idx,t] <- 0
      }
  }

return(X)
}


#simulate the SIR models with parameters given
simulate.sir <- function(N,beta,gammar,I,T,max.time,theta,p,T0,Ds=NA) {

  #N : population size
  #beta : beta in standard sir model
  #gammar : rate of recovery in sir model
  #I : number of initial infections at time T0
  #T : total length of simulation
  #max.time : maximum time we consider for time from infection to death
  #theta : vector giving probability of death at each date post infection for those who do die
  #p : probability of death given infection
  #T0 : date of first infection
  #Ds : optional vector of death counts. If not included, Ds are simulated directly from model. If known, they are incorporated into simulated trajectory.


  if (length(beta)==1) {
    beta <- rep(beta,T-1)
  }


  if (!is.na(Ds[1])) {
    input.Ds<-TRUE
  } else {
    input.Ds <- FALSE
  }

  #initialize
  S <- N-I
  R <- 0
  Nu <- I

  Nus <- c(rep(0,T0-1),Nu, rep(0, T-T0-1)) #newly infected at time t
  Ss <- c(rep(0,T0-1),S, rep(0, T-T0-1)) #susceptible at time t
  Is <- c(rep(0,T0-1),I, rep(0, T-T0-1)) #infectious at time t
  Rs <- c(rep(0,T0-1),R, rep(0, T-T0-1)) #removed at time t = dead at time t + recovered at time t
  if (!input.Ds) {
    Ds <- rep(0, T + max.time + 1) #newly dead at time t
  }
  save.death.times <- NULL

  #simulate modified SIR model to account for known deaths
  for(t in (T0+1):T){
    nu <- Ss[t-1]*Is[t-1]*beta[t-1]/N

    #if deaths unknown simulate them from model
    future.deaths <- rpois(1, p*nu)
    death.times <- sample(0:max.time, future.deaths, prob=theta, replace=TRUE)


    save.death.times <- c(death.times, save.death.times)
    #find out how many deaths were generated at time t that are going to occur on each day in the future
    death.time.counts <- table(c(death.times, 0:max.time)) - 1
    if (!(input.Ds)) {
      Ds[t:(t+max.time)] <- Ds[t:(t + max.time)] + death.time.counts
    }

    #Rs[t] <- Rs[t-1] + gammar*Is[t-1]

    #update SIR model
    Rs[t] <- Rs[t-1] + gammar*(Is[t-1] - Ds[t-1])  + Ds[t-1]
    Ss[t] <- Ss[t-1] - nu
    Is[t] <- Is[t-1] + nu - gammar*(Is[t-1] - Ds[t-1]) - Ds[t-1]
    Nus[t] <- nu
  }
  return(list(Nus=Nus,Ds=Ds,Ss=Ss,Is=Is,Rs=Rs,save.death.times=save.death.times))
}

log.prior <- function(x,xi) {
  #the log of the uniform density (used for gammar | beta)
  #x : value at which to compute log of density
  #xi : vector of (a,b) for a Uniform(a,b) distribution

  if (x>xi[1] & x<xi[2]) {
    return(0)
  } else {
    return(-Inf)
  }
}

log.gprior <- function(x,zeta) {
  #the log of the uniform distribution (used for gammar)
  #x : value at which to compute log of density
  #zeta : vector of (a,b) for a Uniform(a,b) distribution
  #note : redundant with log.prior but want to leave these separate to easily allow possibility to change priors on gamma and beta in the future

  if (x>zeta[1] & x<zeta[2]) {
    return(0)
  } else {
    return(-Inf)
  }
}

log.T0prior <- function(x,zeta) {
  #the log of the uniform distribution (used for T0)
  #x : value at which to compute log of density
  #zeta : vector of (a,b) for a Uniform(a,b) distribution
  #note : redundant with log.prior but want to leave these separate to easily allow possibility to change priors on gamma and beta in the future

  if (x>=zeta[1] & x<=zeta[2]) {
    return(0)
  } else {
    return(-Inf)
  }
}


propose.sir <- function(beta,gammar,I,T,max.time,theta,theta.t,p,Nus,Ds,xi,zeta,
                        k.adapt=0,c0=1,BG=NA,k=NA,mn.bg=NA,SXXt=NA,c1=NA,T0=1,zeta.T0) {
# propose new parameters of SIR model
  #beta : beta from last iteration
  #gammar : gammar from last iteration
  #I : intial number of infections
  #T : length of time series
  #max.time : maximum possible time for death to occur after infection
  #theta : probability of death on each day post-infection
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


     if (k<k.adapt) {
       prop.bg <- rnorm(3,c(beta,gammar,T0),sqrt(c0))
     }
    if (k==k.adapt) {
      mn.bg <- matrix(apply(BG[1:(k-1),],2,mean),3,1)

      SXXt = t(BG[1:(k-1),])%*%BG[1:(k-1),]
      S.bg <- 1/(k-2)*SXXt-(k-1)/(k-2)*mn.bg%*%t(mn.bg)
      S.bg <- c1*S.bg #+ c*Eps*diag(d)

      prop.bg <- matrix(rmvnorm(1,c(beta,gammar,T0),S.bg),3,1)
    }
    if (k>k.adapt) {
      bg <- matrix(c(beta,gammar,T0),3,1)
      mn.bg <- matrix(((k-2)*mn.bg + bg)/(k-1),3,1)

      SXXt <- SXXt + bg%*%t(bg)
      S.bg <- 1/(k-2)*SXXt-(k-1)/(k-2)*mn.bg%*%t(mn.bg)
      S.bg <- c1*S.bg #+ c*Eps*diag(d)

      prop.bg <- matrix(rmvnorm(1,bg,S.bg),3,1)
    }
    log.com.add <- 0
    beta.prop <- prop.bg[1]
    gammar.prop <- prop.bg[2]
    T0.prop <- prop.bg[3]


  lp.prop <- log.prior(beta.prop/gammar.prop,xi) + log.gprior(1/gammar.prop,zeta) + log.T0prior(T0.prop,zeta.T0)

  if (!is.infinite(lp.prop)) {
    sim1 <- simulate.sir(N,beta.prop,gammar.prop,I,T,max.time,theta,p,floor(T0.prop),Ds)
    sim2 <- simulate.sir(N,beta.prop,gammar.prop,I,T,max.time,theta,p,ceiling(T0.prop),Ds)
    Nus.prop <- sim1$Nus*(T0.prop-floor(T0.prop)) + sim2$Nus*(ceiling(T0.prop)-T0.prop)

    ll.prop <- sum(dpois(Ds,apply(Nus.prop*p*theta.t,2,sum),log=T)) + lp.prop
    ll.curr <- sum(dpois(Ds,apply(Nus*p*theta.t,2,sum),log=T)) + log.prior(beta/gammar,xi) + log.gprior(1/gammar,zeta) +
      log.T0prior(T0,zeta.T0)
    log.prob <- ll.prop-ll.curr + log.com.add
    acc <- exp(log.prob)>runif(1)
  } else {
    acc<-FALSE
  }

  if (is.na(acc) | is.infinite(acc)) {
    print('check')
  }

  if (acc) {
    beta <- beta.prop
    gammar <- gammar.prop
    T0 <- T0.prop
    Nus <- Nus.prop
  }
  return(list(beta=beta,Nus=Nus,gammar=gammar,acc=acc,mn.bg=mn.bg,SXXt=SXXt,T0=T0))

}


make.theta.factor <- function(theta,T,i.max) {
  #function to expand theta into a matrix for easier book-keeping
    theta.tmp <- t(matrix(theta, ncol=T, nrow = i.max))
  j <- 2
  for(i in T:(T-i.max+2)){
    theta.tmp[i,j:i.max] <- 0
    j <- j + 1
  }
theta.factor <- apply(theta.tmp, 1, sum)
return(theta.factor)
}




run.mcmc <- function(nmc,theta,beta,gammar,T,I,max.time,Ds,Nus,xi,zeta,disp.int,
                     k.adapt=0,c0=1,c1=1,T0=1,zeta.T0) {

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
  theta.t <- matrix(0,T,T)
  for (j in 1:T) {
    if ((j+length(theta))<=T) {
      theta.t[j,j:(j+length(theta)-1)] <- theta
    } else {
      theta.t[j,j:T] <- theta[1:length(seq(from=j,to=T))]
    }
  }

  ACC <- matrix(0,nmc,1)
  BETA <- matrix(0,nmc,1)
  GAMMAR <- matrix(0,nmc,1)
  BG <- matrix(0,nmc,3)
  NUS <- matrix(0,nmc,length(Nus))
  T0S <- matrix(0,nmc,1)
  Nus.True <- Nus
  mn.bg <- 0
  SXXt <- 0

  for (t in 1:nmc) {
    #print(t)
    # sample X
    #X <- sample.Xs(Ds,Nus,theta,t.max,i.max)
    #rhos <- apply(X,1,sum)

    # sample parameters
    eta <- propose.sir(beta,gammar,I,T,max.time,theta,theta.t,p,Nus,Ds,xi,zeta,
                       k.adapt=k.adapt,c0=c0,BG=BG,k=t, mn.bg=mn.bg,SXXt =SXXt,c1=c1,T0=T0,zeta.T0=zeta.T0)
    mn.bg <- eta$mn.bg
    SXXt <- eta$SXXt



    # save parameters
    Nus <- eta$Nus
    #rhos <- apply(X,1,sum)
    ACC[t] <- eta$acc
    beta <- eta$beta
    gammar <- eta$gammar
    T0 <- eta$T0
    T0S[t] <- eta$T0
    BETA[t] <- beta
    GAMMAR[t] <- gammar
    NUS[t,] <- Nus

    BG[t,] <- c(beta,gammar,T0)


    if (t%%disp.int==0) {
      par(mfrow=c(3,1))
      plot(GAMMAR[1:t],main=expression(gamma))
      plot(BETA[1:t],main=expression(beta))
      plot(T0S[1:t],main=expression(T[0]))

      print(paste('acceptance rate in the last',disp.int,'iterations was',round(mean(ACC[(t-disp.int+1):t]),3)))
    }
  }

  return(list(BETA=BETA,NUS=NUS,ACC=ACC,GAMMAR=GAMMAR,T0S=T0S))


}


simulate.sir.basic <- function(N,beta,gamma,I,T) {
  #SIR without separately accounting for deaths

  S <- N-I
  R <- 0
  Nu <- I

  Nus <- c(Nu, rep(0, T-1)) #newly infected at time t
  Ss <- c(S, rep(0, T-1)) #susceptible at time t
  Is <- c(I, rep(0, T-1)) #infectious at time t
  Rs <- c(R, rep(0, T-1)) #removed at time t = dead at time t + recovered at time t

  for(t in 2:T){
    nu <- Ss[t-1]*Is[t-1]*beta/N

    Rs[t] <- Rs[t-1] + gamma*Is[t-1]
    Ss[t] <- Ss[t-1] - nu
    Is[t] <- Is[t-1] + nu - gamma*Is[t-1]
    Nus[t] <- nu
  }
  return(list(Nus=Nus,Ss=Ss,Is=Is,Rs=Rs))
}


