#
# Authors:     KL
# Maintainers: JJ, PB, MG
# Copyright:  
# =========================================
# Covid19-modeling/analyze/src/2-fit_model_to_real_data.R

#this code fits a model to the data

if (!require('pacman')) {install.packages('pacman')}
pacman::p_load(ggplot2, dplyr, tidyr, readr, lubridate, 
               future, future.apply, here, glue, argparse)

parser <- ArgumentParser()
parser$add_argument('--functions', default=here::here('analyze/src/functions.R'))
parser$add_argument('--theta', default=here::here('analyze/output/theta.Rds'))
parser$add_argument('--pops', default=here::here('analyze/input/country_pops.RData'))
parser$add_argument('--worldometers', default=here::here('analyze/input/worldounts_deaths_us.csv'))
parser$add_argument('--p', default=0.01)
parsre$add_argument('--country', default='us')

args <- parser$parse_args()

args[['mcmc']] <- here::here(glue('analyze/output/mcmc_run_{p}.RData'))
args[['bigDs']] <- here::here(glue('analyze/output/bigDs_{p}.RData'))
args[['bigNus']] <- here::here(glue('analyze/output/bigNus_{p}.RData'))

source(args$functions)

##----- set up
p <- 0.01 # the IFR. This is an input/assumption of the model.
R0.case <- seq(from=0.75, to=2.25, by=0.25) #different R0 cases to consider
n.cases <- length(R0.case)

#read in values for theta (time to death distribution as calculated in previous setup script)
theta <- readRDS(args$theta)
#read in population size data
pops <- readRDS(args$pops)


#maximum days until death (minus one because 0 days is allowed)
max.time <- length(theta) - 1
N <- pops[rownames(pops) == 'us']

#how long since mitigation went into place
days.since.mitigation <- 10
#how many days into the future do we want to predict
n.days.future <- 21

nmc <- 50000
burn <- 10000
n.thin <- 1000
thinned <- round(seq(burn, nmc, length.out=n.thin))
disp.int <- 1000

# worldometers data entered manually from website
dat <- read_csv(args$worldometers)
dat <- dat %>% mutate(date=mdy(date))
Ds <- dat$deaths_us
Ds <- c(rep(0, 14), Ds) #add on some extra 0's to the front so that if the model wants to it can estimate an earlier first infection date

# In a "normal" SIR model, the R0 is beta/gamma
# The average days to recovery is 1/gamma
# Let's choose
gammar <- 1 / 6
zeta <- c(4, 25) # limits of uniform prior on gammar^(-1)
beta <- gammar * 2.5
xi <- c(1.4, 3.9) # prior on beta | gammar  is Uniform(xi[1]*gammar,xi[2]*gammar) on this interval

#length of data time series
T <- length(Ds)

#number of initial infections
I <- 1

#time of first infection
T0 <- 15


##----- run mcmc
#find a reasonable place to initialize parameters
sim1 <- simulate.sir(N,beta,gammar,I,T+max.time + 1,max.time,theta,p,T0)
sim1$Ds
plot(sim1$Ds[1:T])
plot(sim1$Ds[1:T], Ds, xlim=c(0,max(sim1$Ds[1:T], Ds)), ylim=c(0,max(sim1$Ds[1:T], Ds)))
abline(0,1)


mcmc.out <- run.mcmc(nmc=nmc, theta=theta, beta=beta, gammar=gammar, T=T, I=I,
                     max.time=max.time, Ds=Ds, Nus=sim1$Nus[1:T], xi=xi,
                     zeta=zeta, disp.int=disp.int, k.adapt=5000,
                     c0=c(0.01^2, 0.01^2, 0.75^2), c1=0.5 * (2.38 / 3), T0=T0,
                     zeta.T0=c(1, 30))

##save output here from running mcmc
save(file=args$mcmc, mcmc.out)

##post-process MCMC output to be more easily plotted.... in next script (this takes a while to run, like the mcmc)
NUS <- mcmc.out$NUS
tmp <- simulate.sir(N,
                    mean(mcmc.out$BETA), 
                    mean(mcmc.out$GAMMAR),
                    I, 
                    T,
                    max.time, 
                    theta, 
                    p, 
                    floor(mean(mcmc.out$T0S)))
plot(tmp$Ds[1:T], Ds)
abline(0, 1)

#prepare data for plots
n.days.future <- 21
n.thin <- 1000
thinned <- round(seq(burn, nmc, length.out=n.thin)) # here burn should be greater than k.adapt from line 91
big.Ds <-  array(0, dim=c(T + n.days.future, length(thinned), n.cases + 1))
big.Nus <- array(0, dim=c(T + n.days.future, length(thinned), n.cases+1))

BETA.thinned <- mcmc.out$BETA[thinned]
GAMMAR.thinned <- mcmc.out$GAMMAR[thinned]
T0.thinned <- mcmc.out$T0S[thinned]

plan(multisession)
BGt <- split(cbind(BETA.thinned, GAMMAR.thinned, T0.thinned), seq(1:n.thin))

#for each mcmc sample simulate sir model under those model parameters; times are integers so use T0 as floor of estimated  T0,
tmp1 <- future_lapply(BGt, function(x){simulate.sir(N, 
                                                    x[1], 
                                                    x[2], 
                                                    I, 
                                                    T + n.days.future, 
                                                    max.time,
                                                    theta, 
                                                    p, 
                                                    floor(x[3]))})

#for each mcmc sample simulate sir model under those model parameters; times are integers so use T0 as ceiling of estimated  T0,
tmp2 <- future_lapply(BGt, function(x){simulate.sir(N, 
                                                    x[1], 
                                                    x[2], 
                                                    I,
                                                    T + n.days.future, 
                                                    max.time,
                                                    theta, 
                                                    p, 
                                                    ceiling(x[3]))})

#extract the samples of number of deaths under this model
tmp01 <- sapply(tmp1, function(x){return(x$Ds[1:(T + n.days.future)])})
tmp02 <- sapply(tmp2, function(x){return(x$Ds[1:(T + n.days.future)])})
#save the samples of D we just made (in lines 124-125)
big.Ds[,,n.cases+1] <- (tmp01 * (T0.thinned - floor(T0.thinned))
                        + tmp02 * (ceiling(T0.thinned) - T0.thinned)) #put simulation with actual parameters in last position

#extract samples Nu from each sample run
tmp01 <- sapply(tmp1, function(x){return(x$Nu[1:(T + n.days.future)])})
tmp02 <- sapply(tmp2, function(x){return(x$Nu[1:(T + n.days.future)])})
#save the extracted Nus
big.Nus[,,n.cases+1] <- (tmp01 * (T0.thinned - floor(T0.thinned)) 
                         + tmp02 * (ceiling(T0.thinned) - T0.thinned)) #put simulation with actual parameters in last position

days.since.mitigation <- 10

for (k in 1:n.cases) {
  print(k)
  tmp1 <- future_lapply(BGt, function(x){simulate.sir(N, 
                                                      c(rep(x[1], T - days.since.mitigation),
                                                        rep(x[2], days.since.mitigation + n.days.future) * R0.case[k]), #sample SIR as  though the R0 decreased by R0.case[k] days.since.mitigation ago
                                                      x[2], 
                                                      I, 
                                                      T + n.days.future, 
                                                      max.time, 
                                                      theta, 
                                                      p, 
                                                      floor(x[3]))})
  tmp2 <- future_lapply(BGt, function(x){simulate.sir(N,
                                                      c(rep(x[1], T - days.since.mitigation),
                                                        rep(x[2], days.since.mitigation + n.days.future) * R0.case[k]),
                                                      x[2], 
                                                      I, 
                                                      T + n.days.future, 
                                                      max.time, 
                                                      theta, 
                                                      p, 
                                                      ceiling(x[3]))})
  
  #extract samples of deaths
  tmp01 <- sapply(tmp1, function(x){return(x$Ds[1:(T + n.days.future)])})
  tmp02 <- sapply(tmp2, function(x){return(x$Ds[1:(T + n.days.future)])})
  big.Ds[,,k] <- (tmp01 * (T0.thinned - floor(T0.thinned)) 
                  + tmp02 * (ceiling(T0.thinned) - T0.thinned))
  
  #extract Nus
  tmp01 <- sapply(tmp1, function(x){return(x$Nus)})
  tmp02 <- sapply(tmp2, function(x){return(x$Nus)})
  big.Nus[,,k] <- (tmp01 * (T0.thinned - floor(T0.thinned))
                   + tmp02 * (ceiling(T0.thinned) - T0.thinned))
}

save(file=args$bigDs, big.Ds)
save(file=args$bigNus, big.Nus)

# done.
