#
# Authors: JJ,  KL
# Maintainers: PB, MG
# Copyright:
# =========================================
# Covid19-modeling/analyze/src/2a-fit_model_to_real_data.R

#this code fits a model to the data. It is faster than the other script, which we left in case someone is currently working off of it.

if (!require('pacman')) {install.packages('pacman')}
pacman::p_load(ggplot2, dplyr, tidyr, readr, lubridate,
               future, future.apply, here, glue, argparse)

parser <- ArgumentParser()
parser$add_argument('--functions', default=here::here('analyze/src/functions_faster.R'))
parser$add_argument('--theta', default=here::here('analyze/output/theta.Rds'))
parser$add_argument('--pops', default=here::here('analyze/input/country_pops.RData'))
parser$add_argument('--worldometers', default=here::here('analyze/input/worldcounts_deaths_us.csv'))
parser$add_argument('--p', default=0.01)
parser$add_argument('--country', default='us')

args <- parser$parse_args()

args[['mcmc']] <- here::here(glue('analyze/output/mcmc_run_{args$p}.RData'))
args[['bigDs']] <- here::here(glue('analyze/output/bigDs_{args$p}.RData'))
args[['bigNus']] <- here::here(glue('analyze/output/bigNus_{args$p}.RData'))

source(args$functions)
#this code fits a model to the data

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
N <- pops[rownames(pops) == args$country]


#how long since mitigation went into place
days.since.mitigation <- 10
#how many days into the future do we want to predict
n.days.future <- 100


nmc <- 50000
burn <- 10000
n.thin <- 1000
thinned <- round(seq(burn, nmc, length.out=n.thin))
disp.int <- 1000


# worldometers data entered manually from website
dat <- read_csv(args$worldometers)
dat <- dat %>% mutate(date=mdy(date))
Ds <- dat$deaths_us
Ds <- c(rep(0, 14), Ds)


# In a "normal" SIR model, the R0 is beta/gamma
# The average days to recovery is 1/gamma
# Let's choose
gammar <- 1 / 7
zeta <- c(4, 25) # limits of uniform prior on gammar^(-1)
beta <- gammar * 2.5
xi <- c(1.4, 3.9) # prior on beta | gammar  is Uniform(xi[1]*gammar,xi[2]*gammar) on this interval


#length of data time series
T <- length(Ds)

#proportion of initial infections
I <- 1 / N

#time of first infection
T0 <- 15

# cases to consider for percent decrease in R0
per_decrease <- seq(from=0, to=0.8, by=0.1)
n.cases <- length(per_decrease)


##----- run mcmc
#find a reasonable place to initialize parameters

sim1 <- SIR(pars=list(beta=beta, gamma=gammar),
            init=c(S=1 - 1 / N, I=1 / N, R=0),
            time=seq(from=1, to=T, by=1))

plot(sim1$results[,2], type='l')
lines(sim1$results[,3], col=2)
lines(sim1$results[,4], col=3)


S <- sim1$results[,2] * N
Nus0 <- S[1:(length(S) - 1)] - S[2:length(S)]
Nus <- rep(0, T)
Nus[floor(T0):T] <- Nus0[1:(T - T0 + 1)]
sim.Ds <- sample.Ds.ma(Nus, theta, p)

#sim1 <- simulate.sir(N,beta,gammar,I,T+max.time + 1,max.time,theta,p,T0)
plot(sim.Ds)
plot(sim.Ds[1:T], Ds, xlim=c(0, max(sim.Ds[1:T], Ds)))
abline(0, 1)

##----- run MCMC
mcmc.out <- run.mcmc(nmc=nmc, theta=theta, beta=beta, gammar=gammar, T=T, I=I,
                     Ds=Ds, Nus=Nus, xi=xi, zeta=zeta, disp.int=disp.int,
                     k.adapt=5000, c0=c(0.01^2, 0.01^2, 0.75^2),
                     c1=0.5 * (2.38 / 3), T0=T0, zeta.T0=c(1, 30))


save(file=args$mcmc, mcmc.out)

thinned <- seq(burn, nmc, length.out=n.thin)

NUS <- mcmc.out$NUS
SIS <- mcmc.out$SIS

#do things look alright?
death.samples <- apply(NUS[thinned,], 1, sample.Ds.ma, theta=theta, p=p)
plot(apply(death.samples, 1, mean), Ds,
     ylim=c(0, max(Ds) * 1.2), xlim=c(0, max(Ds) * 1.2))
abline(0, 1)
#yes, posterior predictive mean for Ds is good!

#prepare data for plots
big.Ds <-  array(0, dim=c(T + n.days.future, length(thinned), n.cases))
big.Nus <- array(0, dim=c(T + n.days.future, length(thinned), n.cases))

BETA.thinned <- mcmc.out$BETA[thinned]
GAMMAR.thinned <- mcmc.out$GAMMAR[thinned]
T0.thinned <- mcmc.out$T0S[thinned]
SIS.thinned <- t(mcmc.out$SIS[T - days.since.mitigation,,thinned])
NUS.thinned <- NUS[thinned,]

sim1 <- simulate.sir(BETA.thinned[1],
                     GAMMAR.thinned[1],
                     SIS.thinned[1,2],
                     N,
                     1,
                     n.days.future + days.since.mitigation + 1,
                     SIS.thinned[1,1])
# looks good
plot(c(NUS.thinned[1,1:(T - days.since.mitigation - 1)], sim1$Nus))

plan(multisession)
param.list <- split(cbind(BETA.thinned, GAMMAR.thinned, T0.thinned, SIS.thinned,
                          NUS.thinned), seq(1:n.thin))

for (j in 1:n.cases) {
  print(j)
  phi <- 1 - per_decrease[j]
  sir.predict <- future_sapply(param.list,
                               function(x) {c(x[6:(T-days.since.mitigation-1+5)],
                                              simulate.sir(x[1] * phi,
                                                           x[2], 
                                                           x[5],
                                                           N,
                                                           1,
                                                           1 + days.since.mitigation + n.days.future,
                                                           x[4])$Nus)})
  
  Ds.predict <- future_apply(sir.predict, 2, function(x){sample.Ds.ma(x, theta, p)})
  big.Nus[,,j] <- sir.predict
  big.Ds[,,j] <- Ds.predict
}


plot(big.Ds[,1,6])

save(file=args$bigDs, big.Ds)
save(file=args$bigNus, big.Nus)

# done.
