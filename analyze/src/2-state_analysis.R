#!/usr/bin/env Rscript --vanilla
# set expandtab ts=4 sw=4 ai fileencoding=utf-8
#
# Author: KL JJ
# Maintainer(s): PB, MG 
# License: (c) HRDAG 2020, GPL v2 or newer
#
# Covid19-undercount/analyze/src/4-state_analysis.R
# -----------------------------------------------------------
#

if (!require('pacman')){install.packages('pacman')}
pacman::p_load(argparse, readr, tidyr, dplyr, stringr, 
               ggplot2, lubridate, glue, pracma)

stopifnot(str_ends(getwd(), "Covid19-modeling/analyze"))

set.seed(19481210)

#Restaurant: CA 15, NY 16, WA 16, FL 20, MI 16, NJ 16
#School: CA 13, NY 18, WA 17, MI 16, FL 17, NJ 18

##-----GLOBAL vars
##-----set initial conditions to get reasonable starting point for MCMC
States.sd <- data.frame(State=c("NY", "NJ", "CA", "WA", "MI", "FL"),
                        date.sd=c("20200318", "20200318", "20200315",
                                  "20200317", "20200316", "20200320"),
                        T0=rep(30, 6),
                        R0=c(2.8, 2.5, 2.3, 2.3, 2.4, 2.5),
                        phi=rep(0.4, 6))

# State <- "CA"
smooth.deaths <- TRUE

states <- c("CA", "FL", "NY", "WA")
ps <- c(0.005, 0.01, 0.015)

for (State in states) {
  for (p in ps) {
    
    date.sd <- States.sd$date.sd[States.sd$State == State]
    T0 <- States.sd$T0[States.sd$State == State]
    R0.init <- States.sd$R0[States.sd$State == State]
    phi <- States.sd$phi[States.sd$State == State]
    # p <- 0.01
    p_str <- gsub('\\.', '_', p)
    
    print(paste(State, p_str))
    
    getargs <- function() {
      # sets file paths in one place, could be used in command line call
      # works in RStudio, too
      parser <- argparse::ArgumentParser()
      parser$add_argument("--funs", default="src/functions.R")
      parser$add_argument("--pops", default="input/pop_data.csv")
      parser$add_argument("--state_abbreviations", default="input/state_abbrev.csv")
      parser$add_argument("--theta", default="output/theta.rds")
      parser$add_argument("--state_deaths", default="input/us-states-20200418.csv")
      parser$parse_args()
    }
    args <- getargs()
    
    args[["mcmc_output"]] <- glue("output/mcmc_sim_{tolower(State)}_{p_str}.RData")
    args[["posterior_hists"]] <- glue("output/posterior_hists_{tolower(State)}_{p_str}.png")
    args[["posterior_trace"]] <- glue("output/posterior_trace_{tolower(State)}_{p_str}.png")
    
    
    source(args$funs)
    
    # https://github.com/nytimes/covid-19-data
    dat <- read_csv(args$state_deaths)  # snapshot of NYT data as of 18 April
    dat$date <- ymd(as.character(dat$date))
    
    pops <- read_csv(args$pops)
    abbrev <- read_csv(args$state_abbreviations)
    pops <- pops %>% inner_join(abbrev, by=c("State"="state"))
    
    N <- pops$Pop2018[pops$abbreviation == State]
    
    deaths <- dat %>% left_join(abbrev) %>%
      filter(abbreviation == State) %>% 
      select(date, deaths) %>%
      arrange(date) %>% 
      mutate(death_increase=deaths - lag(deaths)) %>%
      replace_na(list(death_increase=0))
    
    n.to.add <- as.numeric(min(deaths$date) - ymd("20200101"))
    Ds <- c(rep(0, n.to.add), deaths$death_increase)
    
    T1 <- as.numeric(n.to.add + ymd(date.sd) - min(deaths$date))
    I <- 1 / N
    gammar <- 1 / 6
    zeta <- c(mean=6.4, sd=1.5, a=3.4, b=9.4) # limits of uniform prior on gammar^(-1)
    beta <- gammar * R0.init
    xi <- c(mean=2.5, sd=1.5, a=1, b=4) # prior on beta | gammar  is Uniform(xi[1]*gammar,xi[2]*gammar) on this interval
    zeta.T0 <- c(1, 60)
    zeta.phi <- c(0.01, 0.99)
    
    theta <- readRDS(args$theta)
    T <- length(Ds)
    t.max <- length(Ds)
    i.max <- length(theta)
    max.time <- length(theta) - 1
    
    n.thin <- 1000
    
    sim1 <- simulate.sir.timebreak(beta, gammar, 1 / N, N, T0, T, T1, phi)
    Nus <- sim1$Nus
    
    Ds.sim <- calculate.death.means(Nus, theta, p)
    
    plot(Ds.sim[1:T])
    plot(Ds.sim[1:T], Ds, xlim=c(0, max(Ds.sim[1:T], Ds)))
    abline(0, 1) #looks pretty good!
    
    # everywhere but new york, new jersey
    if (!is.element(State, c("NY", "NJ"))) {
      c0 <- c(0.005^2, 0.005^2, 0.75^2, 0.005^2)
      c1 <- (2.38 / 4)
      nmc <- 50000
      burn <- 25000
      k.adapt <- 10000
    } else { # new york, new jersey only
      c0 <- c(0.002^2, 0.002^2, 0.75^2, 0.002^2)
      c1 <- (2.38 / 4) * 0.5 
      nmc <- 100000
      burn <- 50000
      k.adapt <- 20000
    }
    thinned <- round(seq(burn, nmc, length.out=n.thin))
    
    adaptive <- TRUE
    
    # if you aren't sure you have the tuning parameters right, then set to TRUE to see trace plots
    # the plotting slows things down, so if you are confident it will mix well, then don't make the plots as it runs
    plotting <- FALSE
    disp.int <- 5000
    
    mcmc.out <- run.mcmc.state(nmc, theta, beta, gammar, T, I, Ds, Nus, xi, zeta,
                               disp.int, k.adapt=k.adapt, c0=c0, c1=c1, T0=T0, 
                               zeta.T0=zeta.T0, phi, T1, zeta.phi, plotting)
    
    save(file=args$mcmc_output, mcmc.out)
    load(file=args$mcmc_output)
    
    df <- data.frame(beta=mcmc.out$BETA,
                     gammar=mcmc.out$GAMMAR,
                     T0=mcmc.out$T0S,
                     phi=mcmc.out$PHI,
                     SIt=mcmc.out$SIS[dim(mcmc.out$SIS)[1],1,])
    df <- df %>% mutate(R0=beta / gammar, RT1=(R0 * phi) * SIt)
    df$iter <- seq(nmc)
    
    df <- df %>% pivot_longer(-iter)
    
    png(args$posterior_hists, width=600, height=400)
    print({ggplot(df %>% filter(iter >= burn, name %in% c("R0", "RT1", "T0")),
                  aes(x=value)) +
        geom_histogram() + 
        facet_wrap(~name, scales='free') +
        ggtitle(pops$State[pops$abbreviation == State])})
    dev.off()
    
    png(args$posterior_trace, width=600, height=400)
    print({ggplot(df %>% filter(name %in% c("beta", "gammar", "R0")),
                  aes(x=iter, y=value)) +
        geom_point(size=0.2) + 
        facet_wrap(~name, scales='free')})
    dev.off()
    
  }
}

# done.
