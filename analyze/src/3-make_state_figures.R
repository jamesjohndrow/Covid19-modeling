#!/usr/bin/env Rscript --vanilla
# set expandtab ts=4 sw=4 ai fileencoding=utf-8
#
# Author: JJ
# Maintainer(s): MG, PB, JJ, KL
# License: (c) HRDAG 2020, GPL v2 or newer
#
# Covid19-undercount/analyze/src/3-make_state_figures.R
# -----------------------------------------------------------
#

if (!require('pacman')){install.packages('pacman')}
pacman::p_load(argparse, ggplot2, dplyr, reshape2, stringr, 
               glue, readr, tidyr, lubridate, future, future.apply)

stopifnot(str_ends(getwd(), "Covid19-modeling/analyze"))

set.seed(19481210)

states <- c("CA", "FL", "NY", "WA")
ps <- c(0.005, 0.01, 0.015)

##-----GLOBAL vars
for (State in states) {
  for (p in ps) {
    n.days.future <- 20
    p_str <- gsub("\\.","_",p)
    
    getargs <- function() {
      runtag <- glue("{tolower(State)}")
      parser <- ArgumentParser()
      parser$add_argument("--funs", default="src/functions.R")
      parser$add_argument("--mcmc_output",
                          default=glue("output/mcmc_sim_{runtag}_{p_str}.RData"))
      parser$add_argument("--theta", default="output/theta.rds")
      parser$add_argument("--pops", default="input/pop_data.csv")
      parser$add_argument("--state_abbreviations", default="input/state_abbrev.csv")
      parser$add_argument("--state_deaths", default="input/us-states-20200418.csv")
      parser$add_argument("--state_forecast", default=glue("output/{State}_Death_forecast.png"))
      parser$add_argument("--state_fit", default=glue("output/{State}_{p_str}.png"))
      parser$parse_args()
    }
    args <- getargs()
    
    source(args$funs)
    
    pops <- read_csv(args$pops)
    abbrev <- read_csv(args$state_abbreviations)
    pops <- pops %>% inner_join(abbrev, by=c("State"="state"))
    
    N <- pops$Pop2018[pops$abbreviation == State]
    
    # https://github.com/nytimes/covid-19-data
    dat <- read_csv(args$state_deaths)  # snapshot of NYT data as of 18 April
    dat$date <- ymd(as.character(dat$date))
    
    deaths <- dat %>% 
      left_join(abbrev) %>%
      filter(abbreviation == State) %>% 
      select(date, deaths) %>%
      arrange(date) %>% 
      mutate(death_increase=deaths-lag(deaths)) %>%
      replace_na(list(death_increase=0))
    
    n.to.add <- as.numeric(min(deaths$date) - ymd("20200101"))
    Ds <- c(rep(0, n.to.add), deaths$death_increase)
    
    theta <- readRDS(args$theta)  
    T <- length(Ds)
    
    load(args$mcmc_output)
    
    nmc <- nrow(mcmc.out$SIS) # check for state
    burn <- 25000   # check for state
    n.thin <- 1000
    thinned <- round(seq(burn, nmc, length.out=n.thin))
    
    # prepare data for plots
    big.Ds <- matrix(0, T + n.days.future, length(thinned))
    
    BETA.thinned <- mcmc.out$BETA[thinned]
    PHI.thinned <- mcmc.out$PHIS[thinned]
    GAMMAR.thinned <- mcmc.out$GAMMAR[thinned]
    T0.thinned <- mcmc.out$T0S[thinned]
    SIS.thinned <- t(mcmc.out$SIS[T, ,thinned])
    NUS.thinned <- mcmc.out$NUS[thinned, ]
    
    param.list <- cbind(
      BETA.thinned,
      GAMMAR.thinned,
      T0.thinned,
      SIS.thinned,
      PHI.thinned,
      NUS.thinned
    )
    
    Nus.predict <- future_apply(param.list, 1, function(x){return(c(x[7:length(x)],
                                                                    simulate.sir(x[1] * x[6],
                                                                                 x[2],
                                                                                 x[5],
                                                                                 N,
                                                                                 1,
                                                                                 1 + n.days.future,
                                                                                 x[4])$Nus[-1]))})
    
    
    Ds.predict <- future_apply(Nus.predict, 2, function(x){sample.Ds.ma(x, theta, p)})
    
    plot(Ds.predict[,1])
    q.hi <- apply(Ds.predict, 1, quantile, probs=0.975)
    q.lo <- apply(Ds.predict, 1, quantile, probs=0.025)
    med <- apply(Ds.predict, 1, median)
    
    q.hi[1:T] <- Ds
    q.lo[1:T] <- Ds
    med[1:T] <- Ds
    
    df.Ds <- data.frame(med=med,
                        q.hi=q.hi,
                        q.lo=q.lo,
                        day=date(ymd("20200101") + seq(T + n.days.future)))
    
    
    plt <- ggplot(df.Ds, aes(x=day, y=med)) +
      geom_line() +
      theme_bw() +
      theme(text=element_text(size=24)) +
      geom_ribbon(data=df.Ds,
                  aes(ymin=q.lo,
                      ymax=q.hi),
                  alpha=0.2) +
      xlab("date") +
      ylab("D") +
      ggtitle(pops$State[pops$abbreviation == State])
    
    png(args$state_forecast, width=600, height=400)
    print(plt)
    dev.off()
    
    
    q.hi <- apply(Ds.predict, 1, quantile, probs=0.975)
    q.lo <- apply(Ds.predict, 1, quantile, probs=0.025)
    df.Ds$q.hi <- q.hi
    df.Ds$q.lo <- q.lo
    
    plt <- ggplot(df.Ds[1:T,],aes(x=day,y=med)) + 
      geom_line() +
      theme_bw() +
      theme(text=element_text(size=24)) +
      geom_ribbon(data=df.Ds[1:T,],
                  aes(ymin=q.lo,
                      ymax=q.hi),
                  alpha=0.2) +
      xlab("date") +
      ylab("D") +
      ggtitle(pops$State[pops$abbreviation == State])
    
    png(args$state_fit, width=600, height=400)
    print(plt)
    dev.off()
  }
}
# done.
