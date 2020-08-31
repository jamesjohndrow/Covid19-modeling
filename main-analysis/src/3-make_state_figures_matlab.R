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

if (!require('pacman')) {install.packages('pacman')}
pacman::p_load(argparse, ggplot2, dplyr, reshape2, stringr, glue,
               readr, tidyr, lubridate, future, future.apply, R.matlab)

stopifnot(str_ends(getwd(), "Covid19-modeling/main-analysis"))

set.seed(19481210)

##-----GLOBAL vars
states <- c("CA", "FL", "NY", "WA")
ps <- c(0.002, 0.005, 0.01, 0.015, 0.025)
multiplier <- data.frame(states=states, multiplier=c(4,3,2,1))
State.long <- cbind(states, c("California", "Florida", "New York", "Washington"))

for (State in states) {

  for (p in ps) {

    if (p == 0.01) {
      scenarios <- c("_", "mult_", "_alt")
    } else {scenarios <- c("_")}

    for (scenario in scenarios){
      if (scenario == "mult_") {
        m <- multiplier[multiplier$state == State,2]
      } else {m <- 1}

      n.days.future <- 20
      p_str <- gsub("\\.", "_", p)
      state.long <- State.long[State.long[,1] == State,2]

      getargs <- function() {
        runtag <- glue("{tolower(State)}")
        parser <- ArgumentParser()
        parser$add_argument("--funs", default="src/functions.R")
        parser$add_argument("--mcmc_output",
                            default=glue("frozen/mcmc_{state.long}_{p_str}_{scenario}.mat"))
        parser$add_argument("--theta", default="output/theta.rds")
        parser$add_argument("--theta2", default="output/theta2.rds")
        parser$add_argument("--pops", default="input/pop_data.csv")
        parser$add_argument("--state_abbreviations", default="input/state_abbrev.csv")
        parser$add_argument("--state_deaths", default="input/us-states.csv")
        parser$add_argument("--state_forecast", default=glue("output/{State}_Death_forecast_{p_str}_{scenario}.png"))
        parser$add_argument("--state_fit", default=glue("output/{State}_{p_str}.png"))
        parser$parse_args()
      }
      args <- getargs()

      source(args$funs)

      # https://github.com/nytimes/covid-19-data
      dat <- read_csv(args$state_deaths)  # snapshot of NYT data as of 13 May 2020

      pops <- read_csv(args$pops)
      abbrev <- read_csv(args$state_abbreviations)
      pops <- pops %>% inner_join(abbrev, by=c("State"="state"))

      N <- pops$Pop2018[pops$abbreviation == State]

      deaths <- dat %>%
        left_join(abbrev) %>%
        filter(abbreviation == State) %>%
        select(date, deaths) %>%
        arrange(date) %>%
        mutate(death_increase=deaths - lag(deaths)) %>%
        replace_na(list(death_increase=0))

      n.to.add <- as.numeric(min(deaths$date) - ymd("20200101"))
      Ds <- c(rep(0, n.to.add), deaths$death_increase) * m

      if(scenario!= "_alt") {
        theta <- readRDS(args$theta)
      } else {
        theta <-  readRDS(args$theta2)
      }

      # time period used to fit data, 01/01/2020 - 04/30-2020 (inclusive on both ends)
      T <- as.integer(as.Date("2020-04-30") - as.Date("2020-01-01")) +  1

      mcmc.out <- readMat(args$mcmc_output)
      mcmc.out$X <- data.frame(mcmc.out$X)
      names(mcmc.out$X) <- c("T0S", "GAMMAR.INV", "R0", "PHIS", "ETA")

      nmc <- nrow(mcmc.out$X) # check for state
      burn <- 25000   # check for state
      n.thin <- 1000
      thinned <- round(seq(burn, nmc, length.out=n.thin))

      # prepare data for plots
      big.Ds <- matrix(0, T + n.days.future, length(thinned))

      # R0 = beta/gamma
      BETA.thinned <- mcmc.out$X$R0[thinned] / mcmc.out$X$GAMMAR.INV[thinned]
      PHI.thinned <- mcmc.out$X$PHIS[thinned]
      GAMMAR.thinned <- 1 / mcmc.out$X$GAMMAR.INV[thinned]
      T0.thinned <- mcmc.out$X$T0S[thinned]
      SIS.thinned <- t(mcmc.out$SIR[T, ,thinned])
      NUS.thinned <- mcmc.out$NU[thinned,] # only take the first T, some extras included as forecasts
      ETAS.thinned <-  mcmc.out$X$ETA[thinned]

      param.list <- cbind(
        BETA.thinned, # col 1
        GAMMAR.thinned, # col 2
        T0.thinned, # col 3
        SIS.thinned, # col 4 &  5
        PHI.thinned, # col 6
        ETAS.thinned, # col 7
        NUS.thinned # col 8 through...
      )

      Nus.predict <- future_apply(param.list, 1, function(x){return(c(x[8:length(x)]*N,
                                                                      simulate.sir(x[1] * x[6] , # new beta = beta*phi*eta
                                                                                   x[2] * x[7], # new gamma = gamma*eta
                                                                                   x[5], # I0
                                                                                   N, # N
                                                                                   1, # start point
                                                                                   1 + n.days.future, # number of days in future
                                                                                   x[4] # S0
                                                                      )$Nus[-1]))})

      plot(Nus.predict[,1])
      Ds.predict <- future_apply(Nus.predict, 2, function(x){sample.Ds.ma(x, theta, p)})

      plot(Ds.predict[,1])
      q.hi <- apply(Ds.predict, 1, quantile, probs=0.975)
      q.lo <- apply(Ds.predict, 1, quantile, probs=0.025)
      med <- apply(Ds.predict, 1, median)

      avg <- stats::filter(Ds, filter=c(1, 1, 1, 1, 1, 1, 1) / 7, method="convolution", sides=2)
      last.day <- length(Ds)
      avg[last.day - 2] <- mean(Ds[last.day:(last.day - 6)])
      avg[last.day - 1] <- mean(Ds[last.day:(last.day - 5)])
      avg[last.day] <- mean(Ds[last.day:(last.day - 4)])

      df.Ds <- data.frame(avg=c(avg, rep(NA, length(med) - length(Ds))),
                          Ds=c(Ds, rep(NA, length(med) - length(Ds))),
                          med=med,
                          q.hi=q.hi,
                          q.lo=q.lo,
                          day=date(ymd("20200101") + seq(0, T + n.days.future-1)))

      plt <- ggplot(filter(df.Ds, day >= as.Date("2020-02-15")),
                    aes(x=day, y=med)) +
        geom_line(linetype = "dashed") +
        theme_bw() +
        geom_line(aes(x=day, y=avg), color = "red") +
        geom_line(aes(x=day, y=Ds), color="gray", alpha=0.6) +
        geom_vline(xintercept=as.Date("2020-04-30")) +
        theme(text=element_text(size=24)) +
        geom_ribbon(data=filter(df.Ds, day >= as.Date("2020-02-15")),
                    aes(ymin=q.lo, ymax=q.hi),
                    alpha=0.2) +
        xlab("date") +
        ylab("D") +
        ggtitle(pops$State[pops$abbreviation == State])

      png(args$state_forecast, width=600, height=400)
      print(plt)
      dev.off()

    }
  }
}

# done.
