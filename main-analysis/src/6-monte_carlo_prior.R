#!/usr/bin/env Rscript --vanilla
# set expandtab ts=4 sw=4 ai fileencoding=utf-8
#
# Author: JJ
# Maintainer(s): MG, PB, JJ, KL
# License: (c) HRDAG 2020, GPL v2 or newer
#
# Covid19-undercount/analyze/src/6-monte_carlo_prior.R
# -----------------------------------------------------------
#

if (!require('pacman')){install.packages('pacman')}
pacman::p_load(xtable, lubridate, dplyr, argparse, glue,
               stringr, readr, tidyr, R.matlab)

stopifnot(str_ends(getwd(), "Covid19-modeling/main-analysis"))

# make prior
p.prob <- dgamma(seq(0.002, 0.025, 0.001), 30, 3000)
p.prob <- p.prob / sum(p.prob)

ps <- seq(0.002, 0.025, 0.001)

# load data
np <- length(ps)
states <- c("NY", "CA", "FL", "WA")
nstate <- length(states)
states <- sort(states)

parser <- ArgumentParser()
parser$add_argument("--pops",
                    default="input/pop_data.csv")
parser$add_argument("--state_abbreviations",
                    default="input/state_abbrev.csv")
parser$add_argument("--state_deaths",
                    default="input/us-states.csv")
parser$add_argument("--prior_table",
                    default="output/prior_table.tex")

args <- parser$parse_args()

args[["mcmc_results"]] <- "frozen/mcmc"

# https://github.com/nytimes/covid-19-data
# snapshot of NYT data as of 13 May 2020
state.dat <- read_csv(args$state_deaths)
abbrev <- read_csv(args$state_abbreviations)
pops <- read_csv(args$pops)

pops <- pops %>% select(State, Pop2018)
state.dat <- state.dat %>%
  inner_join(abbrev, by="state") %>%
  inner_join(pops, by=c("state"="State"))
state.dat <- state.dat %>%
  filter(is.element(abbreviation, states))

State.long <- cbind(states, c("California", "Florida", "New York", "Washington"))
results.mat <- matrix(nrow=length(states), ncol=3)
rownames(results.mat) <- states
colnames(results.mat) <- c("R0", "RT1", "Undercount")

  for (state in states) {
    state.long <- State.long[State.long[,1] == state,2]

      all.data <- list()
      ctr <- 1
      for(p in ps){
        p_str <- gsub('\\.', '_', p)
        if (exists('mcmc.out')) {rm(mcmc.out)}
        file <- glue("{args$mcmc_results}_{state.long}_{p_str}__.mat")
        mcmc.out <- readMat(file)
        mcmc.out$X <- data.frame(mcmc.out$X)
        names(mcmc.out$X) <- c("T0S", "GAMMAR.INV", "R0", "PHIS", "ETA")

        T <- length(mcmc.out$SIR[,1,1])
        RT1 <- (mcmc.out$X$R0
                * mcmc.out$X$PHIS
                * mcmc.out$SIR[dim(mcmc.out$SIR)[1],1,]
                / (mcmc.out$X$ETA))
        R0 <-  mcmc.out$X$R0

        S <- mcmc.out$SIR[121, 1,]

        stats <- list(RT1=RT1, R0=R0, S=S)

        all.data[[ctr]] <- stats
        ctr <- ctr + 1
      }

        sample.ps <- sample(1:length(ps), nrow(mcmc.out$X), prob=p.prob, replace=TRUE)

        combined.RT1 <- matrix(nrow=length(sample.ps), ncol=1)
        combined.R0 <- matrix(nrow=length(sample.ps), ncol=1)
        combined.S <- matrix(nrow=length(sample.ps), ncol=1)

        for(i in 1:length(sample.ps)){
          combined.RT1[i] <- unlist(all.data[[sample.ps[i]]]$RT1[i])
          combined.R0[i] <- unlist(all.data[[sample.ps[i]]]$R0[i])
          combined.S[i] <- unlist(all.data[[sample.ps[i]]]$S[i])
        }

        state.pop <- pops[pops$State == state.long,]$Pop2018
        final.cases <- filter(state.dat, state == state.long, date == as.Date("2020-04-30"))$cases
        ucounts <- (state.pop - combined.S * state.pop) / final.cases

        burn <- floor(length(combined.RT1) / 2)
        nmc <- length(combined.RT1)

        results.mat[state, 2] <- paste(sprintf("%.2f", round(mean(combined.RT1[burn:nmc]),2)), " (",
              sprintf("%.2f", round(quantile(combined.RT1[burn:nmc], 0.025), 2)), ",",
              sprintf("%.2f", round(quantile(combined.RT1[burn:nmc], 0.975), 2)), ")", sep="")

        results.mat[state, 1] <- paste(sprintf("%.2f", round(mean(combined.R0[burn:nmc]), 2)), " (",
              sprintf("%.2f", round(quantile(combined.R0[burn:nmc], 0.025), 2)), ",",
              sprintf("%.2f", round(quantile(combined.R0[burn:nmc], 0.975), 2)), ")", sep="")

        results.mat[state, 3] <- paste(sprintf("%.2f", round(mean(ucounts[burn:nmc])), 2), " (",
                                       sprintf("%.2f", round(quantile(ucounts[burn:nmc], 0.025), 2)), ",",
                                       sprintf("%.2f", round(quantile(ucounts[burn:nmc], 0.975), 2)), ")", sep="")

  }

print(xtable(results.mat,
             label="tab:prior_table",
             caption="$R_0$, $\\rho_T$, and the estimated undercount (95\\% posterior credible intervals) calculated by averaging samples over a prior on $p$."),
      file=args$prior_table)

# done.
