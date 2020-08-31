#!/usr/bin/env Rscript --vanilla
# set expandtab ts=4 sw=4 ai fileencoding=utf-8
#
# Author: JJ
# Maintainer(s): MG, PB, JJ, KL
# License: (c) HRDAG 2020, GPL v2 or newer
#
# Covid19-undercount/analyze/src/4-make_state_tables.R
# -----------------------------------------------------------
#

if (!require('pacman')) {install.packages('pacman')}
pacman::p_load(xtable, lubridate, dplyr, ggplot2,
               argparse, glue, stringr, readr, tidyr, R.matlab)

stopifnot(str_ends(getwd(), "Covid19-modeling/main-analysis"))

set.seed(19481210)

ps <- c(0.002, 0.005, 0.01, 0.015, 0.025)
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
parser$add_argument("--state_RT1",
                    default="output/state_RT1.tex")
parser$add_argument("--state_R0",
                    default="output/state_R0.tex")
parser$add_argument("--state_undercount_graph",
                    default="output/state_undercount_daily.png")
parser$add_argument("--state_undercount_table",
                    default="output/state_undercount.tex")
parser$add_argument("--state_undercount_table_alt",
                    default="output/state_undercount_alt.tex")
parser$add_argument("--state_R0_alt",
                    default="output/state_R0_alt.tex")
parser$add_argument("--state_RT1_alt",
                    default="output/state_RT1_alt.tex")
args <- parser$parse_args()

args[["mcmc_results"]] <- "frozen/mcmc"

# https://github.com/nytimes/covid-19-data
# snapshot of NYT data as of 13 May 2020
state.dat <- read_csv(args$state_deaths)
abbrev <- read_csv(args$state_abbreviations)
pops <- read_csv(args$pops)

pops <- pops %>% select(State,Pop2018)
state.dat <- state.dat %>%
  inner_join(abbrev, by="state") %>%
  inner_join(pops, by=c("state"="State"))
state.dat <- state.dat %>%
  filter(is.element(abbreviation, states))

state.dat$p <- ps[1]
state.dat.tmp <- state.dat

for (p in ps[2:length(ps)]) {
  state.dat.tmp$p <- p
  state.dat <- rbind(state.dat, state.dat.tmp)
}

state.dat <- state.dat %>%
  arrange(p, abbreviation, date)

RT1.mat <- matrix("", np, nstate)
rownames(RT1.mat) <- ps
colnames(RT1.mat) <- states

RT1.mat.alt <- matrix("", 2, nstate)
rownames(RT1.mat.alt) <- c("mult_", "_alt")
colnames(RT1.mat.alt) <- states

R0.mat <- matrix("", np, nstate)
rownames(R0.mat) <- ps
colnames(R0.mat) <- states

R0.mat.alt <- matrix("", 2, nstate)
rownames(R0.mat.alt) <- c("mult_", "_alt")
colnames(R0.mat.alt) <- states

State.long <- cbind(states, c("California", "Florida", "New York", "Washington"))

ctr <- 0
for (p in ps) {
  for (state in states) {
    state.long <- State.long[State.long[,1]==state,2]

    if (p == 0.01) {
      scenarios <- c("_", "mult_", "_alt")
    } else {scenarios <- c("_")}


    for (scenario in scenarios) {

      ctr <- ctr + 1
      print(paste(p, state))
      p_str <- gsub('\\.', '_', p)
      if (exists('mcmc.out')) {rm(mcmc.out)}
      file <- glue("{args$mcmc_results}_{state.long}_{p_str}_{scenario}.mat")
      mcmc.out <- readMat(file)
      mcmc.out$X <- data.frame(mcmc.out$X)
      names(mcmc.out$X) <- c("T0S", "GAMMAR.INV", "R0", "PHIS", "ETA")

      T <- length(mcmc.out$SIR[,1,1])
      RT1 <- (mcmc.out$X$R0
              * mcmc.out$X$PHIS
              * mcmc.out$SIR[dim(mcmc.out$SIR)[1],1,]
              / (mcmc.out$X$ETA))
      R0 <-  mcmc.out$X$R0
      burn <- floor(length(mcmc.out$X$R0) / 2)
      nmc <- length(mcmc.out$X$R0)


      if(scenario == "_") {
        RT1.mat[which(p == ps), which(state == states)] <- paste(sprintf("%.2f", round(mean(RT1[burn:nmc]), 2)), " (",
                                                                 sprintf("%.2f", round(quantile(RT1[burn:nmc], 0.025), 2)), ",",
                                                                 sprintf("%.2f", round(quantile(RT1[burn:nmc], 0.975), 2)), ")", sep="")

        R0.mat[which(p == ps), which(state == states)] <- paste(sprintf("%.2f", round(mean(R0[burn:nmc]), 2)), " (",
                                                                sprintf("%.2f", round(quantile(R0[burn:nmc], 0.025), 2)), ",",
                                                                sprintf("%.2f", round(quantile(R0[burn:nmc], 0.975), 2)), ")", sep="")

        df.tmp <- data.frame(date=date(ymd("20200101") + seq(from=0, to=T - 1)),
                             state=state,
                             p=p,
                             S.mean=apply(mcmc.out$SIR[,1,burn:nmc], 1, mean),
                             S.025=apply(mcmc.out$SIR[,1,burn:nmc], 1, quantile, probs=0.025),
                             S.975=apply(mcmc.out$SIR[,1,burn:nmc], 1, quantile, probs=0.975))
        print(R0.mat)

      } else {
        RT1.mat.alt[scenario, which(state == states)] <- paste(sprintf("%.2f", round(mean(RT1[burn:nmc]),2)), " (",
                                                               sprintf("%.2f", round(quantile(RT1[burn:nmc], 0.025), 2)), ",",
                                                               sprintf("%.2f", round(quantile(RT1[burn:nmc], 0.975), 2)), ")", sep="")

        R0.mat.alt[scenario, which(state == states)] <- paste(sprintf("%.2f", round(mean(R0[burn:nmc]), 2)), " (",
                                                              sprintf("%.2f", round(quantile(R0[burn:nmc], 0.025), 2)), ",",
                                                              sprintf("%.2f", round(quantile(R0[burn:nmc], 0.975), 2)), ")", sep="")

        df.tmp <- data.frame(date=date(ymd("20200101") + seq(from=0, to=T - 1)),
                             state=state,
                             p=p,
                             S.mean=apply(mcmc.out$SIR[,1,burn:nmc], 1, mean),
                             S.025=apply(mcmc.out$SIR[,1,burn:nmc], 1, quantile, probs=0.025),
                             S.975=apply(mcmc.out$SIR[,1,burn:nmc], 1, quantile, probs=0.975))

      }

      state.dat.tmp <- state.dat %>%
        inner_join(df.tmp, by=c("abbreviation"="state", "p"="p", "date"="date")) %>%
        mutate(scenario=scenario)

      if (ctr == 1) {
        all.state.dat <- state.dat.tmp
      } else {
        all.state.dat <- rbind(all.state.dat, state.dat.tmp)
      }
    }
  }
}

print(xtable(RT1.mat,
             caption="Estimated posterior mean of $\\rho_{T}$, with 95 percent posterior credible interval shown in parentheses
             \\label{tab:RT1}"),
      file=args$state_RT1,
      size="footnotesize")

rownames(RT1.mat.alt) <- c("multiplier", "alternative theta")
print(xtable(RT1.mat.alt,
             caption="Estimated posterior mean of $\\rho_{T}$, with 95 percent posterior credible interval shown in parentheses.
Rows correspond to two alternative scenarios considered: death counts adjusted using an excess mortality multiplier and alternative value of $\\theta$ used.
In both cases, $p =0.01$.
             \\label{tab:RT1.alt}"),
      file=args$state_RT1_alt,
      size="footnotesize")

print(xtable(R0.mat,
             caption="Estimated posterior mean of $R_{0}$, with 95 percent posterior credible interval shown in parentheses.
             \\label{tab:R0}"),
      file=args$state_R0,
      size="footnotesize")

rownames(R0.mat.alt) <- c("multiplier", "alternative theta")
print(xtable(R0.mat.alt,
             caption="Estimated posterior mean of $R_{0}$, with 95 percent posterior credible interval shown in parentheses. Rows correspond to two alternative scenarios considered: death counts adjusted using an excess mortality multiplier and alternative value of $\\theta$ used.
In both cases, $p =0.01$.
             \\label{tab:R0.alt}"),
      file=args$state_R0_alt,
      size="footnotesize")

all.state.dat <- all.state.dat %>%
  mutate(S.mean=S.mean * Pop2018,
         S.025=S.025 * Pop2018,
         S.975=S.975 * Pop2018) %>%
  mutate(ucount.mean=(Pop2018 - S.mean) / cases,
         ucount.025=(Pop2018 - S.975) / cases,
         ucount.975=(Pop2018 - S.025) / cases,
         ucount.mean.lead5=(lag(Pop2018, 5) - lag(S.mean, 5)) / cases,
         ucount.mean.lead10=(lag(Pop2018, 10) - lag(S.mean, 10)) / cases)

all.state.dat.tmp <- all.state.dat %>% mutate(p=factor(p))

states.panels <- c("CA", "FL", "NY", "WA")

ulim <- all.state.dat %>%
  filter(date >= ymd("20200328")) %>%
  filter(scenario == "_") %>%
  summarise(ulim=max(ucount.mean))

plt <- ggplot(all.state.dat.tmp %>%
                filter(scenario == "_") %>%
                filter(date >= ymd("20200328"), abbreviation %in% states.panels),
              aes(x=date, y=ucount.mean, col=p)) +
  geom_line() +
  facet_wrap(~state) +
  ylim(0, ulim$ulim) +
  scale_y_log10() +
  theme(text = element_text(size=20)) +
  ylab("estimated undercount factor") +
  geom_ribbon(aes(ymin=ucount.025, ymax=ucount.975), alpha=0.2) +
  theme_bw() +
  theme(text=element_text(size=22))

plt
png(args$state_undercount_graph, width=900, height=600)
print(plt)
dev.off()

ucount_ests <- all.state.dat %>%
  filter(date == max(date)) %>%
  filter(scenario == "_") %>%
  select(abbreviation, p, ucount.mean, ucount.025, ucount.975)

ucount_mns <- ucount_ests %>%
  select(abbreviation, p, ucount.mean) %>%
  pivot_wider(names_from=abbreviation, values_from=ucount.mean) %>%
  select(-p)

ucount_025 <- ucount_ests %>%
  select(abbreviation, p, ucount.025) %>%
  pivot_wider(names_from=abbreviation, values_from=ucount.025) %>%
  select(-p)

ucount_975 <- ucount_ests %>%
  select(abbreviation, p, ucount.975) %>%
  pivot_wider(names_from=abbreviation, values_from=ucount.975) %>%
  select(-p)

ucount.mat <- matrix("", np, nstate)
rownames(ucount.mat) <- ps
colnames(ucount.mat) <- states

for (p in ps) {
  for (state in states) {
    row <- which(p == ps)
    col <- which(state == states)
    ucount.mat[row, col] <- paste(sprintf("%.2f", round(ucount_mns[row,col], 2)), " (",
                                  sprintf("%.2f", round(ucount_025[row,col], 2)), ",",
                                  sprintf("%.2f", round(ucount_975[row,col], 2)), ")", sep="")
  }
}

print(xtable(ucount.mat,
             caption="Estimated posterior mean of undercount, with 95 percent posterior credible interval shown in parentheses
             \\label{tab:undercount}"),
      file=args$state_undercount_table,
      size="footnotesize")

# make table for alternative scenarios
ucount_ests <- all.state.dat %>%
  filter(date == max(date)) %>%
  filter(scenario != "_") %>%
  select(abbreviation, p, scenario, ucount.mean, ucount.025, ucount.975)

ucount_mns <- ucount_ests %>%
  select(abbreviation, scenario, ucount.mean) %>%
  pivot_wider(names_from=abbreviation, values_from=ucount.mean) %>%
  select(-scenario)

ucount_025 <- ucount_ests %>%
  select(abbreviation, scenario, ucount.025) %>%
  pivot_wider(names_from=abbreviation, values_from=ucount.025) %>%
  select(-scenario)

ucount_975 <- ucount_ests %>%
  select(abbreviation, scenario, ucount.975) %>%
  pivot_wider(names_from=abbreviation, values_from=ucount.975) %>%
  select(-scenario)

ucount.mat <- matrix("", 2, nstate)
rownames(ucount.mat) <- c("mult_", "_alt")
colnames(ucount.mat) <- states

scenarios <- rownames(ucount.mat)
for (p in c("mult_", "_alt")) {
  for (state in states) {
    row <- which(p == scenarios)
    col <- which(state == states)
    ucount.mat[row, col] <- paste(sprintf("%.2f", round(ucount_mns[row,col], 2)), " (",
                                  sprintf("%.2f", round(ucount_025[row,col], 2)), ",",
                                  sprintf("%.2f", round(ucount_975[row,col], 2)), ")", sep='')
  }
}
rownames(ucount.mat) <- c("multiplier", "alternative theta")

print(xtable(ucount.mat,
             caption="Estimated posterior mean of undercount, with 95 percent posterior credible interval shown in parentheses for alternative scenarios : deaths adjusted according to excess mortality multipliers or using the alternative $\\theta$. For both alternative scenarios, $p = 0.01$.
             \\label{tab:undercount.alt}"),
      file=args$state_undercount_table_alt,
      size="footnotesize")

# done.
