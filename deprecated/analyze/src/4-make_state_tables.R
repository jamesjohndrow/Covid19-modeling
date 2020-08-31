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

if (!require('pacman')){install.packages('pacman')}
pacman::p_load(xtable, lubridate, dplyr, ggplot2,
               argparse, glue, stringr, readr, tidyr)

stopifnot(str_ends(getwd(), "Covid19-modeling/analyze"))

set.seed(19481210)

ps <- c(0.005, 0.01, 0.015)
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
                    default="input/us-states-20200418.csv")
parser$add_argument("--state_RT1",
                    default="output/state_RT1.tex")
parser$add_argument("--state_R0",
                    default="output/state_R0.tex")
parser$add_argument("--state_undercount_graph",
                    default="output/state_undercount_daily.png")
parser$add_argument("--state_undercount_table",
                    default="output/state_undercount.tex")

args <- parser$parse_args()

args[["mcmc_stub"]] <- "output/mcmc_sim"

# reported cases as of 18 April 2020, for undercount analysis
state.dat <- read_csv(args$state_deaths)
abbrev <- read_csv(args$state_abbreviations)
pops <- read_csv(args$pops)

pops <- pops %>% select(State,Pop2018)
state.dat <- state.dat %>%
  inner_join(abbrev, by="state") %>%
  inner_join(pops, by=c("state"="State"))
state.dat <- state.dat %>% 
  filter(is.element(abbreviation, states))

state.dat$p <- 0.01
state.dat.tmp <- state.dat
state.dat.tmp$p <- 0.005
state.dat <- rbind(state.dat, state.dat.tmp)

state.dat.tmp$p <- 0.015
state.dat <- rbind(state.dat, state.dat.tmp)

state.dat <- state.dat %>%
  arrange(p, abbreviation, date)

RT1.mat <- matrix('', np, nstate)
rownames(RT1.mat) <- ps
colnames(RT1.mat) <- states

R0.mat <- matrix('', np, nstate)
rownames(R0.mat) <- ps
colnames(R0.mat) <- states

ctr <- 0
for (p in ps) {
  for (state in states) {
    ctr <- ctr + 1
    print(paste(p, state))
    p_str <- gsub('\\.', '_', p)
    if (exists('mcmc.out')) {rm(mcmc.out)}
    # load(paste(path, '/output/mcmc_sim_', tolower(state), '_', p_str, '.RData', sep=''))
    load(glue("{args$mcmc_stub}_{tolower(state)}_{p_str}.RData"))
    T <- length(mcmc.out$SIS[,1,1])
    RT1 <- (mcmc.out$BETA 
            * mcmc.out$PHIS 
            * mcmc.out$SIS[dim(mcmc.out$SIS)[1],1,]
            / (mcmc.out$GAMMAR))
    R0 <-  mcmc.out$BETA / mcmc.out$GAMMAR
    burn <- floor(length(mcmc.out$BETA) / 2)
    nmc <- length(mcmc.out$BETA)
    
    RT1.mat[which(p == ps), which(state == states)] <- paste(sprintf("%.2f", round(mean(RT1[burn:nmc]),2)), " (",
                                                             sprintf("%.2f", round(quantile(RT1[burn:nmc], 0.025), 2)), ",",
                                                             sprintf("%.2f", round(quantile(RT1[burn:nmc], 0.975), 2)), ")", sep='')
    
    R0.mat[which(p == ps), which(state == states)] <- paste(sprintf("%.2f", round(mean(R0[burn:nmc]), 2)), " (",
                                                            sprintf("%.2f", round(quantile(R0[burn:nmc], 0.025), 2)), ",",
                                                            sprintf("%.2f", round(quantile(R0[burn:nmc], 0.975), 2)), ")", sep='')
    
    df.tmp <- data.frame(date=date(ymd("20200101") + seq(from=0, to=T - 1)),
                         state=state,
                         p=p,
                         S.mean=apply(mcmc.out$SIS[,1,burn:nmc], 1, mean),
                         S.025=apply(mcmc.out$SIS[,1,burn:nmc], 1, quantile, probs=0.025),
                         S.975=apply(mcmc.out$SIS[,1,burn:nmc], 1, quantile, probs=0.975))
    
    state.dat.tmp <- state.dat %>%
      inner_join(df.tmp, by=c("abbreviation"="state", "p"="p", "date"="date"))
    if (ctr==1) {
      all.state.dat <- state.dat.tmp
    } else {
      all.state.dat <- rbind(all.state.dat, state.dat.tmp)
    }
    
  }
}

all.state.dat <- all.state.dat %>% 
  mutate(S.mean=S.mean * Pop2018, S.025=S.025 * Pop2018, S.975=S.975 * Pop2018) %>%
  mutate(ucount.mean=(Pop2018 - S.mean) / cases,
         ucount.025=(Pop2018 - S.975) / cases,
         ucount.975=(Pop2018 - S.025) / cases,
         ucount.mean.lead5=(lag(Pop2018, 5) - lag(S.mean, 5)) / cases,
         ucount.mean.lead10=(lag(Pop2018, 10) - lag(S.mean, 10)) / cases)

print(xtable(RT1.mat,
             caption="Estimated posterior mean of $R_{T_1}$, with 95 percent posterior credible interval shown in parentheses 
             \\label{tab:RT1}"),
      file=args$state_RT1,
      # file=paste(path,'/output/','state_RT1.tex',sep=''),
      size='footnotesize')

print(xtable(R0.mat,
             caption="Estimated posterior mean of $R_{0}$, with 95 percent posterior credible interval shown in parentheses 
             \\label{tab:R0}"),
      # file=paste(path,'/output/','state_R0.tex',sep=''),
      file=args$state_R0,
      size='footnotesize')

states.panels <- c("CA", "FL", "NY", "WA")

ulim <- all.state.dat %>%
  filter(date >= ymd("20200328")) %>% 
  summarise(ulim=max(ucount.mean))
plt <- ggplot(all.state.dat %>% filter(date >= ymd("20200328"), abbreviation %in% states.panels),
              aes(x=date, y=ucount.mean, col=factor(p))) + 
  geom_point() + 
  geom_line() +
  facet_wrap(~state, scales='free') + ylim(0, ulim$ulim)

# png('output/state_undercount_daily.png', width=900, height=600)
png(args$state_undercount_graph, width=900, height=600)
print(plt)
dev.off()

ucount_ests <- all.state.dat %>% 
  filter(date == max(date)) %>% 
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

ucount.mat <- matrix('', np, nstate)
rownames(ucount.mat) <- ps
colnames(ucount.mat) <- states

for (p in ps) {
  for (state in states) {
    row <- which(p == ps)
    col <- which(state == states)
    ucount.mat[row, col] <- paste(sprintf("%.2f", round(ucount_mns[row,col], 2)), " (",
                                  sprintf("%.2f", round(ucount_025[row,col], 2)), ",",
                                  sprintf("%.2f", round(ucount_975[row,col], 2)), ")", sep='')
  }
}

print(xtable(ucount.mat,
             caption="Estimated posterior mean of undercount, with 95 percent posterior credible interval shown in parentheses 
             \\label{tab:undercount}"),
      # file=paste(path,'/output/','state_undercount.tex',sep=''),
      file=args$state_undercount_table,
      size='footnotesize')

# done.
