#!/usr/bin/env Rscript --vanilla
# set expandtab ts=4 sw=4 ai fileencoding=utf-8
#
# Author: JJ
# Maintainer(s): MG, PB, JJ, KL
# License: (c) HRDAG 2020, GPL v2 or newer
#
# Covid19-undercount/analyze/src/5-MCMC_diagnostics.R
# -----------------------------------------------------------
#

if (!require("pacman")) {install.packages("pacman")}
pacman::p_load(argparse, coda, R.matlab, glue, xtable, stringr)

stopifnot(str_ends(getwd(), "Covid19-modeling/main-analysis"))

parser <- ArgumentParser()
parser$add_argument("--gelman_rubin", default="output/gelman-rubin.tex")
parser$add_argument("--mpsrf", default="output/mpsrf.tex")

args <- parser$parse_args()

args[["mcmc_results"]] <- "frozen/mcmc"

states <- c("California", "Florida", "New York", "Washington")
diags.point <- matrix(0, nrow=5, ncol=length(states))
diags.CI.upper <- matrix(0, nrow=5, ncol=length(states))

multivariate.psrf <- NULL

for (j in 1:length(states)) {
  state <- states[j]

  chains <- list()
  for (i in 1:5) {
    dat <- readMat(glue("{args$mcmc_results}_{state}_0_01___{i}.mat"))
    chains[[i]] <- dat
  }


  Param.list <- list()
  for (i in  1:5) {
    Param.list[[i]] <- mcmc(chains[[i]]$X)
  }
  Param.list <- mcmc.list(Param.list)
  Param.diag <- gelman.diag(Param.list, multivariate=TRUE)

  diags.point[,j] <- Param.diag$psrf[,1]
  diags.CI.upper[,j] <- Param.diag$psrf[,2]

  multivariate.psrf[j] <- Param.diag$mpsrf

}

final.table <- matrix(paste0(round(diags.point, 3), " (" , round(diags.CI.upper, 3), ")"),
                      nrow=5,
                      ncol=4)
rownames(final.table) <- c("T0", "R0", "1/Gamma", "Phi", "Eta")
colnames(final.table) <- states

print(xtable(final.table,
             caption="Gelman-Rubin convergence diagnostic point estimate (upper 95\\% confidence interval) for parameters of model calculated with 5 independent runs of MCMC",
             label="tab:gr"),
      file=args$gelman_rubin)

tmp <- data.frame(multivariate.psrf)
names(tmp) <- c( "MPSRF")
rownames(tmp) <- states
print(xtable(tmp,
             caption = "Multivariate Gelman-Rubin PSRF",
             label="tab:mpsrf"),
      file=args$mpsrf)

# done.
