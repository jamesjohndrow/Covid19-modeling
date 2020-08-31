# -*- coding: utf-8 -*-
#
# Author: PB
# Maintainer(s): MG, KL, JJ
# License: (c) HRDAG 2020, GPL v2 or newer
#
# Covid19-undercount/analyze/Makefile
# -----------------------------------------------------------
#
.PHONY: all clean

STEP1 := \
	output/country_pops.RData \
	output/theta.png \
	output/theta.rds \
	output/theta2.png \
	output/theta2.rds \
	output/example_realization_w_nu.pdf \
	output/example_realization_w_nu.png \
	output/example_realization1.png

STEP3 := \
	output/CA_Death_forecast_0_01__.png \
	output/FL_Death_forecast_0_01__.png \
	output/NY_Death_forecast_0_01__.png \
	output/WA_Death_forecast_0_01__.png

STEP4 := \
	output/state_RT1.tex \
	output/state_R0.tex \
	output/state_undercount_daily.png \
	output/state_undercount.tex \
	output/state_undercount_alt.tex \
	output/state_R0_alt.tex \
	output/state_RT1_alt.tex

STEP5 := \
	output/gelman-rubin.tex \
	output/mpsrf.tex

STEP6 := \
	output/prior_table.tex

all: $(STEP1) $(STEP3) $(STEP4) $(STEP5) $(STEP6)

clean:
	-rm output/*

$(STEP1): \
		src/1-simulate_data.R \
		src/functions.R
	-mkdir output/
	Rscript --vanilla $<

# STEP2 - see 2-sir_death_mle.m; results from one run of this Matlab script are
# stored in frozen/ and will be used in proceeding steps

$(STEP3): \
		src/3-make_state_figures_matlab.R  \
		$(STEP1) \
		frozen/mcmc_*.mat
	Rscript --vanilla $<

$(STEP4): \
		src/4-make_state_tables_matlab.R \
		$(STEP3) \
		frozen/mcmc_*.mat
	Rscript --vanilla $<

$(STEP5): \
		src/5-MCMC_diagnostics.R \
		$(STEP3) \
		frozen/mcmc_*.mat
	Rscript --vanilla $<

$(STEP6): \
		src/6-monte_carlo_prior.R \
		$(STEP3) \
		frozen/mcmc_*.mat
	Rscript --vanilla $<

# done.
