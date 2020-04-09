# -*- coding: utf-8 -*-
#
# Author: PB
# Maintainer(s): JJ, KL, MG
# License: (c) HRDAG 2020, GPL v2 or newer
#
# Covid19-modeling/analyze/Makefile
# -----------------------------------------------------------
#
.PHONY: all clean

all: output/example_realization.png \
	output/mcmc_run_0.01.RData \
	output/bigDs_0.01.RData \
	output/bigNus_0.01.RData

clean:
	rm output/*

STEP_3 :=
		output/trace_plots_fractional_us_0.01.png \
		output/post_hists_fractional_us_0.01.png \
		output/num_deaths_cases_fractional_us_0.01.png \
		output/num_newly_infected_cases_fractional_us_0.01.png \
		output/Predicted_vs_actual_fractional_us_0.01.png \
		output/R0_fractional_us_0.01.png \
		output/infected_10_days_ago_fractional_us_0.01.png \
		output/newly_infected_10_days_ago_fractional_us_0.01.png \
		output/time_to_recovery_fractional_us_0.01.png

output/example_realization.png: \
		src/1-simulate_data_for_examples.R \
		src/functions.R \
		output/theta.Rds
	Rscript --vanilla $<

output/mcmc_run_0.01.RData output/bigDs_0.01.RData output/bigNus_0.01.RData: \
		src/2-fit_model_to_real_data.R \
		src/functions.R \
		output/theta.Rds \
		input/country_pops.RData \
		input/worldcounts_deaths_us.csv
	Rscript --vanilla $<

$(STEP_3): \
		src/3-make_figures.R \
		src/functions.R \
		output/theta.Rds \
		input/country_pops.RData \
		input/worldcounts_deaths_us.csv \
		output/mcmc_run_0.01.RData \
		output/bigDs_0.01.RData \
		output/bigNus_0.01.RData
	Rscript --vanilla $<

# done.
