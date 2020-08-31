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
	output/example_realization_w_nu.pdf \
	output/example_realization_w_nu.png \
	output/example_realization1.png

STEP2 := \
	output/mcmc_sim_ca_0_01.RData \
	output/mcmc_sim_ca_0_015.RData \
	output/mcmc_sim_ca_0_05.RData \
	output/mcmc_sim_fl_0_01.RData \
	output/mcmc_sim_fl_0_015.RData \
	output/mcmc_sim_fl_0_05.RData \
	output/mcmc_sim_ny_0_01.RData \
	output/mcmc_sim_ny_0_015.RData \
	output/mcmc_sim_ny_0_05.RData \
	output/mcmc_sim_wa_0_01.RData \
	output/mcmc_sim_wa_0_015.RData \
	output/mcmc_sim_wa_0_05.RData \
	output/posterior_hists_ca_0_01.png \
	output/posterior_hists_ca_0_015.png \
	output/posterior_hists_ca_0_05.png \
	output/posterior_hists_fl_0_01.png \
	output/posterior_hists_fl_0_015.png \
	output/posterior_hists_fl_0_05.png \
	output/posterior_hists_ny_0_01.png \
	output/posterior_hists_ny_0_015.png \
	output/posterior_hists_ny_0_05.png \
	output/posterior_hists_wa_0_01.png \
	output/posterior_hists_wa_0_015.png \
	output/posterior_hists_wa_0_05.png \
	output/posterior_trace_ca_0_01.png \
	output/posterior_trace_ca_0_015.png \
	output/posterior_trace_ca_0_05.png \
	output/posterior_trace_fl_0_01.png \
	output/posterior_trace_fl_0_015.png \
	output/posterior_trace_fl_0_05.png \
	output/posterior_trace_ny_0_01.png \
	output/posterior_trace_ny_0_015.png \
	output/posterior_trace_ny_0_05.png \
	output/posterior_trace_wa_0_01.png \
	output/posterior_trace_wa_0_015.png \
	output/posterior_trace_wa_0_05.png

STEP3 := \
	output/ca_Death_forecast.png \
	output/fl_Death_forecast.png \
	output/ny_Death_forecast.png \
	output/wa_Death_forecast.png \
	output/ca_0_01.png \
	output/ca_0_015.png \
	output/ca_0_05.png \
	output/fl_0_01.png \
	output/fl_0_015.png \
	output/fl_0_05.png \
	output/ny_0_01.png \
	output/ny_0_015.png \
	output/ny_0_05.png \
	output/wa_0_01.png \
	output/wa_0_015.png \
	output/wa_0_05.png

STEP4 := \
	output/state_RT1.tex \
	output/state_R0.tex \
	output/state_undercount_daily.png \
	output/state_undercount.tex

all: $(STEP1) $(STEP2) $(STEP3) $(STEP4)

clean:
	rm output/*


$(STEP1): \
		src/1-simulate_data.R \
		src/functions.R
	Rscript --vanilla $<


$(STEP2): \
		src/2-state_analysis.R \
		src/functions.R \
		input/pop_data.csv \
		input/state_abbrev.csv \
		output/theta.rds \
		input/us-states-20200418.csv
	Rscript --vanilla $<


$(STEP3): \
		src/3-make_state_figures.R \
		src/functions.R \
		input/pop_data.csv \
		input/state_abbrev.csv \
		output/theta.rds \
		input/us-states-20200418.csv
	Rscript --vanilla $<


$(STEP4): \
		src/4-make_state_tables.R \
		input/state_abbrev.csv \
		input/pop_data.csv \
		input/us-states-20200418.csv
	Rscript --vanilla $<

# done.
