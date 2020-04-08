# -*- coding: utf-8 -*-
#
# Author: PB
# Maintainer(s): PB, JJ, KL, MG
# License: (c) HRDAG 2020, GPL v2 or newer
#
# Covid19-modeling/analyze/Makefile
# -----------------------------------------------------------
#
.PHONY: all clean

all: output/example_realization.png

clean:
	rm output/*

output/example_realization.png: \
		src/1-simulate_data_for_examples.R \
		src/functions.R \
		output/theta.Rds
	Rscript --vanilla $<

# done.
