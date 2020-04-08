# Covid19-modeling

This repository contains data, code, and LaTeX files for the paper "Estimating SARS-CoV-2-positive Americans using deaths-only data" on [arXiv](https://arxiv.org/abs/2004.02605) by James E. Johndrow, Kristian Lum, Maria Gargiulo, and Patrick Ball. 

The abstract:

   We fit a Bayesian model to data on the number of deaths attributable to COVID-19 with the goal of estimating the number of infected individuals. Our model links an underlying Susceptible Infectious Removed (SIR) model of disease dynamics to observed deaths via a time-to-death distribution informed by previous studies. This allows us to actually fit a statistical model to the data, unlike many epidemiological studies in which the SIR model parameters are simply "calibrated" to obtain outputs that look similar to the real data. The main outputs of our model are estimates of the number of infections currently, as well as forecasts of the number of infections and deaths under various scenarios for the effectiveness of social distancing measures. All of our outputs have attached Bayesian credible intervals. An important conclusion is that the confirmed case counts greatly underestimate the total number of infected individuals.

And the goals:

    The goal of this work is not to provide specific forecasts of the infected population and likely deaths, although our initial estimates here fit the observed patterns of deaths. Accurate forecasts will require location-specific measures of containment and treatment efficacy, as well as age- and comorbidity-specific infected fatality rates. We don’t have these data at present, but our model could incorporate them as better information becomes available. Our paper offers a modeling approach using minimal but probably-good data, describes a likelihood and priors, is fitted to data, and is underpinned by a widely-used epidemiological model that is designed to approximate the real dynamics of disease spread.


<!-- done --> 
