# Estimating the number of SARS-CoV-2 infections and the impact of mitigation policies in the United States

This repository contains R and Matlab code for "Estimating the number of SARS-CoV-2 infections and the impact of mitigation policies in the United States" by James Johndrow, Patrick Ball, Maria Gargiulo, and Kristian Lum, forthcoming in the [Harvard Data Science Review](https://hdsr.mitpress.mit.edu/).

> **Abstract**: Knowledge of the number of individuals who have been infected with the novel coronavirus SARS-CoV-2 and the extent to which attempts for mitigation by executive order have been effective at limiting its spread are critical for effective policy going forward. Directly assessing prevalence and policy effects is complicated by the fact that case counts are unreliable. In this paper, we present a model for using death-only data---in our opinion, the most stable and reliable source of COVID-19 information---to estimate the underlying epidemic curves. Our model links observed deaths to an SIR model of disease spread via a likelihood that accounts for the lag in time from infection to death and the infection fatality rate. We present estimates of the extent to which confirmed cases in the United States undercount the true number of infections, and analyze how effective social distancing orders have been at mitigating or suppressing the virus. We provide analysis for four states with significant epidemics: California, Florida, New York, and Washington.
>
> **Keywords**: COVID-19, SARS-CoV-2, compartmental model, SIR model, Bayesian

**Note:** this repository requires [git lfs](https://git-lfs.github.com/)

### Some brief notes
* We use Matlab to fit our model because of the strengths of Matlab's ODE solver. We understand that not everyone has access to Matlab, so we provide the results of one run (the one this paper is based on) in `main-analysis/frozen` and reference those result files in steps 3-6.
* Everything except step 2 (the piece written in Matlab) can be executed sequentially via `make`.
* The code will also run interactively from RStudio as long as the working directory is correctly set (and there's a unit test that will fail if it's not).

<!-- done. -->
