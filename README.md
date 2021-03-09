# Process-Dependent-IPMs

This repository contains example R and nimble code to support the manuscript "Process Dependent Integral Projection Models". 

Soay sheep data is available from https://doi.org/10.1111/j.1600-0706.2012.00035.x on Coulson (2012), and North Atlantic Ocean Index (NAO) data is available on https://crudata.uea.ac.uk/cru/data/nao/nao.dat from the Climate Research Unit at the University of East Anglia.

You will have to download the Soay sheep data from the above link yourself to run the provided .R files. 

For NAO data, we provide you the required data in this repository. “NAO” comprises all the monthly NAO data in [1986,1997] while “NAOc” comprises all the monthly NAO data in[1969,2019].



## Code Description:

fn1.R: it includes the nimble code to estimate the parameters of interest

fn2.R: it includes the code of MCEM to approximate the mle in the correlated random individual effect model. This approximated mle will be used as the initialization of MCMC to reduce the computational cost for convergence.

fn3.R: it includes the code to approximate log lambda, given the posterior samples of parameters. 

run1.R: this is the first .R file you should run. You will obtain the posterior samples of parameters in a new directory named "sample".

run2.R: this is the second .R you should run, after you run run1.R to obtain the posterior samples. Tha approximated log lambda will be stored in a new directory named "result". 
