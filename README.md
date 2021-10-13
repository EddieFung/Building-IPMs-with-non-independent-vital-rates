# Building integral projection models with non-independent vital rates

This repository contains example R and nimble code to support the manuscript "Building integral projection models with non-independent vital rates". 

Soay sheep data is available from https://doi.org/10.1111/j.1600-0706.2012.00035.x on Coulson (2012), and North Atlantic Ocean Index (NAO) data is available on https://crudata.uea.ac.uk/cru/data/nao/nao.dat from the Climate Research Unit at the University of East Anglia.

You will have to download the Soay sheep data from the above link yourself to run the provided .R files. 

For convenience, we copied the NAO data that we used in the manuscript in this repository. “NAO” comprises all the monthly NAO data in [1986,1997] while “NAOc” comprises all the monthly NAO data in [1969,2019].



## Code Description:

fn1.R: this .R includes the nimble code and required libraries of the MCMC algorithms for estimating vital rate parameters.

fn2.R: this .R includes the code of MCEM to approximate the mle in the correlated random individual effect model (M3). This approximated mle will be used as the initialization of MCMC to reduce the computational cost for convergence.

fn3.R: this .R includes the code and required libraries to approximate log lambda given that you have run run1.R to obtain the posterior samples of parameters.

run1.R: this is the first .R file you should run. This .R include code for estimating vital rate parameters by MCMC on nimble. The posterior samples of vital rate parameters will be stored in a new directory name "sample". The computational time is a few minutes per model, except model D1b and D3, which may run around a day

run2.R: this is the second .R file you should run, given that you already ran run1.R to obtain the posterior samples. This .R include code for approximate the log lambda. The estimates of log lambda will be stored in a new directory name "result". The computational time depends on the number of available core. It took us less than 1 hour to finish with 5 cores
