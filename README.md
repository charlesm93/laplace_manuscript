# Supplementary Material: Code

Code for the the article _Hamiltonian Monte Carlo using an adjoint-differentiated Laplace approximation: Bayesian inference for latent Gaussian models and beyond_. 
`supplementary.pdf` contains a proof of theorem 1, and details about our computer experiments.
The rest of the directory contains the code used for the computer experiments.


### Installation 

In order to run the R scripts, you need to install cmdStan. 
The embedded Laplace approximation functions presented in the paper 
have not been released into Stan (yet),
but are available on GitHub on the appropriate branches.
To use the prototype, you need to (i) get the C++ code for the functions
and (ii) expose the functions to the Stan language.
This can be done using the following commands on a (linux based) terminal window
```
# 1) download Stan's transpiler with the relevant function signatures 
git clone https://github.com/stan-dev/stanc3.git
cd stanc3
git checkout try-laplace_approximation

# 2) Download cmdstan and get the relevant C++ code
git clone https://github.com/stan-dev/cmdstan.git
cd cmdstan
make stan-update
cd stan/lib/stan_math
git checkout try-laplace_approximation

# 3) Build the binaries by compiling an example
cd ..
cd ..
cd ..
STANC3=../stanc3 make examples/bernoulli/bernoulli
``` 
More recent prototypes of the code are available as we develop them.
Currently, the include improved function signatures and can be installed by replacing the checked out branch in `stanc3` with `try-laplace_approximation2`, and similarly, the checked out branch in `stan_math` with `try-laplace_approximation2`.
Note that this is _not_ the prototype we used in our computer experiment

It then remains to install cmdStanR (https://mc-stan.org/cmdstanr/),
which is a lightweight wrapper of cmdStan for R.

### R scripts

When running the R scripts, adjust the working, data, model, delivery, and other directories to your setting at the top of the file.

With the following libraries, you can all the R scripts:
```
library(boot)
library(cmdstanr)
library(directlabels)
library(dplyr)
library(latex2exp)
library(parallel)
library(posterior)
library(rstan)
library(scales)
library(stringr)
library(tidyr)
```
The required libraries for each individual script is specified at the top of the script.

The script `bernoulli_logit.r` can be used to run the computer experiment for the general linear model (GLM) with a horseshoe prior and SKIM model (with full HMC or with the embedded Laplace approximation), on the prostate cancer data set. The output is saved into a StanFit object. Depending on your setting (say, whether you are using a personal computer with only 4 cores or a cluster) you may need one or two runs to generate 6 chains.

To run the disease map model, use `disease_map2.r`. You can generate a subsample of the data with the `read_data.r` script.

All the Stan files are in the `model3` directory.

`cluster_analysis.r` contains the code to plot and compare the samples, and other metrics of interest, for the GLM with a horseshoe prior. The code reads in the StanFit objects saved when using `bernoulli_logit.r`. Next, `cluster_analysis_v2.r` is an improved version; as written the code applies to saved fits for the SKIM, but can be adjusted for other models.

`error_analysis.r` contains code to compute and plot the error for all three models. To do the error analysis, you need to generate benchmarks. In the paper, these benchmarks are generated using long runs of MCMC. There are practical difficulties when reading the benchmark for the GLM, since we have ~12,000 parameters over 98,000 samples. To avoid memory overflow in R, the benchmarks are created separately, using `read_cmdstan_output.r`.

Finally, the analysis for ADVI on all three models is in `advi_sample.r`.
