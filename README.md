# laplace_manuscript

Supplementary code for the the article _Hamiltonian Monte Carlo using an embedded Laplace approximation_. When running the R scripts, adjust the working, data, model, delivery, and other directories to your setting.

### Installation 

In order to run the scripts, you need to install cmdStan. 
The embedded Laplace approximation functions presented in the paper 
have not been released into Stan (yet).
To use the prototype, you need to (i) get the C++ code for the functions
and (ii) expose the functions to the Stan language.
This can be done using the following commands on a terminal window
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

It then remains to install cmdStanR (https://mc-stan.org/cmdstanr/),
which is a lightweight wrapper of cmdStan for R.

### Overview

The script `bernoulli_logit.r` can be used to run the computer experiment for the general linear regressions and SKIM models (with and without the embedded Laplace approximation), on both the prostate and ovarian cancer data sets. Once the output is saved into an RStan object, the plots can be generated using the cluster analysis scripts.

To run the disease map model, use `disease_map2.r`. You can generate a subsample of the data with the `read_data.r` script.

The Stan files are stored under the `model3` directory.
