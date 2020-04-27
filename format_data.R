
.libPaths("~/Rlib")
setwd("~/Desktop/Code/laplace_approximation/Script")

library(rstan)

## Ovarian data n=54, p=1536
ovarian <- read.csv("aki_disease_data/ovarian.csv", header=FALSE)

y <- ovarian[,1537]
x <- ovarian[,1:1536]
## faster with just 100 covariates (including the most relevant one)
# x <- ovarian[,1437:1536]

n = nrow(x);
p = ncol(x);

## regularized horseshoe prior
p0 <- 5
sigma <- sqrt(1/mean(y)/(1-mean(y)))
sigma <- 1
tau0 <- p0/(p-p0)*sigma/sqrt(n)
## data
data <- list(n = n, d = ncol(x), x = x, y = y,
             nu_local = 1, nu_global = 1, scale_global = tau0,
             scale_icept=5, slab_scale=2, slab_df=100)

with(data, stan_rdump(ls(data), file.path("data", "ovarian_reduced.data.R")))

## Bernoulli logistic model with compound GLM function
# fitb1 <- stan("bernoulli_logit_glm_rhs.stan", data=data, chains=4,
#               iter=3000, warmup=1000,
#               control=list(adapt_delta = 0.8, max_treedepth=10), cores=4)

## Gaussian linear model with compound GLM function, adapt_delta=0.8
fitg1 <- stan("normal_id_glm_rhs.stan", data=data, chains=4,
              iter=3000, warmup=1000,
              control=list(adapt_delta = 0.8, max_treedepth=10), cores=4)

## Gaussian linear model with compound GLM function, adapt_delta=0.999
fitg1 <- stan("normal_id_glm_rhs.stan", data=data, chains=4,
              iter=3000, warmup=1000,
              control=list(adapt_delta = 0.999, max_treedepth=15), cores=4)

## Gaussian linear model with analytic integration over weights
fitg2 <- stan("normal_id_dotprod_rhs.stan", data=data, chains=4,
              iter=3000, warmup=1000,
              control=list(adapt_delta = 0.8, max_treedepth=10), cores=4)

###############################################################################
## Additional code (by Charles) to output the code for the prostate space.

# Read-in prostate.RData from the aki_disease_data directory.
load("aki_disease_data/prostate.RData")
y <- y[, 1]  # convert matrix to vector
x <- x[, 2501:2700]  # sub-sample (including 1816 and 2586)
n = nrow(x)
p = ncol(x)
p0 <- 5
sigma <- sqrt(1 / mean(y) / (1 - mean(y)))
sigma <- 1  # CHECK -- what's the point of the above line?
tau0 <- p0 / (p - p0) * sigma / sqrt(n)

## data
data <- list(n = n, d = ncol(x), x = x, y = y,
             nu_local = 1, nu_global = 1, scale_global = tau0,
             scale_icept = 5, slab_scale = 2, slab_df = 100)

with(data, stan_rdump(ls(data), file.path("data", "prostate_200.data.R")))

