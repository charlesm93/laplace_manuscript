
## Draft code: do not source file.

modelName <- "bernoulli_logit_glm_rhs"

setwd("~/Desktop/Code/laplace_approximation/Script/")
.libPaths("~/Rlib")

library(rstan)
library(boot)
source("tools/cmdStanTools.r")

###############################################################################
chains <- 1
seed <- 12345
file <- file.path("deliv_cmdstan", paste0(modelName, ".", seed, ".", chains,
                                          ".csv"))

chain_output <- read.table(file, skip = 37, sep = ",")

stanfit <- read_stan_csv(file)
check_div(stanfit)

pars = c("tau", "caux", "lambda[1816]", "lambda[2586]")
parms <- rstan::extract(stanfit, pars = pars)

benchmark <- c(mean(log(parms$tau)),
               mean(log(parms$caux)),
               mean(log(parms$`lambda[1816]`)),
               mean(log(parms$`lambda[2586]`))
               )

write.csv(benchmark, file = "deliv_cmdstan/mean.csv")

print(get_elapsed_time(stanfit))
 
###############################################################################

# Get the log mean
# read in benchmark csv output
modelName <- "bernoulli_logit_glm_rhs"
seed <- "12345"
nChain <- 8

pars <- c("tau", "caux")
pars_lambda <- c("lambda[2586]", "lambda[1816]", "lambda[4960]",
                 "lambda[4238]", "lambda[4843]", "lambda[3381]", "lambda[4647]")

pars_f <- c("f[1]", "f[2]", "f[99]", "f[100]")

n_iter <- 12000
n_parm <- length(pars_lambda)
samples_lambda <- array(NA, c(n_iter, n_parm, nChain))

for (i in 2:nChain) {  # should be 1:nChain
  chain_output <- read_stan_csv(file.path("deliv_cmdstan",
                    paste0(modelName, ".", seed, ".", i, ".csv")))

  # compute the mean of the log parameters
  chain_samples <- rstan::extract(chain_output, pars = pars)
  mean_log_sample <- 
    lapply(chain_samples, mean_log <- function(x) { mean(log(x)) })
  
  chain_samples_lambda <- rstan::extract(chain_output, pars = pars_lambda)
  mean_log_sample_lambda <-
    lapply(chain_samples_lambda, mean_log <- function(x) { mean(log(x))} )

  # store all samples for quantile estimation
  samples_lambda[, , i] <- rstan::extract(chain_output, pars = pars_lambda,
                                          permuted = F)

  # compute the 90th quantile for the lambda parameters
  # NOTE - we cannot average estimates of the quantiles!!
  # quantile.90 <-
  #   lapply(chain_samples_lambda, quant <- function(x)
  #     { quantile(log(x), probs = 0.9) })

  # compute the probability log(lambda) is above a certain threshold.
  threshold <- 2.5
  prob_threshold <- function(x, threshold) {
    sum(x > threshold) / length(x)
  }

  prob <- lapply(chain_samples_lambda, prob <- function(x) { 
    prob_threshold(x, threshold) } )

  # get the mean for the parameters (i.e. no log scale)
  chain_samples_f <- rstan::extract(chain_output, pars = pars_f)
  mean_logit <-
    lapply(chain_samples_f, mean_bernoulli <- function(x)
      { mean(inv.logit(x)) })

  write.csv(c(mean_log_sample, mean_log_sample_lambda, mean_logit),
            file = paste0("deliv_cmdstan/", modelName, ".logmean.",
                          i, ".csv"))

  write.csv(prob, file = paste0("deliv_cmdstan/", modelName, ".prob",
                                i, ".csv"))

  # write.csv(quantile.90, file = paste0("deliv_cmdstan/", modelName,
  #                                      ".90quantile", i, ".csv"))

  print(paste0(i, " / ", nChain))
}

# compute the quantiles using samples from all 8 chains
samples_lambda_all <- samples_lambda[, , 1]
for (i in 2:nChain) samples_lambda_all <- rbind(samples_lambda_all,
                                                samples_lambda[, , i])

n_pars <- length(pars_lambda)
quantiles_90 <- rep(NA, n_pars)
for (i in 1:n_pars) {
  quantiles_90[i] <- quantile(log(samples_lambda_all[, i]), probs = 0.9)
}

write.csv(quantiles_90, file = paste0("deliv_cmdstan/", modelName,
                                      ".90quantile.csv"))
