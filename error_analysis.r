
## Error analysis

setwd("~/Desktop/Code/laplace_approximation/Script")
.libPaths("~/Rlib/")
library(rstan)
library(dplyr)
library(tidyr)
library(scales)
library(boot)
library(directlabels)
library(latex2exp)
source("tools/errorTools.r")
source("tools/analysisTools.r")

################################################################
## Analysis of the disease map
modelName <- "disease_map"
fit_standard <- readRDS(file = file.path("deliv", modelName,
                        paste0(modelName, "Fit.Rsave")))

fit_benchmark <- readRDS(file = file.path("deliv", modelName,
                         paste0(modelName, "100long_Fit.Rsave")))

modelName <- "disease_map_ela"
fit_laplace <- readRDS(file = file.path("deliv", modelName,
                       paste0(modelName, "Fit.Rsave")))

pars <- c("alpha", "rho", "theta[1]", "theta[2]")

# it's worth examining the trace plots
traceplot(fit_standard, pars = pars)
traceplot(fit_laplace, pars = pars)

# Compare the Monte Carlo estimates
summary_standard <- summary(fit_standard, pars = pars)[[1]]
summary_benchmark <- summary(fit_benchmark, pars = pars)[[1]]
summary_laplace <- summary(fit_laplace, pars = pars)[[1]]

# plot Monte Carlo estimate with error
mcse_plot(summary_standard, summary_laplace, summary_benchmark, pars)

# Plot error
iter <- seq(from = 1, to = 500, by = 1)
nChains <- 4
samples_standard <- rstan::extract(fit_standard, pars, permuted = FALSE)  # [, 1, ]
samples_benchmark <- rstan::extract(fit_benchmark, pars, permuted = FALSE)
samples_laplace <- rstan::extract(fit_laplace, pars, permuted = FALSE)  # [, 1, ]

benchmark <- rep(NA, length(pars))
for (i in 1:length(pars)) benchmark[i] <- 
  mean(samples_benchmark[, , i])

error_plot(samples_standard, samples_laplace, benchmark,
           pars = pars, one_chain = F, average_chain = T) +
  scale_y_continuous(trans = "log")


# Plot the error as a function of time.
# The time for each iteration is the mean warmup time,
# plus the estimated time to produce n samples,
# based on the mean total sampling time.
include_warmup <- F
time_standard <- iteration_time(fit_standard, iter, include_warmup)
time_laplace <- iteration_time(fit_laplace, iter, include_warmup)

error_plot(samples_standard, samples_laplace, benchmark,
           pars = pars, plot_time = T,
           time_standard = time_standard,
           time_laplace = time_laplace,
           one_chain = F, average_chain = T) +
  # scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log")


#####################################################################
## Analysis for SKIM model

modelName <- "skim_logit"
fit_standard <- readRDS(file = file.path("deliv", modelName, "cluster",
                                         paste0(modelName, "Fit.Rsave")))
benchmarks <- 1:6
fit_benchmark <- list()
for (i in benchmarks) {
  fit_benchmark[[i]] <- readRDS(file = file.path("deliv", modelName,
                                paste0(modelName, "long_", benchmarks[i],
                                       "_Fit.Rsave")))
}
fit_benchmark <- sflist2stanfit(fit_benchmark)

modelName <- "skim_logit_ela"
fit_laplace <- readRDS(file = file.path("deliv", modelName, "cluster",
                                        paste0(modelName, "Fit.Rsave")))

pars <- c("caux", "lambda[81]", "lambda[86]", "tau", "xi")
# pars <- c("lambda[81]")
summary_standard <- summary(fit_standard, pars = pars)[[1]]
summary_benchmark <- summary(fit_benchmark, pars = pars)[[1]]
summary_laplace <- summary(fit_laplace, pars = pars)[[1]]

# plot Monte Carlo estimate with error
mcse_plot(summary_standard, summary_laplace, summary_benchmark, pars) +
  scale_y_continuous(trans = "log")

## Plot error
# TEST -- see if the plots are more interpretable with one chain
iter <- seq(from = 1, to = 2000, by = 1)
nChains <- 6  # 6
samples_standard <- rstan::extract(fit_standard, pars, permuted = FALSE)
samples_benchmark <- rstan::extract(fit_benchmark, pars, permuted = FALSE)
samples_laplace <- rstan::extract(fit_laplace, pars, permuted = FALSE)

log_parm <- T
benchmark <- rep(NA, length(pars))

if (!log_parm) {
  for (i in 1:length(pars)) benchmark[i] <- mean(samples_benchmark[, , i])
} else {
  for (i in 1:length(pars)) benchmark[i] <- 
                                     mean(log(samples_benchmark[, , i]))
}

error_plot(samples_standard, samples_laplace, benchmark,
           pars = pars, one_chain = F, average_chain = T,
           log_parm = log_parm, f = c(mean)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="lr")

# Plot for Neurips
# first create benchmark for plots
quantile_90 <- function(x) { quantile(x, probs = 0.9) }
benchmark <- c(mean(log(samples_benchmark[, , 1])),
               quantile_90(log(samples_benchmark[, , 2])),
               quantile_90(log(samples_benchmark[, , 3])),
               mean(log(samples_benchmark[, , 4])),
               mean(log(samples_benchmark[, , 5])))

include_warmup <- F
time_standard <- iteration_time(fit_standard, iter, include_warmup)
time_laplace <- iteration_time(fit_laplace, iter, include_warmup)
f <- c(mean, quantile_90, quantile_90, mean, mean)

parm_label <- c(TeX("$E( \\log \ c_{aux})$"), 
                TeX("$Q_{90}( \\log \\lambda_{2516})$"),
                TeX("$Q_{90}( \\log \\lambda_{2586})$"), 
                TeX("$E( \\log \\tau)$"),
                TeX("$E( \\log \\chi)$"))

error_plot(samples_standard, samples_laplace, benchmark,
           pars = pars, one_chain = F, average_chain = T,
           log_parm = T, f = f, plot_time = T,
           time_laplace = time_laplace,
           time_standard = time_standard,
           parm_label = parm_label) +
  theme(legend.position = "none", text = element_text(size = 25)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^floor(x)),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^floor(x)),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="l")


# Check which covariates are softly selected (i.e. inspect 90th quantile)
pars <- c("lambda")
log_lambda <- log(rbind(rstan::extract(fit_benchmark, pars)$lambda))
select_lambda(log_lambda, quant = 0.9, n_select = 6)
# 86 160 179  81 120 151

log_lambda <- log(rbind(rstan::extract(fit_standard, pars)$lambda))
select_lambda(log_lambda, quant = 0.9, n_select = 6)
# 86 160 179  81 120 151

log_lambda <- log(rbind(rstan::extract(fit_laplace, pars)$lambda))
select_lambda(log_lambda, quant = 0.9, n_select = 6)
# 86 179 160  81 120  48

# Examine some of the high scoring lambdas.
pars <- c("lambda[86]", "lambda[160]", "lambda[179]", "lambda[81]")

# First compute the 90th quantile of log lambda
samples_benchmark <- rstan::extract(fit_laplace, pars)
quantile_90_log <- function(x) { quantile(log(x), probs = 0.9) }
benchmark <- lapply(samples_benchmark, quantile_90_log)
benchmark <- unlist(benchmark, use.names = F)

quantile_90 <- function(x) { quantile(x, probs = 0.9) }

samples_standard <- rstan::extract(fit_standard, pars, permuted = F)
samples_laplace <- rstan::extract(fit_laplace, pars, permuted = F)

parm_label <- c(TeX("$\\c_mathrm{aux}"), 
                TeX("$\\rho$"), 
                TeX("$\\theta_1$"), 
                TeX("$\\theta_2$")) 

error_plot(samples_standard, samples_laplace, benchmark,
           pars = pars, one_chain = F, average_chain = T,
           log_parm = log_parm, f = quantile_90) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="lr")

# Next, examine the probability that log lambda is greater than
# a certain threshold
threshold <- 2.5
prob_threshold_log <- function(x) { sum(log(x) > threshold) / length(x) }
benchmark <- lapply(samples_benchmark, prob_threshold_log)
benchmark <- unlist(benchmark, use.names = F)

prob_threshold <- function(x) { sum(x > threshold) / length(x) }

error_plot(samples_standard, samples_laplace, benchmark,
           pars = pars, one_chain = F, average_chain = T,
           log_parm = log_parm, f = prob_threshold) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="lr")


#####################################################################
## Analysis for horseshoe regression
modelName <- "bernoulli_logit_glm_rhs"

fit_standard <- list()
fit_standard[[1]] <- readRDS(file = file.path("deliv", modelName, "cluster",
                             paste0(modelName, "_1_to_4Fit.Rsave")))
fit_standard[[2]] <- readRDS(file = file.path("deliv", modelName, "cluster",
                             paste0(modelName, "_5_to_6Fit.Rsave")))
fit_standard <- sflist2stanfit(fit_standard)

modelName <- "bernoulli_logit_glm_ela"
fit_laplace <- readRDS(file = file.path("deliv", modelName, "cluster",
                       paste0(modelName, "Fit.Rsave")))

pars <- c("tau", "lambda[1816]", "lambda[2586]", "caux")
summary_standard <- summary(fit_standard, pars = pars)[[1]]
summary_laplace <- summary(fit_laplace, pars = pars)[[1]]

# for the benchmark, read in the summary file generated by cmdstan
# Need to skip first and last row
summary_benchmark_full <- read.table("deliv_cmdstan/summary.dat",
                                     skip = 6, nrows = 18108)
summary_benchmark <- summary_benchmark_full[pars, ]

# plot Monte Carlo estimate with error (not for f, since we will
# examine these after applying an inverse logit transformation)
mcse_plot(summary_standard, summary_laplace, summary_benchmark, pars)

# plot error as a function of iterations
iter <- seq(from = 1, to = 2000, by = 1)
nChains <- 6

# WARNING -- make sure the order of parameters are consistent between
# samples and saved benchmarks!!
# pars <- c("tau", "caux", "lambda[2586]", "lambda[1816]")
pars <- c("tau", "caux", "lambda[1816]", "lambda[2586]")
samples_standard <- rstan::extract(fit_standard, pars = pars, permuted = FALSE)
samples_laplace <- rstan::extract(fit_laplace, pars = pars, permuted = FALSE)

# for testing purposes
logs_laplace <- log(samples_laplace)
logs_standard <- log(samples_standard)

log_parm = T

# Let's first start by analyzing the first 4 parameters (as we
# will put these on the log scale).
pars <- pars[1:4]
if (!log_parm) {
  benchmark <- summary_benchmark[, 1]
} else {
  nChains_benchmark <- 8
  mean_log <- array(NA, c(nChains_benchmark, length(pars)))
  fileName <- "bernoulli_logit_glm_rhs.logmean"
  for (i in 1:nChains_benchmark) {
    table <- read.csv(file = paste0("deliv_cmdstan/",
                      fileName, ".", i, ".csv"))
    # TODO -- find a slicker way to do the next operation.
    # mean_log[i, ] <- c(table[[2]])
    elements <- c(table[[2]])
    for (k in 2:length(pars)) elements <- c(elements,
                                         table[[k + 1]])
    
    mean_log[i, ] <- elements
  }

  benchmark <- colMeans(mean_log)
}

error_plot(samples_standard, samples_laplace, benchmark,
           pars, one_chain = F, average_chain = T,
           log_parm = log_parm, f = c(mean)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

# (For Neurips paper), let's examine the error for expectation value
# of caux and tau, and 90th quantile of lambda[1816] and lambda[2586].
# do the same analysis with quantiles.
modelName <- "bernoulli_logit_glm_rhs"
file <- file.path("deliv_cmdstan", paste0(modelName, ".90quantile.csv"))
benchmark_quantile <- read.csv(file = file)
benchmark_quantile <- unname(benchmark_quantile[, 2])

# Recall pars = (tau, caux, lambda[1816], lambda[2586])
benchmark <- c(benchmark[1], benchmark[2], benchmark_quantile[2], 
               benchmark_quantile[1])

# benchmark <- c(benchmark[1], benchmark_quantile[2],
#                benchmark_quantile[1], benchmark[4])

quantile_90 <- function(x) { quantile(x, probs = 0.9) }
f <- c(mean, mean, quantile_90, quantile_90)

include_warmup <- F
time_standard <- iteration_time(fit_standard, iter, include_warmup)
time_laplace <- iteration_time(fit_laplace, iter, include_warmup)

parm_label <- c(TeX("$E( \\log \ c_{aux})$"),
                TeX("$Q_{90}( \\log \\lambda_{1816})$"),
                TeX("$Q_{90}( \\log \\lambda_{2586})$"),
                TeX("$E( \\log \\tau)$"))

error_plot(samples_standard, samples_laplace, benchmark,
           pars, one_chain = F, average_chain = T,
           log_parm = T, f = f, plot_time = T,
           time_standard = time_standard,
           time_laplace = time_laplace) + #,
           # parm_label = parm_label) +
  theme(legend.position = "none", text = element_text(size = 25)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^floor(x)),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("time (s)") + ylab("squared error")


# Let's examine the probability for the first two patients
pars <- c("f[1]", "f[2]", "f[99]", "f[100]")
samples_standard <- rstan::extract(fit_standard, pars = pars, permuted = FALSE)
samples_standard <- inv.logit(samples_standard)
samples_laplace <- rstan::extract(fit_laplace, pars = pars, permuted = FALSE)
samples_laplace <- inv.logit(samples_laplace)

# read in expected value for probability
nChains_benchmark <- 8
mean_prob <- array(NA, c(length(pars), nChains_benchmark))
for (i in 1:nChains_benchmark) {
  chain_prob <- read.csv(paste0("deliv_cmdstan/", fileName, ".", i, ".csv"))
  mean_prob[, i] <- c(chain_prob$`f.1.`, chain_prob$`f.2.`,
                      chain_prob$`f.99.`, chain_prob$`f.100.`)
}
benchmark <- rowMeans(mean_prob)

benchmark
# 0.2024666 0.3154635 0.8985348 0.9599858

c(mean(samples_standard[, , 1]), mean(samples_standard[, , 2]),
  mean(samples_standard[, , 3]), mean(samples_standard[, , 4]))
# 0.2077414 0.3144079 0.8952827 0.9577220

c(mean(samples_laplace[, , 1]), mean(samples_laplace[, , 2]),
  mean(samples_laplace[, , 3]), mean(samples_laplace[, , 4]))
# 0.2588491 0.2721569 0.8141369 0.8930643

error_plot(samples_standard, samples_laplace, benchmark,
           pars = pars, one_chain = F, average_chain = T,
           log_parm = F, f = mean) +
  scale_y_continuous(trans = "log")

##
# Let's take a closer look at the lambdas
pars <- c("lambda[2586]", "lambda[1816]", "lambda[4960]",
          "lambda[4238]", "lambda[4843]", "lambda[3381]",
          "lambda[4647]")

# Read in benchmark for probability of exceeding threshold
# (the prob is computed for each chain, it suffices to
# average out the quantities).
modelName <- "bernoulli_logit_glm_rhs"
file <- file.path("deliv_cmdstan", paste0(modelName, ".prob"))

benchmark <- read.csv(file = paste0(file, "1.csv"))
benchmark <- as.numeric(unname(benchmark[, -1]))


# prob_threshold is defined above
samples_standard <- rstan::extract(fit_standard, pars = pars, permuted = FALSE)
samples_laplace <- rstan::extract(fit_laplace, pars = pars, permuted = FALSE)

prob_threshold_ <- function(x) { prob_threshold(x, threshold)}
error_plot(samples_standard, samples_laplace, benchmark,
           pars = pars, one_chain = F, average_chain = T,
           log_parm = log_parm, f = prob_threshold_) +
  scale_y_continuous(trans = "log")

# do the same analysis with quantiles.
file <- file.path("deliv_cmdstan", paste0(modelName, ".90quantile.csv"))
benchmark <- read.csv(file = file)
benchmark <- unname(benchmark[, 2])

error_plot(samples_standard, samples_laplace, benchmark,
           pars = pars, one_chain = F, average_chain = T,
           log_parm = log_parm, f = quantile_90) +
  scale_y_continuous(trans = "log")
