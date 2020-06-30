
rm(list = ls())
gc()
set.seed(1954)

# modelName <- "bernoulli_logit_glm_ela"
# modelName <- "skim_logit_ela"
modelName <- "skim_logit"

# Adjust paths to your settings
setwd("~/Desktop/Code/laplace_approximation/Script/")
.libPaths("~/Rlib")

scriptDir <- getwd()
modelDir <- file.path(scriptDir, "models3")
dataDir <- file.path(scriptDir, "data")
outDir <- file.path(scriptDir, "deliv", modelName)
delivDir <- file.path("deliv", modelName)

dataFile <- "prostate_full.data.R"
# dataFile <- "prostate_200.data.R"

library(rstan)
library(cmdstanr)
set_cmdstan_path("~/Desktop/Code/laplace_approximation/cmdStan/")
library(parallel)
source(file.path("tools", "cmdStanTools.r"))
source(file.path("tools", "stanTools.r"))

## Take a look at the data
data <- read_rdump(file.path(dataDir, dataFile))


if(FALSE) {
  # Use this sequence to examine a subset of the data
  data$x <- data$x[, (data$d - 99):data$d]

  # If we only retain a subset of x, we need to adjust tau and d.
  data$p <- ncol(data$x)
  p0 <- 5
  sigma <- with(data, sqrt(1 / mean(y) / (1 - mean(y))))
  data$scale_global <- with(data, p0 / (p - p0) * sigma / sqrt(n))
  data$d <- ncol(data$x)

  # Update the Rdump file
  dataFile <- "prostate.data.R"
  with(data, stan_rdump(ls(data), file.path(dataDir, dataFile)))
}

# compile model using cmdstanr
if (TRUE) {
  file <- file.path(modelDir, modelName, paste0(modelName, ".stan"))
  mod <- cmdstan_model(file)
}

nChains <- 6
num_cores <- min(nChains, detectCores())

fit <- mod$sample(
  data = data, num_chains = nChains, num_cores = num_cores,
  num_warmup = 1000, num_samples = 12000, seed = 123,
  adapt_delta = 0.8, term_buffer = 10)

stanfit <- read_stan_csv(fit$output_files())
saveRDS(stanfit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
check_div(stanfit)

parms <- c("log_lambda[1]", "log_lambda[2]", "log_lambda[3]",
           "tau", "caux")

pdf(file = file.path(delivDir, paste(modelName,"Plots%03d.pdf",  sep = "")),
    width = 6, height = 6, onefile = F)

mcmcHistory(stanfit, parms)
mcmcDensity(stanfit, parms, byChain = TRUE)
mcmcDensity(stanfit, parms)
# pairs(fit, parms)

dev.off()
