# Fit the ovarian and prostate cancer model by Aki,
# using cmdstanr (as opposed to the functions in cmdStanTools.R,
# as is done in the ovarian.r script).

rm(list = ls())
gc()
set.seed(1954)

# modelName <- "bernoulli_logit_glm_ela"
# modelName <- "skim_logit_ela"
modelName <- "skim_logit"

setwd("~/Desktop/Code/laplace_approximation/Script/")
.libPaths("~/Rlib")

# Setting for habenero cluster
# setwd("/rigel/home/ccm2172/laplace_approximation/Script")
# .libPaths("/rigel/home/ccm2172/Rlib")

scriptDir <- getwd()
modelDir <- file.path(scriptDir, "models3")
dataDir <- file.path(scriptDir, "data")
outDir <- file.path(scriptDir, "deliv", modelName)
delivDir <- file.path("deliv", modelName)

# dataFile <- "prostate_full.data.R"
dataFile <- "prostate_200.data.R"
# dataFile <- "ovarian_reduced.data.R"

library(rstan)
library(cmdstanr)
set_cmdstan_path("~/Desktop/Code/laplace_approximation/cmdStan/")
# set_cmdstan_path("/rigel/home/ccm2172/laplace_approximation/cmdStan")
library(parallel)
source(file.path("tools", "cmdStanTools.r"))
source(file.path("tools", "stanTools.r"))

## Take a look at the data
data <- read_rdump(file.path(dataDir, dataFile))
## faster with just 100 covariates 
## (for ovarian, we include the most relevant one)
# Comment out to fit the full data set.

if(FALSE) {
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

nChains <- 2
num_cores <- min(nChains, detectCores())

fit <- mod$sample(
  data = data, num_chains = nChains, num_cores = num_cores,
  num_warmup = 500, num_samples = 500, seed = 123,
  adapt_delta = 0.8, term_buffer = 10)

# CHECK -- do we need to create an rstan fit object?
stanfit <- read_stan_csv(fit$output_files())
saveRDS(stanfit, file = file.path(outDir, paste(modelName, "200_Fit.Rsave", sep = "")))
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

ptable <- parameterTable(stanfit, parms)
time <- sum(colSums(get_elapsed_time(stanfit)))
eff <- ptable[, 9] / time
eff_bulk <- ptable[, 21] / time
eff_tail <- ptable[, 22] / time
ptable <- cbind(ptable, eff, eff_bulk, eff_tail)

write.csv(ptable, file = file.path(delivDir, paste(modelName, 
          "ParameterTable.csv", sep = "")))

ptable[, c(1:3, 21, 22, 24, 25)]


ptable_standard <- ptable
ptable_standard[, c(1:3, 21, 22, 24, 25)]
