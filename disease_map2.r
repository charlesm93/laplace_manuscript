# Reproduce disease map model, as developed by Aki and colleagues

rm(list = ls())
gc()
set.seed(1954)

# modelName <- "disease_map"
modelName <- "disease_map_ela"
# modelName <- "skim_logit_ela"

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

dataFile <- "disease_data_full.r"
# dataFile <- "disease_data_100.r"

library(rstan)
library(cmdstanr)
set_cmdstan_path("~/Desktop/Code/laplace_approximation/cmdStan/")
# set_cmdstan_path("/rigel/home/ccm2172/laplace_approximation/cmdStan")
library(parallel)
# library(rootSolve)
# library(invgamma)
source(file.path("tools", "cmdStanTools.r"))
source(file.path("tools", "stanTools.r"))

data <- read_rdump(file.path(dataDir, dataFile))

reduced <- TRUE
n_sample <- 250  # 911
if (reduced) {
  # select_index <- sample(1:911, n_sample)
  bin <- floor(911 / n_sample)
  select_index <- rep(NA, n_sample)
  for (i in 1:n_sample) select_index[i] <-
      sample(((i - 1) * bin + 1):(i * bin), 1)
  # select_index <- seq(from = 1, to = 910, by = 5)
  data$x <- data$x[select_index, ]
  data$y <- data$y[select_index]
  data$ye <- data$ye[select_index]
  data$n_obs <- n_sample
}

##########################################################################
# Construct prior (see Michael's dsicussion on inverse-gamma). Use
# range of the data to find lower and upper bounds for rho.
if (FALSE) {
  fit <- stan(file = "models2/inv_gamma_prior.stan", iter = 1, warmup = 0, chains = 1,
              seed = 5838298, algorithm = "Fixed_param")
}
# a = 2.42393
# b = 14.8171
data$rho_alpha_prior <- 2.42393
data$rho_beta_prior <- 14.8171
with(data, stan_rdump(ls(data), "data/disease_data.r"))

###########################################################################
if (TRUE) {
  file <- file.path(modelDir, modelName, paste0(modelName, ".stan"))
  mod <- cmdstan_model(file)
}

nChains <- 4
num_cores <- min(nChains, detectCores())

fit <- mod$sample(
  data = data, chains = nChains, parallel_chains = num_cores,
  num_warmup = 500, num_samples = 500, seed = 123,
  adapt_delta = 0.8)
stanfit <- read_stan_csv(fit$output_files())
dir.create(outDir)
saveRDS(stanfit, file = file.path(outDir, paste(modelName, n_sample,
                                         "_Fit.Rsave", sep = "")))


# parms <- c("alpha", "rho", "theta_pred[1]", "theta_pred[2]", "lp__")
parms <- c("alpha", "rho", "theta[1]", "theta[2]", "lp__")

dir.create(delivDir)
pdf(file = file.path(delivDir, paste(modelName,"Plots%03d.pdf", 
                                     sep = "")),
    width = 6, height = 6, onefile = F)

mcmcHistory(stanfit, parms)
mcmcDensity(stanfit, parms, byChain = TRUE)
mcmcDensity(stanfit, parms)

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

