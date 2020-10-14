
# Test code for ADVI
setwd("~/Code/laplace_approximation/Script/advi")
.libPaths("~/Rlib/")

library(cmdstanr)
library(rstan)
library(bayesplot)
library(dplyr)
library(tidyr)
library(boot)
library(latex2exp)
# set_cmdstan_path("~/Code/laplace_approximation/cmdStan/")
set_cmdstan_path("~/Rlib/cmdstan/")

scriptDir <- dirname(getwd())
modelDir <- file.path(scriptDir, "models3")
dataDir <- file.path(scriptDir, "data")

# Let's now try fitting one of the horseshoe GLM models.
modelName <- "skim_logit"  # "bernoulli_logit_glm_rhs"  # "disease_map"
dataFile <- "prostate_200.data.R"  # "prostate_full.data.R"  # "disease_data_100.r"

data <- read_rdump(file.path(dataDir, dataFile))

mod <- cmdstan_model(paste0(modelDir, "/", modelName, "/",
                            modelName, ".stan"))

# For bernoulli_logit_glm_rhs, use tol_rel_obj = 0.005
# For skim_logit, use tol_rel_obj = 0.007
# For disease_map, use tol_rel_obj = 0.005
fit_vb <- mod$variational(data = data,
                          algorithm = "meanfield",
                          iter = 1e6,
                          output_samples = 1.2e4,  # 2e3
                          eta = 0.1,
                          tol_rel_obj = 0.007,
                          seed = 123)

pars <- c("lambda", "tau", "caux")
if (modelName == "skim_logit") pars <- c(pars, "xi")

samples_raw_vb <- fit_vb$draws(pars)

## Compare samples from HMC to samples from VB.
## Use the code in cluster_analysis_v2.r
source(file.path(scriptDir, "tools", "analysisTools.r"))

if (modelName == "bernoulli_logit_glm_rhs") {
  nChains <- 6
  nIter <- 2000
  nIter_total <- nChains * nIter
  
  delivDir <- file.path(scriptDir, "deliv", modelName, "cluster")

  stanfit <- readRDS(file = file.path(delivDir, paste0(modelName,
                                        "_1_to_4", "Fit.Rsave")))
  stanfit2 <- readRDS(file = file.path(delivDir, paste0(modelName,
                                         "_5_to_6", "Fit.Rsave")))
  check_hmc_diagnostics(stanfit)  # 0
  check_hmc_diagnostics(stanfit2)  # 13

  samples <- rstan::extract(stanfit, pars = pars)
  samples2 <- rstan::extract(stanfit2, pars = pars)
  log_lambda <- log(rbind(samples$lambda, samples2$lambda))

  index <- select_lambda(log_lambda, quant = 0.9, n_select = 6)
  index

  n_cov <- 5966  # 200
  log_lambda_vb <- log(samples_raw_vb[, 1:n_cov])
  index_vb <- select_lambda(log_lambda_vb, quant = 0.9, n_select = 6)
  index_vb

  quant_select_plot2(log_lambda, log_lambda_vb, quant = 0.9,
                     threshold = 2.5, algo_name1 = "(full) HMC",
                     algo_name2 = "ADVI")

  tau <- c(samples$tau, samples2$tau)
  caux <- c(samples$caux, samples2$caux)
  samples_standard <-
    data.frame(log_lambda_1816 = log_lambda[, 1816],
               log_lambda_2586 = log_lambda[, 2586],
               tau = log(tau), caux = log(caux))

  samples_vb <-
    data.frame(log_lambda_1816 = log_lambda_vb[, 1816],
               log_lambda_2586 = log_lambda_vb[, 2586],
               tau = log(samples_raw_vb[, 5967]),
               caux = log(samples_raw_vb[, 5968]))
  names(samples_vb) <- c("log_lambda_1816", "log_lambda_2586",
                         "tau", "caux")

  n_parm = 4
  samples_all <- rbind(samples_standard, samples_vb)
  samples_all <- gather(samples_all)
  samples_all$method <- rep(rep(c("(full) HMC", "ADVI"), each = 12000),
                                n_parm)
  key_labels <- c(TeX("$\\log \ c_{aux}$"), TeX("$\\log \\lambda_{1816}$"),
                  TeX("$\\log \\lambda_{2586}$"), TeX("$\\log \\tau$")
                  )
  samples_all$key <- factor(samples_all$key, label = key_labels)

  pdf(file = paste("advi_comp","Plots%03d.pdf", sep = ""),
      width = 12, height = 3, onefile = F)
  
  # dimension: 3 x 12
  sample_comparison_plot(samples_all, x = 0.99, y = 0.95)

  dev.off()
}

if (modelName == "skim_logit") {
  nChains <- 4
  nIter <- 3000
  nIter_total <- nChains * nIter

  delivDir <- file.path(scriptDir, "deliv", modelName, "cluster")
  stanfit <- readRDS(file = file.path(delivDir, paste0(modelName, "Fit.Rsave")))
  check_hmc_diagnostics(stanfit)  # 0
  
  samples <- rstan::extract(stanfit, pars = pars)
  log_lambda <- log(rbind(samples$lambda))
  
  index <- select_lambda(log_lambda, quant = 0.9, n_select = 6)
  index
  
  n_cov <- 200
  log_lambda_vb <- log(samples_raw_vb[, 1:n_cov])
  index_vb <- select_lambda(log_lambda_vb, quant = 0.9, n_select = 6)
  index_vb  # 86  26 106  50 194 166

  samples_standard <-
    data.frame(log_lambda_2581 = log_lambda[, 81],
               log_lambda_2586 = log_lambda[, 86],
               tau = log(samples$tau),
               caux = log(samples$caux),
               xi = log(samples$xi))
  
  samples_vb <-
    data.frame(log_lambda_2581 = log_lambda_vb[, 81],
               log_lambda_2586 = log_lambda_vb[, 86],
               tau = log(samples_raw_vb[, 201]),
               caux = log(samples_raw_vb[, 202]),
               xi = log(samples_raw_vb[, 203]))
  names(samples_vb) <- c("log_lambda_2581", "log_lambda_2586",
                         "tau", "caux", "xi")
  
  n_parm <- 5
  samples_all <- rbind(samples_standard, samples_vb)
  samples_all <- gather(samples_all)
  samples_all$method <- rep(rep(c("(full) HMC", "ADVI"), each = 12000),
                            n_parm)
  key_labels <- c(TeX("$\\log \ c_{aux}$"), TeX("$\\log \\lambda_{2581}$"),
                  TeX("$\\log \\lambda_{2586}$"), TeX("$\\log \\tau$"),
                  TeX("$\\log \\chi$")
  )
  samples_all$key <- factor(samples_all$key, label = key_labels)
  
  pdf(file = paste("advi_comp","Plots%03d.pdf", sep = ""),
      width = 12, height = 3, onefile = F)
  sample_comparison_plot(samples_all, x = 0.79, y = 0.99)
  dev.off()
}

if (modelName == "disease_map") {
  nChains <- 4
  nIter <- 500
  nIter_total <- nChains * nIter
  
  delivDir <- file.path(scriptDir, "deliv", modelName)
  stanfit <- readRDS(file = file.path(delivDir, paste0(modelName, "Fit.Rsave")))
  check_hmc_diagnostics(stanfit)  # 0
  
  pars <- c("alpha", "rho", "theta[1]", "theta[2]")
  samples_raw_standard <- rstan::extract(stanfit, pars)
  samples_standard <-
    data.frame(alpha = samples_raw_standard$alpha,
               rho = samples_raw_standard$rho,
               theta_1 = samples_raw_standard$`theta[1]`,
               theta_2 = samples_raw_standard$`theta[2]`)
  
  samples_raw_vb <- fit_vb$draws(pars)
  samples_vb <- data.frame(alpha = samples_raw_vb[, 1],
                           rho = samples_raw_vb[, 2],
                           theta_1 = samples_raw_vb[, 3],
                           theta_2 = samples_raw_vb[, 4])
  names(samples_vb) <- c("alpha", "rho", "theta_1", "theta_2")
  
  n_parm <- 4
  samples_all <- rbind(samples_standard, samples_vb)
  samples_all <- gather(samples_all)
  samples_all$method <- rep(rep(c("(full) HMC", "ADVI"), each = 2000),
                            n_parm)
  key_labels <- c(TeX("$\\alpha$"), TeX("$\\rho$"),
                  TeX("$\\theta_1$"), TeX("$\\theta_2$"))
  samples_all$key <- factor(samples_all$key, label = key_labels)
  
  sample_comparison_plot(samples_all, x = 0.99, y = 0.95)
}
