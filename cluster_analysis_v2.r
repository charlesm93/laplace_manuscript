
# Analyze Stan fit object generated on the cluster
# The here presented code is

rm(list = ls())
gc()

setwd("~/Desktop/Code/laplace_approximation/Script/")
.libPaths("~/Rlib")
scriptDir <- getwd()

library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
library(boot)
source("tools/cmdStanTools.r")
source("tools/stanTools.r")
source("tools/analysisTools.r")

# pars <- c("lambda", "tau", "caux", "f", "lp__")
pars <- c("lambda", "tau", "caux", "xi")

nChains <- 6
nIter <- 2000
nIter_total <- nChains * nIter

modelName <- c("skim_logit", "skim_logit_ela")
delivDir <- file.path("deliv", modelName, "cluster")

stanfit <- readRDS(file = file.path(delivDir[1],
                                    paste0(modelName[1], "Fit.Rsave")))

check_div(stanfit)
samples_standard <- rstan::extract(stanfit, pars = pars)
index_standard <- select_lambda(samples_standard$lambda,
                                quant = 0.9, n_select = 6)

stanfit2 <- readRDS(file = file.path(delivDir[2],
                                     paste0(modelName[2], "Fit.Rsave")))
check_div(stanfit2)
samples_ela <- rstan::extract(stanfit2, pars = pars)
index_ela <- select_lambda(samples_ela$lambda,
                           quant = 0.9, n_select = 6)

pdf(file = file.path("deliv", "cluster_analysis",
                     paste(modelName,"Plots%03d.pdf",  sep = "")),
    width = 6, height = 6, onefile = F)
quant_select_plot2(log(samples_standard$lambda),
                   log(samples_ela$lambda), quant = 0.9,
                   threshold = 2.3, alpha = 0.5, x = 0.98, y = 0.9,
                   index_offset = 2500)

# Let's look at the first two covariates in index_standard, 81 and 86.
# We add 2000, given we've taken a subsample.
# NOTE: we can also look at xi.
samples_select <- function(sample) { 
  data.frame(
  log_lambda_2581 = log(sample$lambda[, 81]),
  log_lambda_2586 = log(sample$lambda[, 86]),
  tau = log(sample$tau),
  caux = log(sample$caux),
  xi = log(sample$xi))
}

posterior_standard <- samples_select(samples_standard)
n_parm <- ncol(posterior_standard)
posterior_ela <- samples_select(samples_ela)
posterior_all <- rbind(posterior_standard, posterior_ela)
posterior_all <- gather(posterior_all)
posterior_all$method <- 
  rep(rep(c("(full) HMC", "HMC + Laplace"), each = 12000), n_parm)
key_labels <- c(TeX("$\\log \ c_{aux}$"), TeX("$\\log \\lambda_{1816}$"),
                TeX("$\\log \\lambda_{2586}$"), TeX("$\\log \\tau$"),
                TeX("$\\log \\chi"))
posterior_all$key <- factor(posterior_all$key, label = key_labels)

sample_comparison_plot(posterior_all)
  
pars = c("tau", "caux", "xi", "lambda[81]", "lambda[86]")
ness_standard <- summary(stanfit, pars = pars)[1]$summary[, 9]
time_standard <- sum(get_elapsed_time(stanfit)) / 6
eff_standard <- ness_standard / time_standard


ness_ela <- summary(stanfit2, pars = pars)[1]$summary[, 9]
time_ela <- sum(get_elapsed_time(stanfit2)) / 6
eff_ela <- ness_ela / time_ela

data_eff <- data.frame(parameter = rep(pars, 2), eff = c(eff_standard, eff_ela),
                       method = rep(c("(full) HMC", "HMC + Laplace"),
                                    each = n_parm))

eff_comparison_plot(data_eff, y = 0.8) +
  scale_x_discrete(labels = c("tau" = parse(text = TeX("$\\tau$")),
                              "caux" = parse(text = TeX("$\\c_{aux}$")),
                              "xi" = parse(text = TeX("$\\chi$")),
                              "lambda[86]" = parse(text = TeX("$\\lambda_{2586}$")),
                              "lambda[81]" = parse(text = TeX("$\\lambda_{2581}$")))) +
  theme(text = element_text(size = 18))

dev.off()

