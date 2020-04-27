
# Analyze Stan fit objects generated on the cluster.

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
library(latex2exp)
source("tools/cmdStanTools.r")
source("tools/stanTools.r")
source("tools/analysisTools.r")

pars <- c("lambda", "tau", "caux", "f", "lp__")

nChains <- 6  # 4 + 2 chains over two runs.
nIter <- 2000  # 2000  # number of samples (so not including warmup)
nIter_total <- nChains * nIter

##########################################################################
## Read in results from Bernoulli logit model with adapt_delta = 0.99
modelName <- "bernoulli_logit_glm_rhs"

delivDir <- file.path("deliv", modelName, "cluster")
# stanfit <- readRDS(file = file.path(delivDir, paste0(modelName,
#                                     "Fit_long_warmup.Rsave")))
stanfit <- readRDS(file = file.path(delivDir, paste0(modelName,
                     "_1_to_4", "Fit.Rsave")))
stanfit2 <- readRDS(file = file.path(delivDir, paste0(modelName,
                      "_5_to_6", "Fit.Rsave")))
check_div(stanfit)
check_div(stanfit2)

# closer examination
sampler_params <- get_sampler_params(stanfit, inc_warmup = FALSE)
# sampler_params_chain1 <- sampler_params[[1]]
if (FALSE) {
divergence_by_chain <- 
  sapply(sampler_params, function(x) sum(x[, "divergent__"]))
divergence_by_chain
step_size_by_chain <- 
  sapply(sampler_params, function(x) mean(x[, "stepsize__"]))
step_size_by_chain
accept_by_chain <- 
  sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
accept_by_chain

mass_matrix <- get_mass_matrix(stanfit, nChains = 4)
str(mass_matrix)
}

# adapt_info <- get_adaptation_info(stanfit2)
# string <- "# Step size = 0.00654987\n# Diagonal elements of inverse mass matrix:\n#"
# start <- nchar(string) + 1
# end <- nchar(adapt_info[1]) - nchar("\n# ")
# mass_matrix_chain1 <- substring(adapt_info[1], start, end)
# mass_matrix_chain1 <- 
#   as.numeric(strsplit(mass_matrix_chain1, split =", ")[[1]])

samples <- rstan::extract(stanfit, pars = pars)
samples2 <- rstan::extract(stanfit2, pars = pars)
log_lambda <- log(rbind(samples$lambda, samples2$lambda))

index <- select_lambda(log_lambda, quant = 0.9, n_select = 6)
log_lambda_select <- log_lambda[, index]
names <- c(paste0("log_lambda[", index, "]"), "iteration", "chain")

posterior.sample <- construct_plot_data(log_lambda_select, nIter, nChains,
                                        names)
# trace_plot(posterior.sample)
density_hist(posterior.sample)
quant_select_plot(log_lambda, quant = 0.9, threshold = 2.5)

#extract other parameters
tau <- c(samples$tau, samples2$tau)
caux <- c(samples$caux, samples2$caux)
summary_table(log_lambda_select, tau, caux, index)

samples_standard <- data.frame(log_lambda_1816 = log_lambda[, 1816],
                               log_lambda_2586 = log_lambda[, 2586],
                               tau = tau,
                               caux = caux)


# inspect the predicted f
f <- rbind(samples$f, samples2$f)
p <- inv.logit(f)
p_expected <- colMeans(p)

##########################################################################
## Read in results from Bernoulli logit using ela
modelName2 <- "bernoulli_logit_glm_ela"
delivDir2 <- file.path("deliv", modelName2, "cluster")
stanfit_ela <- readRDS(file = file.path(delivDir2, paste0(modelName2, "Fit.Rsave")))
check_div(stanfit_ela)

pars <- c("lambda", "tau", "caux", "f", "lp__")
samples_ela <- rstan::extract(stanfit_ela, pars = pars)

# plot estimated probability
# inspect the predicted f
f_ela <- samples_ela$f
p_ela <- inv.logit(f_ela)
p_ela_expected <- colMeans(p_ela)


log_lambda_ela <- log(samples_ela$lambda)

# log_lambda <- log_lambda_ela
index <- select_lambda(log_lambda_ela, quant = 0.9, n_select = 6)
tau <- samples_ela$tau
caux <- samples_ela$caux
parm_select <- log_lambda_ela[, index]
names <- c(paste0("log_lambda[", index, "]"), "iteration", "chain")

# check samples for selected log_lambdas.
posterior.sample <- construct_plot_data(parm_select, nIter, nChains, names)
density_hist(posterior.sample)

# check samples for log tau.
posterior.sample.tau <- construct_plot_data(log(tau), nIter, nChains,
                                            c("log tau", "iteration", "chain"))
density_hist(posterior.sample.tau, bins = 30)
quant_select_plot(log_lambda, quant = 0.9, threshold = 2.5)

caux <- samples_ela$caux
summary_table(log_lambda_select, tau, caux, index)

#####################################################################
## Plots to save

pdf(file = file.path("deliv", "cluster_analysis",
                     paste(modelName,"Plots%03d.pdf",  sep = "")),
    width = 6, height = 6, onefile = F)

# do quantile plot using results from both models
quant_select_plot2(log_lambda, log_lambda_ela, quant = 0.9, threshold = 2.3)

samples_ela <- data.frame(log_lambda_1816 = log_lambda_ela[, 1816],
                          log_lambda_2586 = log_lambda_ela[, 2586],
                          tau = tau,
                          caux = caux) #,
                          # f = f)

samples_all <- rbind(samples_standard, samples_ela)
samples_all$log_tau <- log(samples_all$tau)
samples_all$log_caux <- log(samples_all$caux)
samples_all <- samples_all[, c(1, 2, 5, 6)]

samples_all <- gather(samples_all)
algorithm <- rep(rep(c("(full) HMC", "HMC + Laplace"), each = 12000), 4)
samples_all$method <- algorithm
key_labels <- c(TeX("$\\log \ c_{aux}$"), TeX("$\\log \\lambda_{1816}$"),
                TeX("$\\log \\lambda_{2586}$"), TeX("$\\log \\tau$"))
samples_all$key <- factor(samples_all$key, label = key_labels)

comp_plot <- ggplot(data = samples_all) +
  geom_histogram(aes(x = value, fill = method), alpha = 0.5, color = "black",
                 bins = 30, position = "identity") + theme_bw() +
  facet_wrap(~key, scale = "free", ncol = 1, labeller = "label_parsed") +
  theme(
    legend.position = c(.95, 0.17),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    text = element_text(size = 15)
  )
comp_plot

# next, examine ESS / s
pars = c("tau", "caux", "lambda[1816]", "lambda[2586]")
table_ela <- summary(stanfit_ela, pars = pars)[1]
ness_ela <- table_ela$summary[, 9]
time <- sum(get_elapsed_time(stanfit_ela)) / 6
eff_ela <- ness_ela / time

# table_standard1 <- summary(stanfit, pars = pars)[1]
ness_1 <- summary(stanfit, pars = pars)[1]$summary[, 9]
time1 <- sum(get_elapsed_time(stanfit))

table_standard2 <- summary(stanfit2, pars = pars)[1]
ness_2 <- table_standard2$summary[, 9]
time2 <- sum(get_elapsed_time(stanfit))

ness_standard <- ness_1 + ness_2
time <- (time1 + time2) / 6
eff_standard <- ness_standard / time 

data_eff <- data.frame(parameter = rep(pars, 2), eff = c(eff_standard, eff_ela),
                       method = rep(c("(full) HMC", "HMC + Laplace"),
                                       each = 4))

plot_eff <- ggplot(data = data_eff,
                   aes(x = parameter, y = eff, fill = method)) +
  geom_bar(stat = "identity", width = 0.3, alpha = 0.8, position = "dodge") + 
  # facet_wrap(~ parameter, scale = "free", nrow = 1) +
  theme_bw() + theme(text = element_text(size = 10)) + coord_flip() +
  ylab("ESS / s") + xlab(" ") +
  theme(
    legend.position = c(.95, 0.98),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  ) + theme(text = element_text(size = 18)) +
  scale_x_discrete(labels = c("tau" = parse(text = TeX("$\\tau$")),
                              "caux" = parse(text = TeX("$\\c_{aux}$")),
                              "lambda[1816]" = parse(text = TeX("$\\lambda_{1816}$")),
                              "lambda[2586]" = parse(text = TeX("$\\lambda_{2586}$"))))
plot_eff

plot_prob <- ggplot(data = data.frame(prob = p_expected, 
                                      prob_ela = p_ela_expected),
                    aes(x = prob, y = prob_ela)) +
  geom_point() + theme_bw() + geom_abline(intercept = 0, slope = 1, 
                                          color = "red", 
                                         linetype = "dashed", size = 1.0) +
  xlim(0, 1) + ylim(0, 1) + xlab("Probability (full HMC)") +
  ylab("Probability (HMC + Laplace)") +
  theme(text = element_text(size = 15))
plot_prob

dev.off()

###############################################################################
