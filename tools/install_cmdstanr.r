
# devtools::install_github("stan-dev/cmdstanr")
# remotes::install_github("jgabry/posterior")  # needs to updated

library(cmdstanr)
library(posterior)
set_cmdstan_path("~/Desktop/Code/laplace_approximation/cmdstan/")
cmdstan_path()
# install_cmdstan()

# Bernoulli model example
file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod <- cmdstan_model(file)
mod$print()
mod$exe_file()

data_list <- list(N = 10, y =c(0,1,0,0,0,0,0,0,0,1))
fit <- mod$sample(
  data = data_list,
  seed = 123,
  num_chains = 2,
  num_cores = 2
)

fit$summary()
fit$cmdstan_summary()
fit$cmdstan_diagnose()

stanfit <- rstan::read_stan_csv(fit$output_files())
print(stanfit)

# Further tests
file <- file.path(getwd(), "models3", 
                  "bernoulli_logit_glm_ela",
                  "bernoulli_logit_glm_ela.stan")
mod <- cmdstan_model(file)




