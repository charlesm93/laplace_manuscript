rm(list = ls())
gc()

setwd("~/Desktop/Code/laplace_approximation/Script")
library(stringr)
library(rstan)

# Use function from http://stla.github.io/stlapblog/posts/Numextract.html
# to convert strings to the requisite numbers.
Numextract <- function(string){
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
}

# Read in data
synth <- read.csv("synth_data/testdata.csv")

nrow = 500  # 250
ncol = 3

if (FALSE) {
  synth_matrix <- matrix(nrow = nrow, ncol = ncol)
  for (i in 1:nrow) {
    dat_str <- toString(synth[i, 1])
    synth_matrix[i, ] <- as.numeric(Numextract(dat_str))
  }
  synth_data <- as.data.frame(synth_matrix)
  names(synth_data) <- c("x1", "x2", "y")

  # augment data to get a 500 dimensional latent variable
  x1 <- c(synth_data$x1, synth_data$x1 * runif(nrow, 1, 1.1))
  x2 <- c(synth_data$x2, synth_data$x2 * runif(nrow, 1, 1.1))
  y <- c(synth_data$y, synth_data$y)
}

if (TRUE) {
  synth_data <- synth
  x1 <- synth_data$x1
  x2 <- synth_data$x2
  y <- synth_data$y
  y[y == -1] <- 0  # express -1 as 0
}

sim_data <- data.frame(x1, x2, y)

# Write data for C++ unit tests.
output_directory <- "aki_synth_data"
write.table(t(x1), file = file.path(output_directory, paste0("x1.csv")),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
write.table(t(x2), file = file.path(output_directory, paste0("x2.csv")),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
write.table(t(y), file = file.path(output_directory, paste0("y.csv")),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")


###############################################################################
## Write data for Stan

make_data <- function(n_obs) {
  x1 <- sim_data$x1
  x2 <- sim_data$x2
  y <- sim_data$y

  if (n_obs <= 500) {
    data_stan <- sim_data[sample(1:500, n_obs), ]
  } else {
    for (i in (500 + 1):n_obs) {
      index <- sample(1:500, 1)
      x1 <- c(x1, sim_data$x1[index] * runif(1, 0.9, 1.1))
      x2 <- c(x2, sim_data$x2[index] * runif(1, 0.9, 1.1))
      y <- c(y, sim_data$y[index])
    }
    data_stan <- data.frame(x1, x2, y)
  }
  data_stan_ls <- with(data_stan, list(n_obs = n_obs,
                                       n_covariates = 2,
                                       y = y,
                                       x = cbind(x1, x2)))
  
  with(data_stan_ls, stan_rdump(ls(data_stan_ls),  # FIX ME
    paste0("data/data_gp/data_gp", n_obs, ".R")))
}

dimension <- c(10, 100, 250, 500, 750, 1000, 5000)
for (i in 1:length(dimension)) make_data(dimension[i])

# Write data for Stan (only use original data set, that is first
# 250 entries)
original_data <- sim_data[1:250, ]
data_org_stan <- with(original_data, list(n_obs = 250,
                                          n_covariates = 2,
                                          y = y,
                                          x = cbind(x1, x2)))
with(data_org_stan, stan_rdump(ls(data_org_stan), "data/data_gp/data_gp.R"))


###############################################################################
## Read data for desease map model in Finland
disease_data <- read.table("aki_disease_data/spatial1.txt")
names(disease_data) <- c("x1", "x2", "ye", "y")

# write data for C++ unit tests.
output_directory <- "aki_disease_data"
write.table(t(disease_data$x1), 
            file = file.path(output_directory, paste0("x1.csv")),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
write.table(t(disease_data$x2),
            file = file.path(output_directory, paste0("x2.csv")),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
write.table(t(disease_data$ye),
            file = file.path(output_directory, paste0("ye.csv")),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
write.table(t(disease_data$y),
            file = file.path(output_directory, paste0("y.csv")),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

# output data for Stan file
disease_stan <- with(disease_data, list(n_obs = 911,
                                        n_covariates = 2,
                                        y = y,
                                        ye = ye,
                                        x = cbind(x1, x2)))

with(disease_stan, stan_rdump(ls(disease_stan), "data/disease_data_full.r"))
