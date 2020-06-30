##
# Functions to support error_analysis.r.

mcse_plot <- function(summary_standard,
                      summary_laplace,
                      summary_benchmark, pars) {
  mean <- c(summary_standard[, 1], summary_benchmark[, 1], 
            summary_laplace[, 1])
  mcse <- c(summary_standard[, 2], summary_benchmark[, 2],
            summary_laplace[, 2])

  parm <- rep(pars, 3)
  method <- rep(c("(full) HMC", "benchmark", "HMC + Laplace"),
                each = length(pars))
  
  # method <- rep(c("(full) HMC", "HMC + Laplace", "benchmark"))
  ggplot(data = data.frame(mean, mcse, parm, method),
         aes(x = method, y = mean, ymin = mean - mcse, ymax = mean + mcse,
             color = method)) + # geom_point(size = 3) +
    facet_wrap(~parm, scale = "free") + geom_pointrange(shape = 20) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}

# Return the squared error as a function of iterations.
# By default, the error is computed for each chain and then
# averaged.
# Due to dimensionality, we need to treat the case where we
# have one chain seperately.
# f specifies the functions use for the estimate (for the classic
# Monte Carlo estimate, f is the mean. Other applications include
# quantiles and probability of exceeding a threshold).
# f must specified as a vector of functions, to accomodate the
# case where we apply a different function for each parameter.
squared_error <- function(samples, benchmark, average_chain = T,
                          one_chain = F, f = c(mean)) {
  n_iter <- dim(samples)[1]
  n_chain <- dim(samples)[2]
  n_parm <- dim(samples)[3]
  
  if (length(f) == 1) {
    f_new <- f
    for (i in 2:n_parm) f_new <- c(f_new, f)
    f <- f_new
  }

  strat_err <- array(NA, c(n_iter, n_parm))

  for (i in 1:n_iter) {
    for (j in 1:n_parm) {
      if(!one_chain) {
        if (average_chain) {
          chain_err <- rep(NA, n_chain)
          for (k in 1:n_chain) {
            chain_err[k] <- (f[[j]](samples[1:iter[i], k, j]) - benchmark[j])^2
          }
          strat_err[i, j] <- mean(chain_err)
        } else {  # do not average chain
          strat_err[i, j] <-
            (f[[j]](samples[1:iter[i], , j]) - benchmark[j])^2
        }
      } else {  # only one chain
        strat_err[i, j] <- (f(samples[1:iter[i], j]) - benchmark[j])^2
      }
    }
  }
  strat_err
}

iteration_time <- function(fit, iter, include_warmup = F) {
  elapsed_time <- get_elapsed_time(fit)
  warmup <- mean(elapsed_time[, 1])
  sampling <- mean(elapsed_time[, 2])
  
  time <- (iter / max(iter)) * sampling
  if (include_warmup) time <- time + warmup
  time
}

error_plot <- function(samples_standard, samples_laplace, benchmark,
                       pars,
                       plot_time = FALSE, time_standard = NULL,
                       time_laplace = NULL, one_chain = F,
                       average_chain = T, log_parm = F,
                       f = c(mean), parm_label = NULL) {
  n_iter <- dim(samples_standard)[1]
  iter <- 1:n_iter
  n_chain <- dim(samples_standard)[2]
  n_pars <- length(pars)

  if (log_parm) {
    samples_standard <- log(samples_standard)
    samples_laplace <- log(samples_laplace)
    print("Make sure the benchmark averages the log!")
  }
  
  err_standard <- c(squared_error(samples_standard, benchmark, average_chain,
                                  one_chain, f = f))
  err_laplace <- c(squared_error(samples_laplace, benchmark, average_chain,
                                 one_chain, f = f))
  err <- c(err_standard, err_laplace)

  method <- rep(c("(full) HMC", "HMC + Laplace"),
                each = n_iter * n_pars)
  # TODO - if we don't average chains, iter should be iter * n_chain
  sample_size <- rep(iter, 2 * n_pars)
  
  parm <- rep(rep(pars, each = n_iter), 2)
  if (!is.null(parm_label)) parm <- factor(parm, label = parm_label)

  plot.data <- data.frame(err, parm, method, sample_size)
  
  if (!plot_time) {
    plot <- ggplot(data = plot.data, aes(x = sample_size, 
                                         y = err, color = method))
  }
  
  if (plot_time) {
    plot.data$time <- c(rep(time_standard, length(pars)),
                        rep(time_laplace, length(pars)))
    plot <- ggplot(data = plot.data, aes(x = time, y = err, color = method))
  }
  
  plot <- plot + geom_line() + theme_bw() +
    facet_wrap(~parm, scale = "free", labeller = "label_parsed",
               nrow = 1)

  plot
}
