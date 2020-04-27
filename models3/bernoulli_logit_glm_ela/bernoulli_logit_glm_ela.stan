
functions {
  matrix K_functor (vector parm, matrix X,
                    real[] scale_icept, int[] delta_int) {
    int d = delta_int[1];
    vector[d + 1] k_diag;
    k_diag[1] = square(scale_icept[1]);
    k_diag[2:(d + 1)] = parm[1] * parm[2:(d + 1)];

    return X * diag_pre_multiply(k_diag, X');
  }
}

data {
  int<lower=0> n;				      // number of observations
  int<lower=0> d;             // number of predictors
  int<lower=0,upper=1> y[n];	// outputs
  matrix[n,d] x;				      // inputs
  real<lower=0> scale_icept;	// prior std for the intercept
  real<lower=0> scale_global;	// scale for the half-t prior for tau
  real<lower=1> nu_global;	  // degrees of freedom for the half-t priors for tau
  real<lower=1> nu_local;		  // degrees of freedom for the half-t priors for lambdas
                              // (nu_local = 1 corresponds to the horseshoe)
  real<lower=0> slab_scale;   // for the regularized horseshoe
  real<lower=0> slab_df;
}

transformed data {
  int delta_int[1] = {d};
  real delta[1] = {scale_icept};
  int n_samples[n] = rep_array(1, n);
  vector[n] theta0 = rep_vector(0, n);
  matrix[n, d + 1] X;  // design matrix with intercept
  X[, 1] = rep_vector(1.0, n);
  X[, 2:(d + 1)] = x;
}

parameters {
  real <lower=0> tau;         // global shrinkage parameter
  vector <lower=0>[d] lambda; // local shrinkage parameter
  real<lower=0> caux;
}

transformed parameters {
  vector[d + 1] parm;
  {
    vector[d] lambda_tilde;   // 'truncated' local shrinkage parameter
    real c = slab_scale * sqrt(caux); // slab scale
    lambda_tilde = sqrt( c^2 * square(lambda) ./ (c^2 + tau^2*square(lambda)));
    parm[1] = tau;
    parm[2:(d + 1)] = lambda_tilde;
  }
}

model {
  // half-t priors for lambdas and tau, and inverse-gamma for c^2
  lambda ~ student_t(nu_local, 0, 1);
  tau ~ student_t(nu_global, 0, scale_global*2);
  caux ~ inv_gamma(0.5*slab_df, 0.5*slab_df);

  target += laplace_marginal_bernoulli(y, n_samples, K_functor, parm,
                                       X, delta, delta_int, theta0);
}
generated quantities {
  vector[n] log_lik;
  vector[n] f = laplace_approx_bernoulli_rng(y, n_samples, K_functor, parm,
                                             X, delta, delta_int, theta0);

  for (i in 1:n) log_lik[i] = bernoulli_logit_lpmf(y[i] | f[i]);
}
