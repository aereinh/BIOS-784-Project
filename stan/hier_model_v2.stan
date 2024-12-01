// hier_model.stan
// Hierarchical model
data {
  int<lower=1> N;             // Number of observations
  int<lower=1> I;             // Number of regions
  int<lower=1, upper=I> region[N]; // Region indicator for each observation
  vector[N] x;                // Predictor (univariate)
  int<lower=0, upper=1> y[N]; // Binary response
}

parameters {
  real alpha0;                // Fixed intercept
  real alpha1;                // Fixed slope
  vector[I] mu0;              // Region-level random intercepts
  vector[I] mu1;              // Region-level random slopes
  real<lower=0> sigma_mu0;    // Standard deviation of random intercepts
  real<lower=0> sigma_mu1;    // Standard deviation of random slopes
}

model {
  // Priors
  alpha0 ~ normal(0, 10);
  alpha1 ~ normal(0, 10);
  mu0 ~ normal(0, sigma_mu0);
  mu1 ~ normal(0, sigma_mu1);
  sigma_mu0 ~ cauchy(0, 2.5); // Half-Cauchy prior for sigma_mu0
  sigma_mu1 ~ cauchy(0, 2.5); // Half-Cauchy prior for sigma_mu1

  // Likelihood
  for (n in 1:N) {
    y[n] ~ bernoulli_logit(alpha0 + alpha1 * x[n] + mu0[region[n]] + mu1[region[n]] * x[n]);
  }
}