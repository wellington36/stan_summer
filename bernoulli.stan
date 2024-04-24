
data {
  int<lower=0> N; // Number of observations
  int<lower=0, upper=1> y[N]; // Binary outcome
}
parameters {
  real<lower=0, upper=1> theta; // Probability parameter
}
model {
  theta ~ beta(1, 1); // Prior
  y ~ bernoulli(theta); // Likelihood
}

