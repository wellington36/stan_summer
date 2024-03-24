functions {
  // p: lambda, p, n_i1, n_i2, ..., n_iT
  real logFunction(int k, array[] real param, array[] int ns) {
    vector[size(ns)] bins;
    for (i in 1:size(ns)) bins[i] = binomial_lpmf(ns[i] | k, param[2]);
    return poisson_lpmf(k | param[1]) + sum(bins);
  }
  
  #include "../n_mixture/infiniteSumToThreshold.stan"
}

data {
  int<lower = 1> R;
  int<lower = 1> T;
  array[R, T] int<lower = 0> n;
}

transformed data {
  array[R] int max_ns;
  
  for (i in 1:R) max_ns[i] = max(n[i,:]);
}

parameters {
  real<lower = 0, upper = 1> p;
  real<lower= 0> lambda;
}

transformed parameters {
  array[R] real components;
  array[2] real p_lambda;
  
  p_lambda[1] = p;
  p_lambda[2] = lambda;
  
  for (i in 1:R) components[i] = infiniteSumToThreshold(p_lambda, n[i, ],
    1e-14, 100000, max_ns[i])[1];
}

model {
  target += sum(components);
  
  p ~ beta(1, 1);
  lambda ~ gamma(0.01, 0.01);
}