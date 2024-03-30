functions{
  real logFunction(int k, array[] real p){
    return - 2*log(k+1) - (k+1)*log(p[1]);
  }
  // logFunction must be defined beforehand
  #include infiniteSumToThreshold.stan
  #include infiniteErrorBoundingPairs.stan
  #include infiniteBatches.stan
}

data {
  real<lower=0> p;
  real<lower=0> Epsilon;
  real logL;
  int<lower=0> maxIter;
  int<lower=0> batchSize;
  int<lower=0> n0;
}

transformed data {
  array[1] real params = {p};
}

generated quantities {
  array[2] real estimatedSumToThreshold = infiniteSumToThreshold(params, Epsilon, maxIter, n0);
  array[2] real estimatedErrorBoundingPairs = infiniteErrorBoundingPairs(params, Epsilon, maxIter, logL, n0);
  array[2] real estimatedBatches = infiniteBatches(params, batchSize, Epsilon, maxIter, n0);
}