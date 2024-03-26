functions{
  real logFunction(int k, array[] real p){
    return - 2*log(k+1) - (k+1)*log(p[1]);
  }
  real difference_check(real v1, real v2){
    if (v1 > v2)
      return(exp(log_diff_exp(v1, v2)));
    else
      return(exp(log_diff_exp(v2, v1)));
  }
  // logFunction must be defined beforehand
  #include infiniteSumToThreshold.stan
}

data {
  real<lower=0> p;
  real<lower=0> Epsilon;
  real TrueValue;
  int<lower=0>maxIter;
}

transformed data {
  array[1] real params = {p};
}

generated quantities {
  array[2] real estimatedSumToThreshold = infiniteSumToThreshold(params, Epsilon, maxIter, 0);
  real differenceSumToThreshold;
  real truth_1 = TrueValue;

  differenceSumToThreshold = difference_check(TrueValue, estimatedSumToThreshold[1]);
}