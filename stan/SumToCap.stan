// Sum to a fixed cap
// Requires definition of logFunction with two arguments:
// int k and array[] real parameters
array[] real SumToCap(array[] real p, int nmax, int n0) {
  vector[nmax + 1] storeVal;
  int n = nmax + 1;
  for (i in 1:n){
    storeVal[i] = logFunction(n0 + i - 1, p);
  }
  return {log_sum_exp(sort_asc(storeVal[:n])), 1. * n};
}
