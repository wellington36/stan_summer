functions{
  #include infiniteSumToThreshold.stan
  #include infiniteErrorBoundingPairs.stan
  #include SumToCap.stan
}

// c-folding inifinite sum algorithm
// Requires definition of logFunction with two arguments:
// int k and array[] real parameters
array[] real infiniteBatches2(array[] real p, int batch_size, real epsilon, int maxIter, int n0) {
  logFunction
  
  infiniteErrorBoundingPairs(array[] real p, real epsilon, int maxIter, real logL, int n0)
}