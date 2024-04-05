library(cmdstanr)
library(matrixStats)

Rcpp::cppFunction(code='
  double compute_lterm(long n, Rcpp::NumericVector p)
  {
  long double ans = - 2*log(n+1) - (n+1)* log(p[0]);
  return ans;
  }')

#### Facilitator functions ####
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())
my_dilog <- function(z) {
  if (abs(z) >= 1) {
    return(NA)
  }
  out = stats::integrate(function(y) log1p(-y)/y, lower = 0,
                         upper = z, subdivisions = 1E3,
                         rel.tol = 1E4 * .Machine$double.eps)$value
  return(-out)
}


system("rm stan/comparison_methods")
model <- cmdstan_model("stan/comparison_methods.stan",
                       include_paths = "stan")
####################

Eps <- .Machine$double.eps
M <- 5E5
a1 <- 2
L1 <- 1/a1
B1 <- 1/(a1-1)

lps.1 <- sapply(0:(1E6),
                function(j) compute_lterm(n = j, p = a1))
TV.1 <- matrixStats::logSumExp(sort(lps.1))
TTV.1 <- log(1/12 * (pi^2 - 6*log(2)^2))

cdata <- list(p = a1, Epsilon = Eps, maxIter = M,
              logL = -log(a1), n0 = 0, batchSize = 2)

comparison <- stanfit(
  model$sample(
    data = cdata,
    fixed_param = TRUE,
    iter_sampling = 1,
    adapt_engaged = FALSE,
    chains = 1
  )
)

( ext <- rstan::extract(comparison) )

naive.1 <- ext$estimatedSumToThreshold
adaptive.1 <- ext$estimatedErrorBoundingPairs
batches.1 <- ext$estimatedBatches

lapply(ext[grep("estimated", names(ext))],
       function(x) Rmpfr::mpfr(x, 1000))

####################
a2 <- 1.1
L2 <- 1/a2
B2 <- 1/(a2-1)

lps.2 <- sapply(0:(1E6),
                function(j) compute_lterm(n = j, p = a2))
TV.2 <- matrixStats::logSumExp(sort(lps.2))
TTV.2 <- log(my_dilog(1/a2))

cdata <- list(p = a2, Epsilon = Eps, maxIter = M,
              logL = -log(a2), n0 = 0, batchSize = round(B2/2))

comparison <- stanfit(
  model$sample(
    data = cdata,
    fixed_param = TRUE,
    iter_sampling = 1,
    adapt_engaged = FALSE,
    chains = 1
  )
)

( ext <- rstan::extract(comparison) )

naive.2 <- ext$estimatedSumToThreshold
adaptive.2 <- ext$estimatedErrorBoundingPairs
batches.2 <- ext$estimatedBatches

################

cdata <- list(p = a2, Epsilon = Eps, maxIter = M,
              logL = -log(a2), n0 = 0, batchSize = B2)

comparison <- stanfit(
  model$sample(
    data = cdata,
    fixed_param = TRUE,
    iter_sampling = 1,
    adapt_engaged = FALSE,
    chains = 1
  )
)

( ext <- rstan::extract(comparison) )

batches.2.correct <- ext$estimatedBatches

#################
source("aux/aux.r")

out1a <- tibble::tibble(
  method = c("Threshold", "BoundingPair", "Batches"),
  target_error = Eps,
  L = L1,
  a = a1,
  error = c(
    robust_difference(x = TV.1, y = naive.1[1]),
    robust_difference(x = TV.1, y = adaptive.1[1]),
    robust_difference(x = TV.1, y = batches.1[1])
  ),
  n_iter = c(
    naive.1[2],
    adaptive.1[2],
    batches.1[2]
  ),
  true_method = "Huge"
)

out1b <- tibble::tibble(
  method = c("Threshold", "BoundingPair", "Batches"),
  target_error = Eps,
  L = L1,
  a = a1,
  error = c(
    robust_difference(x = TTV.1, y = naive.1[1]),
    robust_difference(x = TTV.1, y = adaptive.1[1]),
    robust_difference(x = TTV.1, y = batches.1[1])
  ),
  n_iter = c(
    naive.1[2],
    adaptive.1[2],
    batches.1[2]
  ),
  true_method = "Closed-form"
)

out1 <- rbind(out1a, out1b)

out2a <- tibble::tibble(
  method = c("Threshold", "BoundingPair",
             "Batches_wrong", "Batches_right"),
  target_error = Eps,
  L = L2,
  a = a2,
  error = c(
    robust_difference(x = TV.2, y = naive.2[1]),
    robust_difference(x = TV.2, y = adaptive.2[1]),
    robust_difference(x = TV.2, y = batches.2[1]),
    robust_difference(x = TV.2, y = batches.2.correct[1])
  ),
  n_iter = c(
    naive.2[2],
    adaptive.2[2],
    batches.2[2],
    batches.2.correct[2]
  ),
  true_method = "Huge"
)

out2b <- tibble::tibble(
  method = c("Threshold", "BoundingPair",
             "Batches_wrong", "Batches_right"),
  target_error = Eps,
  L = L2,
  a = a2,
  error = c(
    robust_difference(x = TTV.2, y = naive.2[1]),
    robust_difference(x = TTV.2, y = adaptive.2[1]),
    robust_difference(x = TTV.2, y = batches.2[1]),
    robust_difference(x = TTV.2, y = batches.2.correct[1])
  ),
  n_iter = c(
    naive.2[2],
    adaptive.2[2],
    batches.2[2],
    batches.2.correct[2]
  ),
  true_method = "dilogarithm"
)

out2 <- rbind(out2a, out2b)

##################

out1 
robust_difference(x = TV.1, y = TTV.1)

out2
Eps * exp(-log(a2)-log1p(-L2))  ## This should be theoretical upper bound on the error
robust_difference(x = TV.2, y = TTV.2)

subset(rbind(out1, out2), true_method == "Huge")

##################

lapply(ext[grep("estimated", names(ext))],
       function(x) Rmpfr::mpfr(x, 1000))
