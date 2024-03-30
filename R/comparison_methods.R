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
B1<- 1/(a1-1)

lps.1 <- sapply(0:(1E6),
                function(j) compute_lterm(n = j, p = a1))
TV.1 <- matrixStats::logSumExp(sort(lps.1))
TTV.1 <- log(1/12 * (pi^2 - 6*log(2)^2))

cdata <- list(p = a1, Epsilon = Eps, maxIter = M,
              TrueValue = TV.1)

comparison <- stanfit(
  model$sample(
    data = cdata,
    chains = 1,
    iter_warmup = 1,
    iter_sampling = 1,
    fixed_param = TRUE
  )
)

( ext <- rstan::extract(comparison) )

Rmpfr::mpfr(ext$estimatedSumToThreshold[1], 1000)
Rmpfr::mpfr(cdata$TrueValue, 1000)
Rmpfr::mpfr(ext$truth_1, 1000)
Rmpfr::mpfr(ext$difference, 1000)

