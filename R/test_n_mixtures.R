library(cmdstanr)

system("rm n_mixture/n_mixtures")
model <- cmdstan_model("n_mixture/n_mixtures.stan",
                       include_paths = "stan")
