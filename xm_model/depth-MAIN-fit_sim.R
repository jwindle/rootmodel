# -*- ess-style: RStudio; -*-
# Jesse Windle, 2023

source(file.path("shared", "functions.R"))

library("rstan")
library("bayesplot")
library("yaml")


config = read_yaml(file.path("xm_model", "depth-MAIN-fit_sim.yaml"))

MODELS_TO_RUN = config$models_to_run
MODELS_TO_RUN

USE = config$use
ITER = config$mcmc$iter
WARMUP = config$mcmc$warmup
CHAINS = config$mcmc$chains
CORES = config$mcmc$cores


# SESSION_ID = "b7f0512c-f0d6-11ed-af6b-8c8590d39d9a"
# SESSION_ID = "e0236604-f0fc-11ed-ab75-8c8590d39d9a"
SESSION_ID = "TEST"

for (model_name in MODELS_TO_RUN) {

  fit_config = config$models[[model_name]]
  
  # Stan files
  stan_fit_file = file.path("xm_model", sprintf("depth-%s-fit.stan", model_name))
  load_file = file.path("data", sprintf("%s-%s-sim.RData", SESSION_ID, model_name))

  # Load data
  if (!file.exists(load_file)) {
    cat("FILE", load_file, "DOES NOT EXISTS\n")
    break
  }
  load(load_file) # Load samples, depths_list and param_df --- assuming these exist
  
  rand_idx = sample(1:nrow(param_df), size=1, replace=FALSE)

  param_df[rand_idx,]

  depths_ds = downsample_depths(depth_list[[rand_idx]], use=USE)
  
  y_data = depths_ds$depth
  N = length(y_data)

  hist(y_data)
  
  # Param
  dat_i = list(
    y = y_data,
    N = length(y_data),
    M = 1,
    group = rep(1, N),
    mu_dx = param_df$mu_dx[rand_idx],
    K = param_df$K[rand_idx]
  )
  dat = c(dat_i, fit_config$prior)
  # datenv = list2env(dat)
  
  # Simulate
  samp = stan(
    stan_fit_file,
    data = dat,
    chains = CHAINS,
    cores = CORES,
    iter = ITER,
    warmup = WARMUP,
    init_r = 1,
    control = list(max_treedepth=10),
    verbose = TRUE,
    model_name = sprintf("depth-fit-%s", model_name)
  )

  param_df[rand_idx,]
  
  # Summary
  pars = fit_config$pars
  
  print(summary(samp, pars))

  # Stan versions of plotting routines
  pairs(samp, pars=pars)

  traceplot(samp, pars=pars)

  # Plot
  param = as.numeric(param_df[rand_idx, pars])
  samp_array = extract(samp, pars, permute=FALSE)

  p_hist = mcmc_recover_hist(samp_array, param) + ggtitle(sprintf("Posteriors %s", model_name))
  p_hist

  p_hist_file = sprintf("%s-%s-p_hist.png", SESSION_ID, model_name)
  ggsave(p_hist, file=file.path("..", "images", "xm_model", p_hist_file))

  p_pairs = mcmc_pairs(samp_array, grid_args = list(top=sprintf("Pairs %s", model_name)))
  p_pairs

  p_pairs_file = sprintf("%s-%s-p_pairs.png", SESSION_ID, model_name)
  ggsave(p_pairs, file=file.path("..", "images", "xm_model", p_pairs_file))


}
