# -*- ess-style: RStudio; -*-

source("../shared/functions.R")

library("rstan")
library("bayesplot")
library("yaml")




MODELS_TO_RUN = c(2, 5, 6)

USE = 0.1
ITER = 3000
WARMUP = 1500
CHAINS = 2
CORES = 2


config = read_yaml("depth-MAIN-fit_sim.yaml")

fit_config = config[["prior_data"]]

pars_of_interest = list()
pars_of_interest[[2]] = c("mu_dm", "sig_dm")
pars_of_interest[[5]] = c("mu_dm", "sig_dm", "shape_dm")
pars_of_interest[[6]] = c("mu_dm", "sig_dm")


# SESSION_UUID = "b7f0512c-f0d6-11ed-af6b-8c8590d39d9a"
SESSION_UUID = "e0236604-f0fc-11ed-ab75-8c8590d39d9a"

for (model_num in MODELS_TO_RUN) {
  
  # Stan files
  stan_fit_file = sprintf("depth-m%02d-fit.stan", model_num)
  load_file = file.path("..", "data", sprintf("%s-m%02d-sim.RData", SESSION_UUID, model_num))

  # Load data
  if (!file.exists(load_file)) {
    cat("FILE", load_file, "DOES NOT EXISTS\n")
    break
  }
  load(load_file) # Load samples, depths_list and param_df --- assuming these exist
  
  rand_idx = sample(1:nrow(param_df), size=1, replace=FALSE)

  depths_ds = downsample_depths(depth_list[[rand_idx]], use=USE)
  
  y_data = depths_ds$depth
  N = length(y_data)
  
  # Param
  dat_i = list(
    y = y_data,
    N = length(y_data),
    M = 1,
    group = rep(1, N),
    mu_dx = param_df$mu_dx[rand_idx],
    K = param_df$K[rand_idx]
  )
  dat = c(dat_i, fit_config[[model_num]])
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
    verbose = FALSE,
    model_name = sprintf("depth-fit-m%02d", model_num)
  )

  # Summary
  pars = pars_of_interest[[model_num]]
  print(summary(samp, pars))

  # Plot
  # pars = c("mu_dm", "sig_dm")
  param = as.numeric(param_df[rand_idx, pars])
  samp_array = extract(samp, pars, permute=FALSE)

  p_hist = mcmc_recover_hist(samp_array, param) + ggtitle(sprintf("Posteriors M%02d", model_num))
  p_hist

  p_hist_file = sprintf("%s-m%02d-p_hist.png", SESSION_UUID, model_num)
  ggsave(p_hist, file=file.path("..", "images", "xm_model", p_hist_file))

  p_pairs = mcmc_pairs(samp_array, grid_args = list(top=sprintf("Pairs M%02d", model_num)))
  p_pairs

  p_pairs_file = sprintf("%s-m%02d-p_pairs.png", SESSION_UUID, model_num)
  ggsave(p_pairs, file=file.path("..", "images", "xm_model", p_pairs_file))

  # # Stan versions of plotting routines
  # pairs(samp, pars=pars)
  # traceplot(samp, pars=pars)

}
