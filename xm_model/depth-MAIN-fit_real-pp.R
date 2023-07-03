# -*- ess-style: RStudio; -*-

# Generate posterior predictive

source("shared/functions.R")

library("rstan")
library("bayesplot")
library("yaml")


# Config
config = read_yaml("xm_model/depth-MAIN-fit_real-log-config.yaml")

SIMS_TO_RUN = config$simulations_to_run
# SIMS_TO_RUN = c("m09_2k", "m09_3k", "m10_1k", "m10_2k", "m10_3k")
SIMS_TO_RUN

# Load
load(file.path("cache", "detections.RData"))
EPOCHS = levels(df_time_elec$epoch)


# Simulate
samp_list = list()

for (sim_name in SIMS_TO_RUN) {

  sim_config = config$simulations[[sim_name]]
  model_config = config$models[[sim_config$model]]

  base_dat = c(
    list(
      N = 10
    ),
    model_config$constants,
    sim_config$sim_info
  )

  base_dat$K = with(base_dat, get_K(r, mu_dx))

  model_name = sim_config$model
  
  stan_sim_file = sprintf("depth-%s-sim.stan", model_config$stan_file_id)

  ## # This should load samp
  ## if (model_name %in% names(samp_list)) {
  ##   samp = samp_list[[model_name]]
  ## } else {
  ##   samp_list[[model_name]] = samp
  ## }
  # samp = samp_list[[sim_name]]
  
  samp_file = sprintf("%s-fit.RData", model_name)
  load(file=file.path("cache", samp_file))
  
  pars = sim_config$pars

  samp_array_2 = simplify2array(extract(samp, pars))

  sim_model = stan_model(file.path("xm_model", stan_sim_file), verbose=FALSE)

  ## temp_test = sampling(
  ##   sim_model,
  ##   data=c(base_dat, as.list(samp_array_2[1,1,])),
  ##   algorithm="Fixed_param",
  ##   iter=1,
  ##   chains=1
  ## )

  samp_array_sub = samp_array_2[with(model_config$mcmc, seq(1, (iter - warmup) * chains, by=10)),,]
  
  sim_paths = posterior_predictive_paths(
    sim_model,
    samp_array_sub,
    base_dat
  )

  sim_paths_refined = refine_paths(sim_paths, base_dat$mu_dx, 10, model_config$transform)

  sim_depths_df = discretized_depth_distribution_3(
    sim_paths_refined,
    DEPTH_BINS,
    EPOCHS
  )

  if (config$write) {
    paths_file = sprintf("%s-paths.RData", sim_name)
    save(
      sim_paths,
      sim_paths_refined,
      sim_depths_df,
      file=file.path("cache", paths_file)
    )
  }

}

