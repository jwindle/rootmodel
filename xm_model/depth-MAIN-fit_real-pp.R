# -*- ess-style: RStudio; -*-

# Generate posterior predictive

source("shared/functions.R")

library("rstan")
library("bayesplot")
library("yaml")
library("moments")


# Config
config = read_yaml("xm_model/depth-MAIN-fit_real-config.yaml")

SIMS_TO_RUN = config$simulations_to_run
SIMS_TO_RUN

# Load
load(file.path("cache", "detections.RData"))
EPOCHS = levels(df_time_elec$epoch)


# Simulate
samp_list = list()

for (sim_name in SIMS_TO_RUN) {

  sim_config = config$simulations[[sim_name]]
  
  sim_info = sim_config$sim_info
  sim_info$K = with(sim_info, get_K(r, mu_dx))

  model_name = sim_config$model
  stan_sim_file = sprintf("depth-%s-sim.stan", model_name)

  # This should load samp
  if (model_name %in% names(samp_list)) {
    samp = samp_list[[model_name]]
  } else {
    samp_file = sprintf("%s-fit.RData", model_name)
    load(file=file.path("cache", samp_file))
    samp_list[[model_name]] = samp
  }

  # samp = samp_list[[sim_name]]
  
  pars = sim_config$pars

  samp_array_2 = simplify2array(extract(samp, pars))

  base_dat = c(
    list(
      N = 10
    ),
    sim_info
  )

  sim_model = stan_model(file.path("xm_model", stan_sim_file), verbose=FALSE)

  simulated_paths = posterior_predictive_paths(
    sim_model,
    samp_array_2[seq(1, (ITER - WARMUP) * CHAINS, by=1),,],
    base_dat
  )

  simulated_depths_df = discretized_depth_distribution_2(
    simulated_paths,
    DEPTH_BINS,
    EPOCHS
  )

  if (config$write) {
    paths_file = sprintf("%s-paths.RData", sim_name)
    save(simulated_paths, simulated_depths_df, file=file.path("cache", paths_file))
  }

}

