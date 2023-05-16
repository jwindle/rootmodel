# -*- ess-style: RStudio; -*-

source("../shared/functions.R")

library("uuid")
library("yaml")


# Config
N = 1000 # N number of samples

RAD = 8.
MU_DX = c(8/2, 8/3, 8/4)
K = ceiling(RAD / MU_DX) - 1

MODELS_TO_RUN = c(2, 5, 6)


param_config = read_yaml("depth-MAIN-sim_and_plot.yaml")

## param_config = list()

## param_config[[2]] = list(
##   N = N,
##   r = RAD,
##   lb = -20.,
##   ub = 0.,
##   mu_dm = list(from=-1.0, to=-0.1, length.out=5),
##   sig_dm = list(from=0.1, to=0.4, length.out=3),
##   mu_dx = MU_DX,
##   K = K
## )

## param_config[[5]] = list(
##   N = N,
##   r = RAD,
##   mu_dm = list(from=-1.0, to=-0.1, length.out=5),
##   sig_dm = list(from=0.1, to=0.4, length.out=3),
##   shape_dm = list(from=-2.0, to=2.0, length.out=3),
##   mu_dx = MU_DX,
##   K = K
## )

## param_config[[6]] = list(
##   N = N,
##   r = RAD,
##   min_depth = -20.,
##   mu_dm = list(from=-1.0, to=-0.1, length.out=5),
##   sig_dm = list(from=0.1, to=0.4, length.out=3),
##   mu_dx = MU_DX,
##   K = K
## )


# Generate UUID to use for filenames
SESSION_UUID = "TEST"
# SESSION_UUID = UUIDgenerate(use.time=TRUE)


for (model_num in MODELS_TO_RUN) {
  # File names
  stan_sim_file = sprintf("depth-m%02d-sim.stan", model_num)
  save_file = file.path("..", "data", sprintf("%s-m%02d-sim.RData", SESSION_UUID, model_num))
  # Param
  fixed_param = make_base_data(param_config[[model_num]])
  param_df = make_parameter_df(param_config[[model_num]])
  # Simulate
  model = stan_model(stan_sim_file, model_name=stan_sim_file, verbose=TRUE)
  samp_list = simulate_over_parameter_df(model, param_df, fixed_param)
  # Extract summary and depths
  depth_list = extract_depths_list(samp_list)
  summary_df = summarize_depth_distributions(depth_list, param_df)
  depths_df = make_depths_df(depth_list, param_df)
  # PLOT - how?
  # Save
  save(depth_list, param_df, file=save_file)
}
