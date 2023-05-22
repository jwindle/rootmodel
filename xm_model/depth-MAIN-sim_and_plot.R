# -*- ess-style: RStudio; -*-

source("../shared/functions.R")

library("uuid")
library("yaml")

# Config
MODELS_TO_RUN = c(2, 5, 6)

config = read_yaml("depth-MAIN-sim_and_plot.yaml")

param_config = config[["stan_data"]]



# Generate UUID to use for filenames
SESSION_ID = "TEST"
# SESSION_ID = UUIDgenerate(use.time=TRUE)


for (model_num in MODELS_TO_RUN) {
  # File names
  stan_sim_file = sprintf("depth-m%02d-sim.stan", model_num)
  save_file = file.path("..", "data", sprintf("%s-m%02d-sim.RData", SESSION_ID, model_num))
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
