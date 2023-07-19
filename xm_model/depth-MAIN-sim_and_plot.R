# -*- ess-style: RStudio; -*-

source("shared/functions.R")

library("uuid")
library("yaml")

config = read_yaml("xm_model/depth-MAIN-sim_and_plot.yaml")
param_config = config[["stan_data"]]

MODELS_TO_RUN = config$models_to_run
MODELS_TO_RUN



# Generate UUID to use for filenames
SESSION_ID = "TEST"
# SESSION_ID = UUIDgenerate(use.time=TRUE)


for (model_name in MODELS_TO_RUN) {

  # File names
  stan_sim_file = file.path("xm_model", sprintf("depth-%s-sim.stan", model_name))
  save_file = file.path("data", sprintf("%s-%s-sim.RData", SESSION_ID, model_name))

  # Param
  fixed_param = make_base_data(param_config[[model_name]])
  param_df = make_parameter_df(param_config[[model_name]])
  param_df$mu_dx = fixed_param$r / (param_df$K + 1)
  
  # Compile
  model = stan_model(stan_sim_file, model_name=stan_sim_file, verbose=TRUE)

  # Sample
  samp_list = simulate_over_parameter_df(model, param_df, fixed_param)
  
  # Extract summary and depths
  depth_list = extract_depths_list(samp_list)
  summary_df = summarize_depth_distributions(depth_list, param_df)
  depths_df = make_depths_df(depth_list, param_df)
  
  # PLOT - how?
  # Save
  save(depth_list, param_df, file=save_file)
  
}



## depths_ss = depths_df %>%
##   group_by(K, mu_yp, sig_yp) %>%
##   summarize(
##     m_depth = mean(depth),
##     s_depth = sd(depth)
##   )

## print(depths_ss, n=40)


## # Extract param names
## p_depths = depths_df %>%
##   filter(K == 3) %>%
##   ggplot() +
##   geom_histogram(aes(x=depth, y=..density.., fill=as.factor(mu_yp)), bins=10) +
##   facet_grid(col=vars(mu_yp), row=vars(sig_yp)) +
##   xlab("depth (cm)") +
##   ggtitle(TeX("Distribution of depths")) +
##   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
##   coord_flip()

## p_depths
