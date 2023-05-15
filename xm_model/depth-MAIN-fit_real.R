# -*- ess-style: RStudio; -*-

source("../shared/functions.R")

library("rstan")
library("bayesplot")


# CONFIG

RAD = 8.
MU_DX = 8./2
K = ceiling(RAD / MU_DX) - 1

USE = 0.2
ITER = 800
WARMUP = 400
CHAINS = 1
CORES = 1

MODELS_TO_RUN = c(2, 5, 6)


fit_config = list()

fit_config[[2]] = list(
  r = RAD,
  lb = -20.,
  ub = 0.,
  mu_dm_prior_mean = -0.5,
  mu_dm_prior_std = 0.5,
  sig_dm_prior_mean = 0.3,
  sig_dm_prior_std = 0.3
)

fit_config[[5]] = list(
  r = RAD,
  mu_dm_prior_mean = -0.5,
  mu_dm_prior_std = 0.5,
  sig_dm_prior_mean = 0.3,
  sig_dm_prior_std = 0.3,
  shape_dm_prior_mean = 0.0,
  shape_dm_prior_std = 4
)

fit_config[[6]] = list(
  r = RAD,
  min_depth = -20.,
  mu_dm_prior_mean = -0.5,
  mu_dm_prior_std = 0.5,
  sig_dm_prior_mean = 0.3,
  sig_dm_prior_std = 0.3
)

pars_of_interest = list()
pars_of_interest[[2]] = c("mu_dm", "sig_dm")
pars_of_interest[[5]] = c("mu_dm", "sig_dm", "shape_dm")
pars_of_interest[[6]] = c("mu_dm", "sig_dm")

SESSION_UUID = "REAL-DATA"

EPOCHS = c(0, 10, 18, 22, 28)


# Depth metadata
depths = -1 * (seq(0, 21, 1) * (16.2/21) + 0.5)
# depths
width = abs(mean(diff(depths)))
depth_ivals = c(depths[1] + 0.5 * width, depths - 0.5 * width)
# depth_ivals
pairwise_depths = 0.5 * (depths[seq(1, 21, 2)] + depths[seq(2, 22, 2)])
# pairwise_depths


# LOAD
df_all = load_rt_data("../data/acc-1-v3.csv") %>%
  filter(
    uptime > 0,
    electrode_pair > 1,
    days_since_start > 4,
    crop == "Corn"
  ) %>%
  mutate(
    old_count = count,
    count = 1 * (old_count > 0),
    device = as.factor(device_id),
    device_group = as.integer(device),
    time_group = days_since_start - min(days_since_start) + 1
  )

df_all$epoch = cut(df_all$days_since_start, breaks=EPOCHS)
df_all$epoch_group = as.integer(df_all$epoch)



# Sanity check
if (FALSE) {
  
  depth_bins_coarse = depth_ivals[seq(3, 23, by=4)]

  df_coarse = coarsen_rt_data(
    df_all,
    days_since_start = seq(0, 28, 4),
    depth_cm = depth_bins_coarse
  )

  df_time_elec = agg_rt_data(df_coarse, days_since_start, depth_cm) %>%
    group_by(days_since_start) %>%
    mutate(sum_rate = sum(rate)) %>%
    ungroup() %>%
    mutate(prop = rate / sum_rate)

  head(df_time_elec)

  df_time_elec %>%
    ggplot() +
    geom_tile(aes(days_since_start, depth_cm, fill=rate)) + 
    scale_fill_gradient2()

  df_time_elec %>%
    ggplot() +
    geom_line(aes(-as.integer(depth_cm), rate)) +
    facet_wrap(~ days_since_start)

  df_time_elec %>%
    ggplot() +
    geom_line(aes(-as.integer(depth_cm), prop)) +
    facet_wrap(~ days_since_start)

}
  

# Sample (for testing routines)
depths_df = df_all %>%
  filter(count > 0) %>%
  slice_sample(prop=USE, replace=FALSE)

# And get data for Stan
y_data = depths_df$depth_cm
group = depths_df$epoch_group
N = length(y_data)
M = length(levels(depths_df$epoch))

# Run models
samp_list = list()

for (model_num in MODELS_TO_RUN) {
  
  # Stan files
  stan_fit_file = sprintf("depth-m%02d-fit.stan", model_num)
  
  # Param
  dat_i = list(
    y = y_data,
    N = length(y_data),
    M = M,
    group = group,
    mu_dx = MU_DX,
    K = K
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
    control = list(max_treedepth=12),
    verbose = FALSE,
    model_name = sprintf("depth-fit-m%02d", model_num)
  )
  samp_list[[model_num]] = samp

  # Summary
  pars = pars_of_interest[[model_num]]
  print(summary(samp, pars))

  # Plot
  samp_array = extract(samp, pars, permute=FALSE)

  p_hist = mcmc_hist(samp_array) + ggtitle(sprintf("Posteriors M%02d", model_num))
  # p_hist

  p_hist_file = sprintf("%s-m%02d-real-p_hist.png", SESSION_UUID, model_num)
  ggsave(p_hist, file=file.path("..", "images", "xm_model", p_hist_file))

  p_pairs = mcmc_pairs(samp_array, grid_args = list(top=sprintf("Pairs M%02d", model_num)))
  # p_pairs

  p_pairs_file = sprintf("%s-m%02d-real-p_pairs.png", SESSION_UUID, model_num)
  ggsave(p_pairs, file=file.path("..", "images", "xm_model", p_pairs_file))

  # # Stan versions of plotting routines
  # pairs(samp, pars=pars)
  # traceplot(samp, pars=pars)

  # Save sample for later?

}
