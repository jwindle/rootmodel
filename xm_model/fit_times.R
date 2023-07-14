# -*- ess-style: RStudio -*- #
# Jesse Windle

library("ggplot2")
library("gridExtra")

library("rstan")

source("shared/functions.R")


# Config
CROP = "Wheat"


# Load
df_all = load_rt_data("data/acc-1-v3.csv") %>%
    filter(uptime > 0, electrode_pair > 1, days_since_start > 4, crop == CROP)

df_all$old_count = df_all$count
df_all$count = 1 * (df_all$old_count > 0)

time_grid = unique(df_all$days_since_start)
df_all$time_group =  match(df_all$days_since_start, time_grid)


# Simplify and plot data

# By device and time
df_dev_time = agg_rt_data(df_all, device_id, time_group)

# 12 paddles 20 electrodes 12 obs / hour --- so this is the number of
# 12 x 20 (or 12 x 22) matrices we observe.
df_dev_time = df_dev_time %>%
  mutate(
    uptime_hours = uptime / (12 * 20 * 12),
    uptime_days = uptime_hours / 24,
    n_trials = uptime_hours * 2,
    rate = count / uptime_days
  )

head(df_dev_time)

# Just time
df_time = agg_rt_data(df_all, time_group) %>%
  mutate(
    uptime_hours = uptime / (12 * 20 * 12),
    uptime_days = uptime_hours / 24,
    n_trials = uptime_hours * 2,
    rate = count / uptime_days
  )
df_time$time_grid = time_grid[df_time$time_group]


time_plot <- df_time %>%
  ggplot() +
  geom_point(aes(time_grid, rate))

time_plot


# Model
fit_times_mdl = stan_model(file.path("xm_model", "fit_times.stan"), model_name="fit_times", verbose=TRUE)

# Data
max_trials = 48
dat = list(
  n_obs = nrow(df_dev_time),
  n_times = length(time_grid),
  max_trials = max_trials,
  time_grid = time_grid,
  time_group = df_dev_time$time_group,
  n_trials = df_dev_time$n_trials,
  y = df_dev_time$count,
  per_low = 7,
  per_high = 21
)

# Fit - had to up adapt_delta to 0.95 to git rid of divergent
# transitions.
fit_times_out = sampling(
  fit_times_mdl,
  data=dat,
  iter=7000,
  warmup=4000,
  chains=3,
  cores=3,
  control = list(
    adapt_delta = 0.99,
    adapt_t0 = 10
  )
)

# Summary and pairs plot
poi = c("sigma1", "rho1", "per1", "mu", "sigma_plant", "sigma_nug")

poi_samp_df = extract(fit_times_out, poi) %>% as_tibble()

summary(fit_times_out, pars=poi)

pairs(fit_times_out, pars=poi)

# Make sure you have more than 1 chain for these
stan_diag(fit_times_out, information="divergence")

stan_diag(fit_times_out, information="stepsize")

stan_diag(fit_times_out, information="sample")

stan_mcse(fit_times_out, pars=poi)


# Extract posterior samples of interest
# time_var_samp = extract(fit_times_out, c("p1", "p2"))
time_var_samp = extract(fit_times_out, c("p2"))

time_var_q_list = lapply(
  time_var_samp,
  function(x, tg) {
    temp_df = as.data.frame(t(apply(x, 2, quantile, c(0.25, 0.5, 0.75))))
    temp_df$time_grid = tg
    temp_df
  },
  time_grid
)

time_var_df = bind_rows(time_var_q_list, .id = "var")
head(time_var_df)


ggplot(time_var_df, aes(time_grid, `50%`)) +
  geom_line(aes(color=var)) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`, color=var), alpha=0.25)


# Compute posterior predictive sample for `n_plants`
sigma_plant = extract(fit_times_out, "sigma_plant")[[1]]
sigma_nug = extract(fit_times_out, "sigma_nug")[[1]]

n_plants = 100
Z1 = extract(fit_times_out, "lods2")[[1]]
dim_Z1 = dim(Z1)

sigma_error_array = array(sigma_plant, dim=c(dim_Z1, n_plants))
# sigma_error_array = array(sqrt(sigma_plant^2 + sigma_nug^2), dim=c(dim_Z1, n_plants))

sigma_error_array_dim = dim(sigma_error_array)

Z2 = array(rnorm(prod(sigma_error_array_dim), mean=0.0, sd=sigma_error_array), dim=sigma_error_array_dim)

binom_probs = Z1

plant_binom_probs = 1 / (1 + exp(-sweep(Z2, 1:2, Z1, FUN="+")))

plant_binom_means_25 = apply(plant_binom_probs, 2, quantile, 0.25) * max_trials
plant_binom_means_75 = apply(plant_binom_probs, 2, quantile, 0.75) * max_trials


# Put all of this together in a summary data frame
posterior_df = data.frame(
  time_grid = time_grid,
  log_odds_med = apply(Z1, 2, median),
  # p1_25 = time_var_q_list[["p1"]][["25%"]],
  # p1_50 = time_var_q_list[["p1"]][["50%"]],
  # p1_75 = time_var_q_list[["p1"]][["75%"]],
  p2_25 = time_var_q_list[["p2"]][["25%"]],
  p2_50 = time_var_q_list[["p2"]][["50%"]],
  p2_75 = time_var_q_list[["p2"]][["75%"]],
  plant_m_25 = plant_binom_means_25,
  plant_m_75 = plant_binom_means_75
) %>% mutate(
  # m1_25 = max_trials * p1_25,
  # m1_50 = max_trials * p1_50,
  # m1_75 = max_trials * p1_75,
  m2_25 = max_trials * p2_25,
  m2_50 = max_trials * p2_50,
  m2_75 = max_trials * p2_75
)
  

head(posterior_df)


# Plot empirical means versus log odds median
p_log_odds <-
  ggplot(posterior_df) +
  geom_point(aes(time_grid, log_odds_med))

grid.arrange(time_plot, p_log_odds, ncol=1)


## # Plot binomial mean as well as plant-based uncertainty
## p_binom_mean_1 = ggplot(posterior_df, aes(x=time_grid, y=m1_50)) +
##   geom_line() +
##   geom_ribbon(aes(ymin=m1_25, ymax=m1_75), col="gray", alpha=0.5) +
##   geom_ribbon(aes(ymin=plant_m_25, ymax=plant_m_75), col="gray", alpha=0.25) +
##   ggtitle(paste(
##     "Posterior median (black) and interquartile range (dark gray) of ",
##     "time-varying\nmean counts per day, interquartile range (light gray) ",
##     "of individual plant\nmean counts per day, and daily empirical mean root count (red)",
##     sep=""
##   )) +
##   ylab("count") +
##   xlab("days since start") +
##   geom_point(data=df_time, aes(time_grid, rate), color="red")

## p_binom_mean_1


p_binom_mean_2 = ggplot(posterior_df, aes(x=time_grid, y=m2_50)) +
  geom_line() +
  geom_ribbon(aes(ymin=m2_25, ymax=m2_75), col="gray", alpha=0.5) +
  geom_ribbon(aes(ymin=plant_m_25, ymax=plant_m_75), col="gray", alpha=0.25) +
  ggtitle(paste(
    "Posterior median (black) and interquartile range (dark gray) of ",
    "time-varying\nmean counts per day, interquartile range (light gray) ",
    "of individual plant\nmean counts per day, and daily empirical mean root count (red)",
    sep=""
  )) +
  ylab("count") +
  xlab("days since start") +
  geom_point(data=df_time, aes(time_grid, rate), color="red")

p_binom_mean_2

p_binom_mean_file = sprintf("fit_times_real-p_binom_mean-%s.png", CROP)
ggsave(p_binom_mean_2, file=file.path("images", "xm_model", p_binom_mean_file))


p_period_hist = ggplot(poi_samp_df, aes(x=per1)) +
  geom_histogram(aes(y=after_stat(density)), bins=20, alpha=0.5) +
  geom_density() +
  xlab("period") +
  ggtitle("Posterior distribution of period")


p_period_hist


# grid.arrange(p_binom_mean_2, p_period_hist)

p_binom_mean_and_per_hist = grid.arrange(p_binom_mean_2, p_period_hist, widths=c(2,1))

p_both_file = sprintf("p_binom_mean_and_per_hist-%s.png", CROP)
ggsave(p_binom_mean_and_per_hist, file=file.path("images", "xm_model", p_both_file), width=10, height=4, units="in")


# In case you want to examine residuals...
ep_plant = extract(fit_times_out, "ep_plant")[[1]]

dim(ep_plant)

ep_plant_med = apply(ep_plant, 2, median)

df_dev_time$ep_plant_med = ep_plant_med


ggplot(df_dev_time) +
  geom_line(aes(time_group, ep_plant_med, color=device_id))


ggplot(df_dev_time) +
  geom_histogram(aes(x=ep_plant_med), bins=10) +
  facet_wrap(~ device_id)


# Save posterior predictive for later use

data_file = sprintf("fit_times-%s.RData", CROP)
save(time_grid, fit_times_out, binom_probs, plant_binom_probs, max_trials, file=file.path("cache", data_file))
