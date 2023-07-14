# -*- ess-style: RStudio -*- #
# Jesse Windle

library("ggplot2")
library("gridExtra")
library("yaml")
library("rstan")

source("shared/functions.R")


# CONFIG
config = read_yaml("xm_model/depth-m23-fit-corn.yaml")
CROP = config$crop
r = config$r
K = config$K
mu_dx = r / (K+1)


# Load
df_all = load_rt_data(file.path("data", "acc-1-v3.csv")) %>%
    filter(uptime > 0, electrode_pair > 1, days_since_start > 4, crop == CROP)

df_all$old_count = df_all$count
df_all$count = 1 * (df_all$old_count > 0)

time_grid = unique(df_all$days_since_start)
df_all$time_group =  match(df_all$days_since_start, time_grid)


# Data for fitting
depths_df = df_all %>%
    filter(
        count > 0,
        crop == CROP
    )

y_data = depths_df$depth_cm
n_obs = length(y_data)
time_group = depths_df$time_group
n_times = max(time_group)


# Model
mdl_23 = stan_model("xm_model/depth-m23-fit.stan", model_name="depth-m23-fit", verbose=TRUE)


# Data
base_dat = list(
  n_obs = n_obs,
  n_times = n_times,
  time_grid = 1:n_times,
  time_group = time_group,
  r = r,
  mu_dx = mu_dx,
  K = K,
  y = y_data
)

dat = c(base_dat, config$prior)
dat


# Fit - had to up adapt_delta to 0.95 to git rid of divergent
# transitions.
out_23 = sampling(
  mdl_23,
  data=dat,
  iter=config$mcmc$iter,
  warmup=config$mcmc$warmup,
  chains=config$mcmc$chains,
  cores=config$mcmc$cores,
  control = list(
    adapt_delta = 0.95,
    adapt_t0 = 10
  )
)

# Summary and pairs plot
static_pars = c("mu_y_gbl_mean", "sig_y_gbl_mean", "shape_dm_gbl_mean")

proc_pars = c("mu_y", "sig_y", "shape_dm")

# samp_df = extract(out_23, static_pars) %>% as_tibble()

summary(out_23, pars=static_pars)

pairs(out_23, pars=static_pars)


# Diagnositcs
stan_diag(out_23, information="divergence")

stan_diag(out_23, information="stepsize")

stan_diag(out_23, information="sample")

stan_mcse(out_23, pars=poi)


# Proccesses
f_gp = extract(out_23, "f_gp")[[1]]

dim(f_gp)
f_gp_mean = apply(f_gp, c(2, 3), mean)

par(mfrow=c(1,3))
plot(f_gp_mean[1,], main="mu_y deviations", type="l")
plot(f_gp_mean[2,], main="sig_y deviations", type="l")
plot(f_gp_mean[3,], main="shape_dm deviations", type="l")


# Extract posterior samples of interest
param_procs = extract(out_23, c("mu_y", "sig_y", "shape_dm"))

mu_dm = extract(out_23, "mu_dm")[[1]]
sig_dm = extract(out_23, "sig_dm")[[1]]
shape_dm = param_procs[["shape_dm"]]

time_var_q_list = lapply(
  param_procs,
  function(x, tg) {
    temp_df = as.data.frame(t(apply(x, 2, quantile, c(0.25, 0.5, 0.75))))
    temp_df$time_grid = tg
    temp_df
  },
  time_grid
)

time_var_df = bind_rows(time_var_q_list, .id = "var")
head(time_var_df)


# Plot mean process
time_var_df %>% filter(var == "mu_y") %>%
  ggplot(aes(time_grid, `50%`)) +
  geom_line() +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), alpha=0.25)


# Plot sd process
time_var_df %>% filter(var == "sig_y") %>%
  ggplot(aes(time_grid, `50%`)) +
  geom_line() +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), alpha=0.25)


# Plot shape process
time_var_df %>% filter(var == "shape_dm") %>%
  ggplot(aes(time_grid, `50%`)) +
  geom_line() +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), alpha=0.25)



# Save for generating plant roots

file_name = sprintf("m23-fit-%s.RData", CROP)
save(mu_dm, sig_dm, shape_dm, param_procs, out_23, file=file.path("cache", file_name))
