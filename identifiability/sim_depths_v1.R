# -*- ess-style: RStudio; -*-

library("rstan")
library("ggplot2")
library("latex2exp")
library("dplyr")
library("moments")

source(file.path("shared", "functions.R"))



####################
# Simulate
####################


# Model - make sure it compiles
sim_depths_one = stan_model(file.path("identifiability", "sim_depths_v1.stan"), model_name="sim_depths_v1", verbose=TRUE)

# Data
dat = list(
  N = 100,
  K=3,
  r=8,
  mu_dx=2.0,
  sig_dx=1.5,
  mu_dm=-0.25,
  sig_dm=0.1
)

# Make sure model runs
out = sampling(sim_depths_one, data=dat, algorithm="Fixed_param", iter=1, chains=1)

# Plot
y = extract(out, "y")[[1]]

hist(y)


# Parameter exploration

# mu_dx and mu_dm are more important probably, since there is a clear
# tradeoff between these two for getting roots in a reasonable range.
grid_info = list()
grid_info[["mu_dx"]] = list(from=1, to=5, length.out=5)
grid_info[["sig_dx"]] = list(from=0.1, to=0.3, length.out=3)
grid_info[["mu_dm"]] = list(from=-1.0, to=-0.1, length.out=5)
grid_info[["sig_dm"]] = list(from=0.1, to=0.3, length.out=3)

grid = list()
for (nm in names(grid_info)) {
  grid[[nm]] = with(grid_info[[nm]], seq(from, to, length.out=length.out))
}

grid_df = expand.grid(grid)
grid_df


# Sample over parameter space
N = 1000
y_list = list()
sim_depths_one_samp_list = list()


for (i in 1:nrow(grid_df)) {
  if (TRUE) {
    # if (summary_df$sd[i] < 1e-2) {
    dat = with(
      grid_df,
      list(N=N, K=3, r=8, mu_dx=mu_dx[i], sig_dx=sig_dx[i], mu_dm=mu_dm[i], sig_dm=sig_dm[i])
    )
    sim_depths_one_samp_list[[i]] = sampling(
      sim_depths_one,
      data=dat,
      algorithm="Fixed_param",
      iter=1,
      chains=1,
      refresh=FALSE
    )
  }
}


# Construct sample paths
path_list = list()
path_grid = seq(0, 8, length.out = 10)

for (i in 1:nrow(grid_df)) {
  temp = lapply(extract(sim_depths_one_samp_list[[i]], c("y", "x", "dm_0", "dm")), drop)
  y_list[[i]] = temp[["y"]]
  d = dim(temp[["x"]])
  X = array(temp[["x"]], dim=c(1, d[1], d[2]))
  DM_0 = array(temp[["dm_0"]], dim=c(1, d[1]))
  DM = array(temp[["dm"]], dim=c(1, d[1], d[2]))
  paths = drop(reconstruct_path(X, DM_0, DM, path_grid))
  path_list[[i]] = paths
}


# Summary statistics
summary_df = cbind(grid_df, as.data.frame(t(sapply(y_list, summary))))

summary_df$sd = sapply(y_list, sd)
sum(summary_df$sd < 1e-2)

head(summary_df)

# Come up with a preliminary def of what is "ok"
summary_df$ok = with(summary_df, (`1st Qu.` > -16) & (Min. > -25) & (`3rd Qu.` < 0))

dim(summary_df)

# (Re-)define what constitutes reasonable data, set indicator
summary_df$ok = with(summary_df, (`1st Qu.` > -16) & (Min. > -25) & (`3rd Qu.` < 0))

summary_df$sig_dx_dm = with(summary_df, paste(sig_dx, sig_dm, sep=","))

head(summary_df, n=5)


# Save the data for fitting
save(y_list, param, file=file.path("cache", "ident-sim.RData"))


####################
# Plot
####################


# Plot for mu_dx and mu_dm by min or what we have defined as "ok"
# ggplot(summary_df, aes(mu_dx, mu_dm, color=Min.)) +
#   geom_jitter()

p_ok = ggplot(summary_df, aes(mu_dx, mu_dm, color=ok, shape=sig_dx_dm)) +
  geom_jitter(width=0.1, height=0.01) +
  ggtitle(TeX("Acceptable $\\mu_{dx}$ and $\\mu_{dm}$ parameter combinations")) +
  xlab(TeX("$\\mu_{dx}$")) +
  ylab(TeX("$\\mu_{dm}$")) +
  guides(shape=guide_legend(TeX("$\\sigma_{dx},\\sigma_{dm}$")))

p_ok

ggsave(p_ok, file=file.path("images", "xm_model", "sim_depths_one_plots-ok.png"))


# Analysis --- what would constitute reasonable prior values?

# From `p_ok` above we see there is a pretty clear boundary between
# what is acceptable or not.  As mu_dx becomes smaller, mu_dm must
# becomes smaller in magnitude.

# The model below quantifies that boundary.
param_mod_logit = glm(
  ok ~ mu_dx + mu_dm,
  family=binomial(link="logit"),
  data=subset(summary_df, mu_dm < -.2)
)

summary(param_mod_logit)

# Extract the hyperplane boundary we are interested in.  In z =
# logit(p) space, we have z = m3c[1] + m3c[2] * mu_dx + m3c[3] *
# mu_dm.  So for the level set z=0, we have mu_dm = -m3c[1]/m3c[3] -
# m3c[2]/m3c[3] * mu_dx.
m3c = coef(param_mod_logit)

p_ok_slope = -m3c[2]/m3c[3]
p_ok_slope

p_ok_int = -m3c[1]/m3c[3]
p_ok_int

p_ok_with_line = p_ok + 
  geom_abline(slope=p_ok_slope, intercept=p_ok_int)

p_ok_with_line

ggsave(p_ok, file=file.path("images", "xm_model", "sim_depths_one_plots-ok_with_line.png"))


# We could also do a regressionon some quantile.
param_mod_min = lm(Min. ~ mu_dx + mu_dm, data=summary_df)

summary(param_mod_min)


# Distributions by depth... first set up a dataframe for plotting.

ylim = range(sapply(y_list, range))
breaks = seq(ylim[1], ylim[2], length.out=100)

y_plot_list = list()
for (i in 1:length(y_list)) {
  y_plot_list[[i]] = data.frame(
    y = y_list[[i]],
    mu_dx = grid_df$mu_dx[i],
    mu_dm = grid_df$mu_dm[i],
    sig_dx = grid_df$sig_dx[i],
    sig_dm = grid_df$sig_dm[i],
    ok = summary_df$ok[i]
  )
}

y_plot_df = do.call(rbind, y_plot_list)

y_plot_df$sig_dx_dm = with(y_plot_df, paste(sig_dx, sig_dm, sep=","))

head(y_plot_df)


# Plot all together --- at the end of the file are examples of how to subset this.
p_y_all = y_plot_df %>%
  filter(ok) %>%
  ggplot() +
  geom_histogram(aes(x=y, y=..density.., fill=sig_dx_dm)) +
  facet_grid(col=vars(mu_dm), row=vars(mu_dx)) +
  xlab("depth (cm)") +
  guides(fill=guide_legend(TeX("$\\sigma_{dx},\\sigma_{dm}$"))) +
  ggtitle(TeX("Distribution of depths by $(\\mu_{dm},\\mu_{dx})$")) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_y_all

ggsave(p_y_all, file=file.path("images", "xm_model", "sim_depths_one_plots-y_all.png"))


# When do distributions look the same?  Let's look at the mean and
# standard deviation.  You can see a tradeoff in mu_dx for mu_dm can
# result in similar values.
y_plot_df_ss = y_plot_df %>%
  group_by(mu_dx, mu_dm, sig_dx, sig_dm) %>%
  summarize(
    mean = mean(y),
    std = sd(y),
    m3 = skewness(y),
    m4 = kurtosis(y),
    n = n(),
  ) %>%
  ungroup() %>%
  ## mutate(
  ##   mean_cut = cut(mean, breaks=5, ordered_result=TRUE),
  ##   std_cut = cut(std, breaks=5, ordered_result=TRUE)
  ## ) %>%
  arrange(mean, std)

print(y_plot_df_ss, n=45)


# Plot mean by mu_dx and mu_dm
p_mean_depth = y_plot_df_ss %>%
  filter(mean < 0.5) %>%
  ggplot(aes(mu_dx, mu_dm)) +
  geom_tile(aes(fill=mean)) +
  scale_fill_gradient2(midpoint=-10) +
  facet_grid(col=vars(sig_dx), row=vars(sig_dm)) +
  ggtitle(TeX("Mean depth for $\\mu_{dx}, \\mu_{dm}$ stratified by $\\sigma_{dx},\\sigma_{dm}$"))

p_mean_depth

ggsave(p_mean_depth, file=file.path("images", "xm_model", "sim_depths_one_plots-mean_depth.png"))


# Plot std by sig_dx and sig_dm
p_std_depth = y_plot_df_ss %>%
  filter(mean < 0) %>%
  # ggplot(aes(mu_dx, mu_dm)) +
  ggplot(aes(sig_dx, sig_dm)) +
  geom_tile(aes(fill=std)) +
  scale_fill_gradient2(midpoint=median(y_plot_df_ss$std)) +
  # facet_grid(col=vars(sig_dx), row=vars(sig_dm)) +
  facet_grid(col=vars(mu_dx), row=vars(mu_dm)) +
  # ggtitle(TeX("Std. dev. depth for $\\mu_{dx}, \\mu_{dm}$ stratified by $\\sigma_{dx},\\sigma_{dm}$"))
  ggtitle(TeX("Std. dev. depth for $\\sigma_{dx}, \\sigma_{dm}$ stratified by $\\mu_{dx},\\mu_{dm}$"))

p_std_depth

ggsave(p_std_depth, file=file.path("images", "xm_model", "sim_depths_one_plots-std_depth.png"))


# Plt std by mu_dx and sig_dm alternative
y_plot_df_ss %>%
  filter(mean < 0) %>%
  # ggplot(aes(mu_dx, mu_dm)) +
  ggplot(aes(mu_dx, sig_dm)) +
  geom_tile(aes(fill=std)) +
  scale_fill_gradient2(midpoint=0.53) +
  # facet_grid(col=vars(sig_dx), row=vars(sig_dm)) +
  facet_grid(col=vars(sig_dx), row=vars(mu_dm)) +
  # ggtitle(TeX("Std. dev. depth for $\\mu_{dx}, \\mu_{dm}$ stratified by $\\sigma_{dx},\\sigma_{dm}$"))
  ggtitle(TeX("Std. dev. depth for $\\mu_{dx}, \\sigma_{dm}$ stratified by $\\sigma_{dx},\\mu_{dm}$"))



# Can we quantify those level sets via lines? As we saw with our "ok"
# parameters, there are level sets for the mean that trade mu_dx for
# mu_dm.  For the sd, mu_dx, sig_dx, and sig_dm play a large role.
# While this is only based on two moments, this seems to be telling us
# that there will be flat portions in our posterior.  Again, we see a
# situation where a change in 12 units of mu_dx is equivalent to a
# change in 1 unit of mu_dm.

mean_mod_1 = lm(mean ~ mu_dx + mu_dm + sig_dx + sig_dm, data=y_plot_df_ss)

summary(mean_mod_1)

mean_mod_2 = lm(mean ~ mu_dx + mu_dm, data=y_plot_df_ss)

summary(mean_mod_2)

std_mod_1 = lm(std ~ mu_dx + mu_dm + sig_dx + sig_dm, data=y_plot_df_ss)

summary(std_mod_1)

# Add back the fitted values --- this doesn't completely clear things up because the 
y_plot_df_ss$fitted_mean = fitted(mean_mod_2)
y_plot_df_ss$fitted_std = fitted(std_mod_1)

print(y_plot_df_ss %>% arrange(fitted_mean, fitted_std), n=45)



# These plots DEPEND ON PARAMETERS SET IN grid_df!!!

# For varying mu_dx only
p_y_mu_dx = y_plot_df %>%
  filter(ok, mu_dm == -0.55, sig_dx == 0.1, sig_dm == 0.2) %>%
  ggplot() +
  geom_histogram(aes(y)) +
  coord_flip() +
  facet_grid(col=vars(mu_dx)) +
  xlab("depth (cm)")

p_y_mu_dx


# For varying mu_dm only
p_y_mu_dm = y_plot_df %>%
  filter(ok, mu_dx == 3, sig_dx == 0.3, sig_dm == 0.2) %>%
  ggplot() +
  geom_histogram(aes(y)) +
  coord_flip() +
  facet_grid(col=vars(mu_dm)) +
  xlab("depth (cm)")

p_y_mu_dm


# For varying sig_dx, sig_dm only
p_y_sig_dx_dm = y_plot_df %>%
  filter(ok, mu_dx == 3, mu_dm == -0.6) %>%
  ggplot() +
  geom_histogram(aes(y)) +
  coord_flip() +
  facet_grid(row=vars(sig_dm), col=vars(sig_dx)) +
  xlab("depth (cm)")

p_y_sig_dx_dm



