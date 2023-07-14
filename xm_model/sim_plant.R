# -* ess-style: RStudio -*-

library("plot3D")
library("rgl")
# You need tidyr, but do not load it

source("shared/functions.R")

# Load
load(file.path("cache", "fit_times.RData"))

load(file.path("cache", "depth-m23-fit.RData"))


# Transform
plant_binom_probs_array = aperm(plant_binom_probs, c(1, 3, 2))

plant_binom_probs_mat = matrix(
  as.numeric(plant_binom_probs_array),
  ncol=dim(plant_binom_probs_array)[3],
  byrow=FALSE
)

dim(plant_binom_probs_mat)


# Need to define some stuff or have it defined
mu_dx = 4
K = 1
r = 8


# Make sure model compiles
sim_plant_mdl = stan_model(file.path("xm_model", "sim_plant.stan"), model_name="sim_plant", verbose=TRUE)


# Let's test this
p_ave = apply(plant_binom_probs, 2, mean)

num_roots = rbinom(23, p=p_ave, size=max_trials)
root_times = rep(1:23, times=num_roots)

# Check
all(table(root_times) == num_roots[num_roots > 0])

# Data
max_trials = 48
dat = list(
  mu_dx = mu_dx,
  K = K,
  r = r,
  n_times = length(time_grid),
  n_observed_times = length(root_times),
  time_index = root_times,
  mu_dm = mu_dm[1,],
  sig_dm = sig_dm[1,],
  shape_dm = shape_dm[1,]
)

sim_plant_out = sampling(
  sim_plant_mdl,
  data=dat,
  iter=1,
  warmup=0,
  chains=1,
  cores=1,
  algorithm="Fixed_param"
)



# Simulate many different plants

plants_df = sim_plants(
  sim_plant_mdl,
  N = 10,
  p = plant_binom_probs_mat,
  mu_dm = mu_dm,
  sig_dm = sig_dm,
  shape_dm = shape_dm,
  time_grid = time_grid,
  max_trials = max_trials,
  r = r,
  mu_dx = mu_dx
)

plants_df$depth_0 = 0.0
plants_df$time_bin = cut(plants_df$time_index, breaks=4)

head(plants_df)


plants_long_df = plants_df %>%
  filter(depth_1 < 0, depth_2 < 0) %>%
  tidyr::pivot_longer(cols = starts_with("depth"), names_prefix = "depth_", names_to="kink", values_to = "depth")

plants_long_df$kink = as.integer(plants_long_df$kink)

head(plants_long_df)

# 2D view of roots over time
ggplot(plants_long_df, aes(kink, depth)) +
  geom_line(aes(group=root_id)) +
  facet_grid(rows=vars(plant), cols=vars(time_bin))



# Simulate canonical plant

p_mean = apply(plant_binom_probs_mat, 2, mean)
mu_dm_mean = apply(mu_dm, 2, mean)
sig_dm_mean = apply(sig_dm, 2, mean)
shape_dm_mean = apply(shape_dm, 2, mean)

canonical_df = sim_plants(
  sim_plant_mdl,
  N = 10,
  p = p_mean,
  mu_dm = mu_dm_mean,
  sig_dm = sig_dm_mean,
  shape_dm = shape_dm_mean,
  time_grid = time_grid,
  max_trials = max_trials * 0.5,
  r = r,
  mu_dx = mu_dx
)

canonical_df$depth_0 = 0.0
canonical_df$time_bin = cut(canonical_df$time_index, breaks=4)

head(canonical_df)


canonical_long_df = canonical_df %>%
  filter(depth_1 < 0, depth_2 < 0) %>%
  tidyr::pivot_longer(cols = starts_with("depth"), names_prefix = "depth_", names_to="kink", values_to = "depth")

canonical_long_df$kink = as.integer(canonical_long_df$kink)


# 2-D view of roots over time
ggplot(canonical_long_df, aes(kink, depth)) +
  geom_line(aes(group=root_id)) +
  facet_grid(cols=vars(plant), rows=vars(time_bin))


# 2-D view of roots over time
ggplot(canonical_long_df, aes(kink, depth)) +
  geom_line(aes(group=paste(plant, root_id))) +
  facet_grid(cols=vars(plant))



# Series of 3D images to write out


rad = c(rad_0=0, rad_1=4, rad_2=8)

temp_df = cbind(canonical_df, outer(rep(1, nrow(canonical_df)), rad))

movie_df = temp_df %>%
  filter(depth_1 < 0, depth_2 < 0) %>%
  tidyr::pivot_longer(
    cols = starts_with("depth"),
    names_prefix = "depth_",
    names_to="kink",
    values_to = "depth"
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("rad"),
    names_prefix = "rad_",
    names_to="kink_rad",
    values_to = "rad"
  ) %>%
  filter(kink == kink_rad) %>%
  mutate(
    x = rad * cos(angle),
    y = rad * sin(angle),
    kink = as.integer(kink),
    time_hours = time_index * 24 + starting_hour + kink
  ) %>%
  arrange(plant, root_id, time_index, kink)

movie_df$kink_rad = NULL

head(movie_df)


# First kink
df_cur = movie_df %>% filter(kink == 1)
df_prev = movie_df %>% filter(kink == 0)

par(mfrow = c(1,1))
segments3D(
  x0 = df_prev$x,
  y0 = df_prev$y,
  z0 = df_prev$depth,
  x1 = df_cur$x,
  y1 = df_cur$y,
  z1 = df_cur$depth
)


total_hours = max(movie_df$time_hours)
total_hours



# par(mfrow = c(1,1))
# open3d(windowRect=c(100,100,700,700))
segments3D(
  x0 = 0,
  y0 = 0,
  z0 = 0,
  x1 = 0,
  y1 = 0,
  z1 = 0,
  xlim=c(-8,8),
  ylim=c(-8,8),
  zlim=c(-16, 0)
)
one_plant_df = movie_df %>% filter(plant == 4)
for (i in 1:total_hours) {
  file_name = sprintf("plant-%04d.png", i)
  df_cur = one_plant_df %>% filter(time_hours == i, kink > 0)
  df_prev = one_plant_df %>% filter(time_hours == (i-1), kink < (K+1))
  if (nrow(df_cur) > 0) {
    segments3D(
      x0 = df_prev$x,
      y0 = df_prev$y,
      z0 = df_prev$depth,
      x1 = df_cur$x,
      y1 = df_cur$y,
      z1 = df_cur$depth,
      add = TRUE
    )
  }
  # rgl.snapshot(file.path("images", "frames", file_name), fmt = "png")
  # rgl.postscript(file.path("images", "frames", file_name), fmt = "ps")
}

