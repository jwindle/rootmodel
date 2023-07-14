# -* ess-style: RStudio -*-

library("plot3D")
library("yaml")
# You need tidyr, but do not load it

source("shared/functions.R")


# Config
EPOCHS = c(0, 10, 18, 22, 28)
n_exampls = 10

config = read_yaml("xm_model/depth-m23-fit-corn.yaml")
CROP = config$crop
r = config$r
K = config$K
refine_factor = config$refine_factor
mu_dx = r / (K+1)


# Load
fit_times_name = sprintf("fit_times-%s.RData", CROP)
load(file.path("cache", fit_times_name))

depth_m23_name = sprintf("m23-fit-%s.RData", CROP)
load(file.path("cache", depth_m23_name))


# Transform
plant_binom_probs_array = aperm(plant_binom_probs, c(1, 3, 2))

plant_binom_probs_mat = matrix(
  as.numeric(plant_binom_probs_array),
  ncol=dim(plant_binom_probs_array)[3],
  byrow=FALSE
)

dim(plant_binom_probs_mat)


# Build model --- we will pass this to a function.
sim_plant_mdl = stan_model(file.path("xm_model", "sim_plant.stan"), model_name="sim_plant", verbose=TRUE)


# Let's try simulating
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


# Simulate many plants

plants_df = sim_plants(
  sim_plant_mdl,
  N = n_examples,
  p = plant_binom_probs_mat,
  mu_dm = mu_dm,
  sig_dm = sig_dm,
  shape_dm = shape_dm,
  time_grid = time_grid,
  max_trials = max_trials * 0.5,
  r = r,
  mu_dx = mu_dx
)

plants_df$depth_0 = 0.0
plants_df$time_bin = cut(plants_df$time, breaks=EPOCHS, include.lowest=TRUE, ordered_result=TRUE)

head(plants_df)

plants_long_df = plants_df %>%
  filter(depth_1 < 0, depth_2 < 0) %>%
  tidyr::pivot_longer(
    cols = starts_with("depth"),
    names_prefix = "depth_",
    names_to="kink",
    values_to = "depth"
  ) %>%
  mutate(
    kink = as.integer(kink),
    rad = mu_dx * kink
  )

head(plants_long_df)


# 2-D view of roots over time
p_plants_over_time =
  ggplot(plants_long_df, aes(rad, depth)) +
  geom_line(aes(group=root_id)) +
  facet_grid(cols=vars(plant), rows=vars(time_bin)) +
  ggtitle("Example root trajectory by epoch of emergence") +
  xlab("radius (cm)") +
  ylab("depth (cm)")

p_plants_over_time

ggsave(
  p_plants_over_time,
  file = file.path("images", "xm_model", sprintf("p_plants_over_time-%s.png", CROP)),
  width = 6,
  height = 6,
  units = "in"
)


# 2-D view of roots over time
p_plants_final =
  ggplot(plants_long_df, aes(rad, depth)) +
  geom_line(aes(group=paste(plant, root_id))) +
  facet_grid(cols=vars(plant)) +
  ggtitle("Example root growth") +
  xlab("radius (cm)") +
  ylab("depth (cm)")

p_plants_final

ggsave(
  p_plants_final,
  file = file.path("images", "xm_model", sprintf("p_plants_final-%s.png", CROP)),
  width = 6,
  height = 3,
  units = "in"
)



# Simulate canonical plants

p_mean = apply(plant_binom_probs_mat, 2, mean)
mu_dm_mean = apply(mu_dm, 2, mean)
sig_dm_mean = apply(sig_dm, 2, mean)
shape_dm_mean = apply(shape_dm, 2, mean)

canonical_df = sim_plants(
  sim_plant_mdl,
  N = n_examples,
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
canonical_df$time_bin = cut(canonical_df$time, breaks=EPOCHS, include.lowest=TRUE, ordered_result=TRUE)

head(canonical_df)

canonical_long_df = cdf %>%
  filter(depth_1 < 0, depth_2 < 0) %>%
  tidyr::pivot_longer(
    cols = starts_with("depth"),
    names_prefix = "depth_",
    names_to="kink",
    values_to = "depth"
  ) %>%
  mutate(
    kink = as.numeric(kink),
    rad = mu_dx * kink / refine_factor,
    x = rad * cos(angle),
    y = rad * sin(angle),
    kink = as.integer(kink),
    time_hours = time_index * 24 + starting_hour + kink
  ) %>%
  arrange(plant, root_id, time_index, kink)

head(canonical_long_df, n=10)


# 2-D view of roots over time
p_canon_plants_over_time =
  ggplot(canonical_long_df, aes(rad, depth)) +
  geom_line(aes(group=root_id)) +
  facet_grid(cols=vars(plant), rows=vars(time_bin)) +
  ggtitle("Example canonical root trajectory by epoch of emergence") +
  xlab("radius (cm)") +
  ylab("depth (cm)")

p_canon_plants_over_time

ggsave(
  p_canon_plants_over_time,
  file = file.path("images", "xm_model", sprintf("p_canon_plants_over_time-%s.png", CROP)),
  width = 6,
  height = 6,
  units = "in"
)


# 2-D view of roots over time
p_canon_plants_final = 
  ggplot(canonical_long_df, aes(rad, depth)) +
  geom_line(aes(group=paste(plant, root_id))) +
  facet_grid(cols=vars(plant)) +
  ggtitle("Example canonical root growth") +
  xlab("radius (cm)") +
  ylab("depth (cm)")
  
p_canon_plants_final

ggsave(
  p_canon_plants_final,
  file = file.path("images", "xm_model", sprintf("p_canon_plants_final-%s.png", CROP)),
  width = 6,
  height = 3,
  units = "in"
)



# Refine paths and make data frame for constructing animations

temp_copy = canonical_df
temp_copy$depth_0 = NULL

cdf_1 = temp_copy %>% select(-starts_with("depth_"))
cdf_2 = temp_copy %>% select(starts_with("depth_"))

paths_rough = cbind(0.0, as.matrix(cdf_2))
head(paths_rough)

paths_fine = refine_paths_matrix(paths_rough, mu_dx, refine_factor)
colnames(paths_fine) = paste("depth", 0:(refine_factor*(K+1)), sep="_")
head(paths_fine)

cdf = cbind(cdf_1, paths_fine)
head(cdf)


movie_df = cdf %>%
  filter(depth_1 < 0, depth_2 < 0) %>%
  tidyr::pivot_longer(
    cols = starts_with("depth"),
    names_prefix = "depth_",
    names_to="rad_index",
    values_to = "depth"
  ) %>%
  mutate(
    rad_index = as.integer(rad_index),
    kink = rad_index / refine_factor,
    rad = mu_dx * rad_index / refine_factor,
    x = rad * cos(angle),
    y = rad * sin(angle),
    time_hours = time_index * 24 + starting_hour + rad_index
  ) %>%
  arrange(plant, root_id, time_index, kink)

head(movie_df)


# Plotting the first kink, as an example.
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


# Now we plot the root trajectories.  Here, we just view things.
# Below, we actually write files to create an animated GIF.
total_hours = max(movie_df$time_hours)
total_hours

# par(mfrow = c(1,1))
par(mar = c(0, 0, 3, 0))
segments3D(
  x0 = 0,
  y0 = 0,
  z0 = 0,
  x1 = 0,
  y1 = 0,
  z1 = 0,
  xlim=c(-r,r),
  ylim=c(-r,r),
  zlim=c(-18, 0),
  main = sprintf("Simulations of the canonical %s root system from experiment", tolower(CROP)),
  ticktype="detailed"
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
}


# We have to do a double loop because I can't seem to get rgl.snapshot
# or rgl.postscript to work.
for (ell in 1:n_examples) {
  for (j in 1:total_hours) {
    file_name = sprintf("%s-%02d-%04d.png", CROP, ell, i)
    png(file.path("images", "frames", file_name))
    # par(mfrow = c(1,1))
    par(mar = c(0, 0, 3, 0))
    segments3D(
      x0 = 0,
      y0 = 0,
      z0 = 0,
      x1 = 0,
      y1 = 0,
      z1 = 0,
      xlim=c(-r,r),
      ylim=c(-r,r),
      zlim=c(-18, 0),
      main = sprintf("Simulations of the canonical %s root system from experiment", tolower(CROP)),
      ticktype="detailed"
    )
    one_plant_df = movie_df %>% filter(plant == ell)
    for (i in 1:j) {
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
    }
    dev.off()
  }
}



# write.csv(movie_df, file=file.path("cache", "sim_plants_can.csv"), row.names=FALSE)


