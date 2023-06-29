# -*- ess-style: RStudio; -*-


# Make data used for fitting models.

WRITE = TRUE
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
df_all = load_rt_data("data/acc-1-v3.csv") %>%
  filter(
    uptime > 0,
    electrode_pair > 1,
    days_since_start > 4
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


# Sanity check & df_time_elec used for comparisons below
depth_bins_no_pad = depth_ivals[seq(3, 23, by=4)]
depth_bins = depth_bins_no_pad
DEPTH_BINS = depth_bins

df_coarse = coarsen_rt_data(
  df_all,
    # days_since_start = seq(0, 28, 4),
  days_since_start = EPOCHS,
  depth_cm = depth_bins
)

df_time_elec = agg_rt_data(df_coarse, crop, days_since_start, depth_cm) %>%
  group_by(crop, days_since_start) %>%
  mutate(sum_rate = sum(rate)) %>%
  ungroup() %>%
  mutate(prop = rate / sum_rate)

df_time_elec$epoch = df_time_elec$days_since_start
df_time_elec$depth_bin = df_time_elec$depth_cm
df_time_elec$days_since_start = NULL
df_time_elec$depth_cm = NULL
df_time_elec$sum_rate = NULL

head(df_time_elec)

df_time_elec %>%
  ggplot() +
  geom_tile(aes(epoch, depth_bin, fill=rate)) + 
  scale_fill_gradient2() +
  facet_wrap(~ crop)

df_time_elec %>%
  filter(crop == "Corn") %>%
  ggplot() +
  geom_line(aes(-as.integer(depth_bin), prop)) +
  facet_wrap(~ epoch)

df_time_elec %>%
  ggplot(aes(depth_bin, prop, color=crop)) +
  geom_point() +
  geom_line(aes(group=crop)) +
  coord_flip() +
  facet_wrap(~ epoch)


# Moments
df_moments = df_all %>%
  filter(count > 0) %>%
  group_by(epoch, crop) %>%
  summarize(
    m1 = mean(depth_cm),
    m2 = sd(depth_cm),
    m3 = skewness(depth_cm),
    m4 = kurtosis(depth_cm),
    m4_cns = moment(depth_cm, order=4, central=TRUE)
  )
  
df_moments


if (WRITE) {
  write.csv(df_all, "data/acc-1-v3-cleaned.csv")
  write.csv(df_time_elec, "data/acc-1-v3-summary.csv")
  save(df_all, df_time_elec, DEPTH_BINS, file=file.path("cache", "detections.RData"))
}
