# -*- ess-style: RStudio; -*-

# MUST SET STARTING POINT FOR LOG TRANSFORM

source("shared/functions.R")

library("rstan")
library("bayesplot")
library("yaml")
library("gridExtra")


config = read_yaml("xm_model/depth-MAIN-fit_real-log-config.yaml")

SIMS_TO_RUN = config$simulations_to_run


# Load
# 
df_time_elec_corn = df_time_elec %>% filter(crop == "Corn")


# Make comparisons

df_comp_list = list()
frac_positive_list = list()

for (sim_name in SIMS_TO_RUN) {

  paths_file = sprintf("%s-paths.RData", sim_name)
  load(file=file.path("cache", paths_file))

  frac_positive_list[[sim_name]] = sim_depths_df %>%
    group_by(epoch) %>%
    summarize(
      frac_pos = mean(shallowest_cm > 1.0)
    )
  frac_positive_list[[sim_name]]$sim_name = sim_name

  threshes = c(1., Inf)
  temp_list = list()
  for (i in 1:2) {
    temp_list[[i]] = sim_depths_df %>%
      filter(shallowest_cm < threshes[i]) %>% # Eliminate unrealistic paths
      group_by(epoch, depth_bin) %>%
      summarize(n = n()) %>%
      ungroup() %>%
      filter(!is.na(depth_bin)) %>%
      group_by(epoch) %>%
      mutate(
        n_epoch = sum(n)
      ) %>%
      ungroup() %>%
      mutate(
        prop = n / n_epoch
      ) %>%
      select(epoch, depth_bin, prop)
    temp_list[[i]]$sim = sim_name
    temp_list[[i]]$thresh = threshes[i]
  }
  simulated_depths_summary = do.call(rbind, temp_list)
  
  df_comp_list[[sim_name]] = merge(
    simulated_depths_summary,
    df_time_elec_corn[, c("epoch", "depth_bin", "prop")],
    by=c("epoch", "depth_bin"),
    suffixes=c("", ".emp"),
    all.x=TRUE
  ) %>%
    mutate(
      diff = prop - prop.emp,
      klt = prop.emp * (log(prop.emp) - log(prop))
    )
  
}


# Summarize differences

frac_pos_df = do.call(rbind, frac_positive_list)
frac_pos_df

df_comp_null = df_time_elec_corn[,c("epoch", "depth_bin", "prop")] %>%
  mutate(
    thresh=Inf,
    sim = "empirical",
    prop.emp = prop,
    diff = 0,
    klt = 0
  ) %>% as.data.frame()


df_comp = do.call(rbind, c(df_comp_list, list(df_comp_null)))
df_comp$sim = as.factor(df_comp$sim)

head(df_comp)
tail(df_comp)

df_comp_summary = df_comp %>%
  group_by(sim, thresh, epoch) %>%
  summarize(
    mae = sum(abs(diff)),
    kl = sum(klt),
    check = sum(prop)
  )

print(df_comp_summary, n=20)

df_comp_summary %>%
  group_by(sim, thresh) %>%
  summarize(
    mmae = mean(mae),
    mkl = mean(kl)
  )


if (config$write) {
  comp_summary_file = sprintf("%s-fit_real-comp_summary.csv", SESSION_ID)
  write.csv(df_comp_summary, file=comp_summary_file)
  frac_pos_file = sprintf("%s-frac_pos-summary.csv", SESSION_ID)
  write.csv(frac_pos_df, file=frac_pos_file)
}


p_comp_prop_no_thresh = df_comp %>%
  filter(thresh == Inf) %>%
  ggplot(aes(depth_bin, prop, color=sim)) +
  geom_point() +
  geom_line(aes(group=sim)) +
  coord_flip() +
  facet_wrap(~ epoch) +
  ggtitle("Distribution of depths by model or observed")

p_comp_prop_no_thresh


p_comp_prop_thresh = df_comp %>%
  filter((thresh < Inf) | (sim == "empirical")) %>%
  ggplot(aes(depth_bin, prop, color=sim)) +
  geom_point() +
  geom_line(aes(group=sim)) +
  coord_flip() +
  facet_wrap(~ epoch) +
  ggtitle("Distribution of depths by model or observed")

p_comp_prop_thresh

grid.arrange(p_comp_prop_no_thresh, p_comp_prop_thresh)

## if (config$write) {
##   p_comp_prop_file = sprintf("%s-fit_real-p_comp_prop.png", SESSION_ID)
##   ggsave(p_comp_prop, file=file.path("..", "images", "xm_sim", p_comp_prop_file))
## }
