# -*- ess-style: RStudio; -*-

# MUST SET STARTING POINT FOR LOG TRANSFORM

source("shared/functions.R")

library("rstan")
library("bayesplot")
library("yaml")
library("gridExtra")


config = read_yaml("xm_model/depth-MAIN-fit_real-log-config.yaml")

SIMS_TO_RUN = config$simulations_to_run
SIMS_TO_RUN

PATHS_TO_PLOT = 50


# Load
load("cache/detections.RData")

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
    temp_list[[i]]$model = sim_name
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

  # Plot paths
  if (PATHS_TO_PLOT) {

    d = dim(sim_paths_refined)
    n_epoch = d[3]
    n_steps = d[2]
    x_grid = seq(1/n_steps, 1, length.out=n_steps) * 8.0
    col = "#80808040"
    path_idc = floor(seq(1, d[1], length.out = PATHS_TO_PLOT))

    png(file.path("images", "xm_model", sprintf("%s-paths.png", sim_name)))
    par(mfrow=c(2,2))
    for (i in 1:4) {
      jj = path_idc[1]
      plot_title = sprintf("Paths epoch %d", i)
      plot(x_grid, sim_paths_refined[jj,,i], ylim=c(-18., 1.), col=col, type="l", main=plot_title, xlab="x (cm)", ylab="y (cm)")
      for (j in 2:PATHS_TO_PLOT) {
        jj = path_idc[j]
        lines(x_grid, sim_paths_refined[jj,,i], col=col)
      }
    }
    dev.off()

  }
  
}


# Summarize differences

frac_pos_df = do.call(rbind, frac_positive_list)
frac_pos_df

df_comp_null = df_time_elec_corn[,c("epoch", "depth_bin", "prop")] %>%
  mutate(
    thresh=Inf,
    model = "empirical",
    prop.emp = prop,
    diff = 0,
    klt = 0
  ) %>% as.data.frame()


df_comp = do.call(rbind, c(df_comp_list, list(df_comp_null)))
df_comp$model = as.factor(df_comp$model)

head(df_comp)
tail(df_comp)

df_comp_summary = df_comp %>%
  group_by(model, thresh, epoch) %>%
  summarize(
    mae = sum(abs(diff)),
    kl = sum(klt),
    check = sum(prop)
  )

print(df_comp_summary, n=20)

df_comp_summary_summary = df_comp_summary %>%
  # filter(epoch %in% c("(10,18]", "(22,28]")) %>%
  group_by(model, thresh) %>%
  summarize(
    mmae = mean(mae),
    mkl = mean(kl)
  )

print(df_comp_summary_summary, n=30)


if (config$write) {
  comp_summary_file = sprintf("fit_real-comp_summary.csv")
  write.csv(df_comp_summary, file=file.path("cache", comp_summary_file))
  summary_summary_file = sprintf("fit_real-comp_summary-summary.csv")
  write.csv(df_comp_summary_summary, file=file.path("cache", summary_summary_file))
  frac_pos_file = sprintf("frac_pos-summary.csv")
  write.csv(frac_pos_df, file=file.path("cache", frac_pos_file))
}


p_comp_prop_no_thresh = df_comp %>%
  filter(thresh == Inf) %>%
  ggplot(aes(depth_bin, prop, color=model)) +
  geom_point() +
  geom_line(aes(group=model)) +
  coord_flip() +
  facet_wrap(~ epoch) +
  ggtitle("Distribution of depths by model or observed, all paths")

p_comp_prop_no_thresh


p_comp_prop_thresh = df_comp %>%
  filter((thresh < Inf) | (model == "empirical")) %>%
  ggplot(aes(depth_bin, prop, color=model)) +
  geom_point() +
  geom_line(aes(group=model)) +
  coord_flip() +
  facet_wrap(~ epoch) +
  ggtitle("Distribution of depths by model or observed, unrealistic paths removed")

p_comp_prop_thresh

grid.arrange(p_comp_prop_no_thresh, p_comp_prop_thresh)


if (config$write) {
  p_comp_no_thresh_file = sprintf("fit_real-p_comp_prop-no_thresh.png")
  ggsave(p_comp_prop_no_thresh, file=file.path("images", "xm_model", p_comp_no_thresh_file))
  p_comp_thresh_file = sprintf("fit_real-p_comp_prop-tresh.png")
  ggsave(p_comp_prop_thresh, file=file.path("images", "xm_model", p_comp_thresh_file))
}
