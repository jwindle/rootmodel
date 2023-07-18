# -*- ess-style: RStudio; -*-

# MUST SET STARTING POINT FOR LOG TRANSFORM

source("shared/functions.R")

library("rstan")
library("bayesplot")
library("yaml")
library("gridExtra")


config = read_yaml("xm_model/depth-MAIN-fit_real-config.yaml")

SIMS_TO_RUN = config$simulations_to_run
SIMS_TO_RUN

CROP = config$crop
SESSION_ID = tolower(CROP)

PATHS_TO_PLOT = 50


# Load
load("cache/detections.RData")

df_time_elec_corn = df_time_elec %>% filter(crop == CROP)


# Make comparisons

df_comp_list = list()
frac_positive_list = list()

for (sim_name in SIMS_TO_RUN) {

  paths_file = sprintf("%s-%s-paths.RData", SESSION_ID, sim_name)
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

    png(
      file.path("images", "xm_model", sprintf("%s-%s-paths.png", SESSION_ID, sim_name)),
      width=700,
      height=400,
      units="px"
    )
    
    par(mfrow=c(1,4))
    for (i in 1:4) {
      jj = path_idc[1]
      plot_title = sprintf("Paths epoch %d, %s", i, SESSION_ID)
      plot(x_grid, sim_paths_refined[jj,,i], ylim=c(-18., 1.), col=col, type="l", main=plot_title, xlab="x (cm)", ylab="y (cm)")
      for (j in 2:PATHS_TO_PLOT) {
        jj = path_idc[j]
        lines(x_grid, sim_paths_refined[jj,,i], col=col)
      }
      mtext(
        sim_name,
        side = 3,
        line = - 2,
        outer = TRUE
      )
    }
    
    dev.off()

  }
  
}



# COMPARE BY TABLES

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

print(df_comp_summary, n=44)

df_comp_summary_summary = df_comp_summary %>%
  # filter(epoch %in% c("(10,18]", "(22,28]")) %>%
  group_by(model, thresh) %>%
  summarize(
    mmae = mean(mae),
    mkl = mean(kl)
  )

print(df_comp_summary_summary %>% filter(thresh == Inf), n=40)


if (config$write) {

  comp_summary_file = file.path("cache", "%s-fit_real-comp_summary.%s")
  comp_sum_sum_file = file.path("cache", "%s-fit_real-comp_sum_sum.%s")
  frac_pos_file = file.path("cache", "%s-fit_real-frac_pos.%s")

  write.csv(df_comp_summary, file=sprintf(comp_summary_file, SESSION_ID, "csv"))
  print(
    xtable(df_comp_summary, digits=3),
    include.rownames=FALSE,
    file=sprintf(comp_summary_file, SESSION_ID, "tex")
  )  

  write.csv(df_comp_summary_summary, file=sprintf(comp_sum_sum_file, SESSION_ID, "csv"))
  print(
    xtable(df_comp_summary_summary, digits=3),
    include.rownames=FALSE,
    file=sprintf(comp_sum_sum_file, SESSION_ID, "tex")
  )
  
  write.csv(frac_pos_df, file=sprintf(frac_pos_file, SESSION_ID, "csv"))
  print(
    xtable(frac_pos_df, digits=3),
    include.rownames=FALSE,
    file=sprintf(frac_pos_file, SESSION_ID, "tex")
  )

}


# COMPARE BY PLOTS

p_comp_prop_no_thresh = df_comp %>%
  filter(thresh == Inf) %>%
  filter(grepl("(3k|empirical)", model)) %>%
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

# grid.arrange(p_comp_prop_no_thresh, p_comp_prop_thresh)


if (config$write) {
  p_comp_no_thresh_file = sprintf("%s-fit_real-p_comp_prop-no_thresh.png", SESSION_ID)
  ggsave(p_comp_prop_no_thresh, file=file.path("images", "xm_model", p_comp_no_thresh_file))
  p_comp_thresh_file = sprintf("%s-fit_real-p_comp_prop-thresh.png", SESSION_ID)
  ggsave(p_comp_prop_thresh, file=file.path("images", "xm_model", p_comp_thresh_file))
}
