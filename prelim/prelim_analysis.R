# Jesse Windle, 2023
# jesse@bayesfactor.net


library("ggplot2")
library("readr")
library("dplyr")

# library("lme4")
library("lmerTest")
library("kernlab")

source("../shared/functions.R")

options(pillar.subtle = FALSE)

# Load the data.  We have likely aggregated the data before hand to 24 hour
# periods, for instance.  To deal with missing data, we may also have tried to
# identify hours for which we had sufficiently large uptime across all paddles
# and electrodes (e.g. 10 and above) and then redeclared all those paddles and
# uptime to be 12.  In that way we have a complete observation of 12 x 22
# paddles for an entire hour.  

df_all = load_rt_data("../data/acc-1-v3.csv") %>%
    filter(uptime > 0, electrode_pair > 1, days_since_start > 4)

df_all$old_count = df_all$count
df_all$count = 1 * (df_all$old_count > 0)
df_all$rate = with(df_all, count / uptime)

print(df_all, n=30)

sapply(df_all, summary)


# Crop by locations
df_all %>%
    select(x, y, crop) %>%
    distinct() %>%
    ggplot() +
        geom_tile(aes(x, y, fill=as.factor(crop)))


# For data that we have preprocessed, we should have everything be a multiple of
# 2*12.  The preprocessing though, ensures that we have (most) all paddles and
# electrodes within a 5 minute interval.
mean(df_all$uptime %% (2*12) == 0)



depths = -1 * (seq(0, 21, 1) * (16.2/21) + 0.5)
# depths
width = abs(mean(diff(depths)))
depth_ivals = c(depths[1] + 0.5 * width, depths - 0.5 * width)
# depth_ivals
pairwise_depths = 0.5 * (depths[seq(1, 21, 2)] + depths[seq(2, 22, 2)])
# pairwise_depths


depth_bins_coarse = depth_ivals[seq(3, 23, by=4)]
df_coarse = coarsen_rt_data(
    df_all,
    days_since_start = seq(0, 28, 4),
    depth_cm = depth_bins_coarse
)


## # Coarsen grid on which we have RT data
## df_coarse = coarsen_rt_data(
##     df_all,
##     days_since_start = seq(0, 28, 2),
##     electrode_pair = seq(1, 11, 1)
## )



# Aggregate over paddles
df_dev_time_elec = agg_rt_data(df_coarse, crop, device_id, rep, days_since_start, depth_cm)

head(df_dev_time_elec, n=20)




# Aggregate over paddles and devices
df_time_elec = agg_rt_data(df_coarse, crop, rep, days_since_start, depth_cm)

head(df_time_elec)


# Plot
df_time_elec %>%
    filter(
        crop %in% c("Corn", "Wheat"),
    ) %>%
    ggplot() +
    geom_tile(aes(days_since_start, depth_cm, fill=rate)) + 
    scale_fill_gradient2() +
    # scale_fill_gradientn(colours = terrain.colors(10)) +
    facet_wrap(rep ~ crop)


# Aggregate over paddles and electrodes
df_dev_time = agg_rt_data(df_coarse, crop, device_id, days_since_start)

tab = with(subset(df_dev_time, crop=="Corn"), table(count))
tab

plot(log(tab))

tab_df = as.data.frame(tab)
tab_df$count = as.numeric(tab_df$count)

tab_lm = lm(log(Freq) ~ count, data=tab_df)
summary(tab_lm)

df_dev_time %>%
    filter(crop == "Corn") %>%
    ggplot() +
    geom_boxplot(aes(as.factor(days_since_start), count))

df_dev_time %>%
    filter(crop == "Corn") %>%
    ggplot() +
    geom_boxplot(aes(as.factor(days_since_start), uptime))


df_dev_time %>%
  filter(crop == "Corn") %>%
  ggplot() +
  geom_bar(aes(count)) +
  facet_wrap(~ days_since_start)



# Just time
df_time = agg_rt_data(df_all, crop, days_since_start)

df_time %>%
    filter(crop == "Corn") %>%
    ggplot() +
    geom_point(aes(days_since_start, log(rate)))



# Just device
df_dev = agg_rt_data(df_coarse, crop, device_id)

df_dev %>%
  filter(crop == "Corn") %>%
  ggplot() +
  geom_histogram(aes(count), bins=10)


# Device a 2 times
df_very_coarse = coarsen_rt_data(
    df_all,
    days_since_start = seq(0, 28, 14),
    electrode_pair = seq(1, 11, 1)
)


df_dev_2t = agg_rt_data(df_very_coarse, crop, device_id, days_since_start)

head(df_dev_2t)

df_dev_2t %>%
  filter(crop == "Corn") %>%
  ggplot() +
  geom_histogram(aes(count), bins=10) +
  facet_wrap(~ days_since_start)



# All
df_dev_day = agg_rt_data(df_all, crop, device_id, days_since_start)

head(df_dev_day)


with(df_dev_day, table(count, days_since_start))




# Pick crop and fit
CROP = "Corn"

## Fit mixed model (no smoothing over time / elec) on finest scale

df_all_one_crop = df_all %>%
    filter(crop == CROP, uptime > 0)

mm1 = lmer(rate ~ (1|days_since_start:electrode_pair), data=df_all_one_crop)

summary(mm1)



## Fit mixed model (smoothing by coarsening grid) on coarser scale

df_coarse_one_crop = df_coarse %>%
    filter(crop == CROP, uptime > 0)

head(df_coarse_one_crop)

mm2 = lmer(rate ~ (1|days_since_start:electrode_pair), data=df_coarse_one_crop)

summary(mm2)

## Get random effects for plotting
tmp = df_coarse %>% select(days_since_start, electrode_pair)
df_coarse_fit = tmp[!duplicated(tmp),]

df_coarse_fit$fit = ranef(mm2)[[1]][[1]]

head(df_coarse_fit)

### Plot coarser fit
ggplot(df_coarse_fit) + 
    geom_tile(aes(days_since_start, electrode_pair, fill=fit)) +
    scale_fill_gradient2()


## Fit model on coarser scale after aggregating over paddle

df_dte_1 = df_dev_time_elec %>%
    filter(crop == CROP, uptime > 0) %>%
    ungroup()

df_dte_1$days_fact_num = as.numeric(df_dte_1$days_since_start) 
df_dte_1$elec_fact_num = as.numeric(df_dte_1$electrode_pair)

head(df_dte_1)

mm3 = lmer(rate ~ (1|days_since_start:electrode_pair), data=df_dte_1)

summary(mm3)

tmp = df_dte_1 %>% select(days_since_start, electrode_pair, days_fact_num, elec_fact_num)
df_agg_fit = tmp[!duplicated(tmp),]

df_agg_fit

df_agg_fit$fit = ranef(mm3)[[1]][[1]]


### Plot agg fit

ggplot(df_agg_fit) + 
    geom_tile(aes(days_since_start, electrode_pair, fill=fit)) +
    scale_fill_gradient2()



# Spatial mixed model

# This did not work well with Matern kernel
# sp1 = fitme(rate ~ 1 + Matern(1|days_fact_num + elec_fact_num), data=df_dte_1)


# Gaussian process

my = mean(df_dte_1$rate)
ygp = df_dte_1$rate - my

gp = gausspr(x = ~ days_fact_num + elec_fact_num, y = ygp, data = df_dte_1)

gp

df_agg_fit$gp_fit = predict(gp, df_agg_fit)

ggplot(df_agg_fit) + 
    geom_tile(aes(days_since_start, electrode_pair, fill=gp_fit)) +
    scale_fill_gradient2()

gp2 = gausspr(x = ~ days_fact_num + elec_fact_num, y = ygp, data = df_dte_1, kpar = list(sigma=0.2))

gp2

df_agg_fit$gp2_fit = predict(gp2, df_agg_fit)

ggplot(df_agg_fit) + 
    geom_tile(aes(days_since_start, electrode_pair, fill=gp2_fit)) +
    scale_fill_gradient2()



# Plot distributions by depth
ggplot(df_agg_fit) +
  geom_line(aes(elec_fact_num, gp2_fit)) +
  facet_wrap(~ days_fact_num)




# Plot distribution by depth

# Coarsen grid on which we have RT data

for_depth_plots = coarsen_rt_data(
    df_all,
    days_since_start = c(0, 8, 16, 22, 28),
    electrode_pair = seq(1, 11, 1)
)


# Aggregate over paddles
df_depth_plots = agg_rt_data(for_depth_plots, crop, device_id, rep, days_since_start, electrode_pair)


head(df_depth_plots)


df_depth_plots %>%
    filter(crop == "Corn") %>%
    ggplot() +
    geom_boxplot(aes(electrode_pair, rate)) +
    facet_wrap(~ days_since_start)


df_depth_plots %>%
    filter(crop == "Corn") %>%
    ggplot() +
    geom_boxplot(aes(electrode_pair, count)) +
    facet_wrap(~ days_since_start)




