# -*- ess-style: RStudio; -*-

source("shared/functions.R")


library("rstan")
library("bayesplot")
library("yaml")


# CONFIG
config = read_yaml("xm_model/depth-MAIN-fit_real-log-config.yaml")

# EPOCHS = config$epochs
MODELS_TO_RUN = config$models_to_run
# MODELS_TO_RUN = "m10_1k"
MODELS_TO_RUN

USE = config$data_frac
CROP = config$crop
c(USE, CROP)

SESSION_ID = "TEST"

# Load and take subset of data
# df_all = read.csv("data/acc-1-v3-cleaned.csv", stringsAsFactors=TRUE)
load(file.path("cache", "detections.RData"))

head(df_all)

# Sample (for testing routines)
depths_df = df_all %>%
  filter(count > 0, crop == CROP) %>%
  slice_sample(prop=USE, replace=FALSE)

## # Or sample from epoch
## depths_df = df_all %>%
##   filter(count > 0, epoch_group == 2)

# Get data for Stan
y_data = depths_df$depth_cm
group = depths_df$epoch_group
N = length(y_data)
M = length(levels(depths_df$epoch))


# Run models
samp_list = list()

for (model_name in MODELS_TO_RUN) {

  model_config = config$models[[model_name]]
  
  prior_info = model_config$prior_info
  if ("mu_dx" %in% names(prior_info)) {
    prior_info$K = with(prior_info, get_K(r, mu_dx))
  }
  if (!("offset" %in% names(model_config))) {
    model_config$offset = 0.
  }

  # Stan files
  stan_fit_file = sprintf("depth-%s-fit.stan", model_config$stan_file_id)
  
  # Param & pars of interest
  dat_i = list(
    y = y_data,
    N = length(y_data),
    M = M,
    group = group
  )
  dat = c(dat_i, prior_info)
  # datenv = list2env(dat)

  # Transforms, if needed
  if (model_config$log_transform) {
    dat$y = log(-dat$y)
  }
  dat$y = dat$y - model_config$offset

  # Simulate
  mcmc_config = model_config$mcmc
  pars = model_config$keep_pars

  samp = stan(
    file.path("xm_model", stan_fit_file),
    data = dat,
    chains = mcmc_config$chains,
    cores = mcmc_config$cores,
    iter = mcmc_config$iter,
    warmup = mcmc_config$warmup,
    init_r = 1,
    control = list(max_treedepth=10, metric="diag_e"),
    verbose = TRUE,
    model_name = sprintf("depth-fit-%s", model_name),
    pars = pars
  )

  samp_copy = samp
  # samp_list[[model_name]] = samp

  # Summary
  print(summary(samp, pars))

  # Plot
  samp_array = extract(samp, pars, permute=FALSE)

  p_hist = mcmc_hist(samp_array) + ggtitle(sprintf("Posteriors model %s", model_name))
  p_hist

  if (config$write) {
    p_hist_file = sprintf("%s-%s-real-p_hist.png", SESSION_ID, model_name)
    ggsave(p_hist, file=file.path("images", "xm_model", p_hist_file))
  }

  # p_pairs = mcmc_pairs(samp_array, grid_args = list(top=sprintf("Pairs model %s", model_name)))
  # p_pairs
  
  ## if (config$write) {
  ##   p_pairs_file = sprintf("%s-%s-real-p_pairs.png", SESSION_ID, model_name)
  ##   ggsave(p_pairs, file=file.path("..", "images", "xm_model", p_pairs_file))
  ## }

  # # Stan versions of plotting routines
  # pairs(samp, pars=pars)
  # traceplot(samp, pars=pars)

  if (config$write) {
    
    samp_file = sprintf("%s-fit.RData", model_name)
    save(samp, file=file.path("cache", samp_file))
    
  }

}

