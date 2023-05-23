# -*- ess-style: RStudio; -*-

library("readr")
library("ggplot2")
library("dplyr")
library("latex2exp")
library("gridExtra")
library("rstan")



load_rt_data_old <- function(filename) {
    col_types = cols(
        electrode_bins = col_integer(),
        paddle_bins = col_character(),
        time_bins = col_date(),
        device_id = col_character(),
        count = col_double(),
        uptime = col_double(),
        x = col_double(),
        y = col_double(),
        crop = col_character(),
        rep = col_integer(),
        bench = col_integer()
    )
    df = read_csv(filename, col_types = col_types)
    df$days_since_start = as.numeric(with(df, time_bins - min(time_bins)))
    return(df)
}


load_rt_data <- function(filename) {
    col_types = cols(
        electrode_pair = col_integer(),
        paddle = col_character(),
        date = col_date(),
        device_id = col_character(),
        count = col_double(),
        uptime = col_double(),
        depth_cm = col_double(),
        x = col_double(),
        y = col_double(),
        crop = col_character(),
        rep = col_integer(),
        bench = col_integer()
    )
    df = read_csv(filename, col_types = col_types)
    df$days_since_start = as.numeric(with(df, date - min(date)))
    return(df)
}


coarsen_rt_data <- function(df, ...) {
    args = list(...)
    col_names = colnames(df)
    for (nm in names(args)) {
        if (nm %in% col_names) {
            df[[nm]] = cut(df[[nm]], breaks=args[[nm]], include.lowest=TRUE, ordered_result=TRUE)
        } else {
            cat("Name", nm, "not found.\n")
        }
    }
    return(df)
}


agg_rt_data <- function(df_all, ...) {
    df_new = df_all %>%
        group_by(...=...) %>%
        summarize(
                  count=sum(count),
                  uptime=sum(uptime)
                  ) %>%
        mutate(rate = count / uptime)
    return(df_new)
}


coarsen_and_agg_rt_data <- function(df, ...) {

}



reconstruct_path <- function(X, DM_0, DM, grid) {
    ngrid = length(grid)
    X_dim = dim(X)
    path = array(0.0, dim=c(X_dim[1], X_dim[2], ngrid))
    for (i in 1:ngrid) {
        grid_r = grid[i]
        X_left = ifelse(grid_r > X, grid_r - X, 0.)
        Z = apply(DM * X_left, c(1, 2), sum) + DM_0 * grid_r
        path[,,i] = Z
    }
    return(path)
}


log_normal_param <- function(m, s) {
  # Mean and std dev of positive random variable to parameters for
  # log-normal distribution with the same moments.
  v = s**2
  ln_sig2 = log(v / m^2 + 1)
  ln_mu = log(m) - 0.5 * ln_sig2
  ln_sig = sqrt(ln_sig2)
  return(c(ln_mu, ln_sig))
}


file_helper <- function(dir, base, suffix) {
  file_name = paste(base, suffix, sep="")
  path = file.path(dir, file_name)
  return(path)
}


downsample_depths <- function(depths, use=0.1) {
  # We want to use something non-random here to minimize sampling
  # variance and ensure a replicable result.  We return the index so
  # that we can extract paths too, if we want.
  nd = length(depths)
  df = data.frame(depth=depths, idc=1:nd)
  N = ceiling(nd * use)
  idc_to_use = as.integer(round(seq(1, nd, length.out=N), 0))
  df_sort = df[order(df$depth),]
  return(df_sort[idc_to_use,])
}

x_left_for_fixed_dx <- function(r, mu_dx) {
  k = ceiling(r / mu_dx) - 1
  dx = rep(mu_dx, k+1)
  dx[1] = 0.0
  x = cumsum(dx)
  x_left = r - x
  return(x_left)
}


depth_fit_wrapper <- function(filename) {


}


make_base_data <- function(param_config) {
  base = list()
  for (nm in names(param_config)) {
    val = param_config[[nm]]
    if (is.numeric(val) && length(val) == 1) {
      base[[nm]] = val
    }
  }
  return(base)
}


make_parameter_df <- function(param_config) {
  grid = list()
  for (nm in names(param_config)) {
    val = param_config[[nm]]
    if (is.numeric(val) && length(val) > 1) {
      grid[[nm]] = val
    } else if (is.list(val)) {
      grid[[nm]] = with(param_config[[nm]], seq(from, to, length.out=length.out))
    }
  }
  grid_df = expand.grid(grid)
  return(grid_df)
}


simulate_over_parameter_df <- function(model, param_df, base_data) {
  nr = nrow(param_df)
  param_names = colnames(param_df)
  samp_list = list()
  for (i in 1:nr) {
    dat = c(base_data, as.list(param_df[i,]))
    samp_list[[i]] = sampling(model, data=dat, chains=1, iter=1, algorithm="Fixed_param")
  }
  return(samp_list)
}


extract_depths_list <- function(samp_list) {
  depths_list = lapply(samp_list, function(x){drop(extract(x, "y")[[1]])})
  return(depths_list)
}


summarize_depth_distributions <- function(depths_list, param_df) {
  summary_df = cbind(param_df, as.data.frame(t(sapply(depths_list, summary))))
  return(summary_df)
}


make_depths_df <- function(depths_list, param_df) {
  df_components = list()
  param_names = names(param_df)
  for (i in 1:length(depths_list)) {
    df_components[[i]] = data.frame(depth = depths_list[[i]])
    for (nm in param_names) {
      df_components[[i]][[nm]] = param_df[[nm]][i]
    }
  }
  depths_df = do.call(rbind, df_components)
  return(depths_df)
}


posterior_predictive <- function(model, par_array, base) {
  # Extract par_array here?
  ndim = length(dim(par_array))
  Y = apply(par_array, 1:(ndim-1), function(x){
    datx = c(base, as.list(x))
    out = sampling(model, dat=datx, algorithm="Fixed_param", iter=1, chains=1)
    return(extract(out, "y")[[1]])
  })
  return(Y)
}


discretized_depth_distribution <- function(depths, depth_breaks, epoch_levels) {
  df_list = list()
  dim_depths = dim(depths)
  n_epochs = dim_depths[length(dim_depths)]
  depths_matrix = matrix(depths, ncol=n_epochs, byrow=FALSE)
  for (i in 1:n_epochs) {
    df_list[[i]] = data.frame(depth_cm = depths_matrix[,i])
    df_list[[i]]$epoch = i
  }
  df = do.call(rbind, df_list)
  df$epoch = as.ordered(df$epoch)
  levels(df$epoch) = epoch_levels
  df$depth_bin = cut(df$depth_cm, breaks=depth_breaks, include.lowest=TRUE, ordered_result=TRUE)
  return(df)
}
