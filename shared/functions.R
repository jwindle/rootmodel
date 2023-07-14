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


get_K <- function(r, mu_dx) {
  K = ceiling(r / mu_dx) - 1
  return(K)
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
    dat$K = with(dat, ceiling(r / mu_dx) - 1)
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


make_depths_plot <- function(depths_df, facet1, facet2) {

  primary_cols = c("depth", facet1, facet2)
  remaining_cols = setdiff(colnames(depths_df), primary_cols)

  alt_depth_df = depths_df["depth"]
  alt_depth_df[,"f1"] = depths_df[[facet1]]
  alt_depth_df[,"f2"] = depths_df[[facet2]]
  alt_depth_df[,"others"] = ""

  remain_name = paste(remaining_cols, collapse=",")
  n_remain = length(remaining_cols)
  
  if (n_remain == 1) {
    alt_depth_df[["others"]] = depths_df[,remaining_cols[1]]
  } else if (n_remain > 1) {
    alt_depth_df[["others"]] = apply(depths_df[,remaining_cols], 1, paste, collapse=",")
  }
  
  p_y_all = alt_depth_df %>%
    ggplot() +
    geom_histogram(
      aes(x=depth, y=..density.., fill=as.factor(others)),
      bins=10,
      # position="identity",
      alpha=1.0
    ) +
    facet_grid(col=vars(f1), row=vars(f2)) +
    xlab("depth (cm)") +
    ggtitle(sprintf("Distribution of depths %s, %s", facet1, facet2)) +
    guides(fill=guide_legend(title=remain_name)) +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(p_y_all)
}


posterior_predictive_depths <- function(model, par_array, base) {
  # Extract par_array here?
  ndim = length(dim(par_array))
  Y = apply(par_array, 1:(ndim-1), function(x){
    datx = c(base, as.list(x))
    out = sampling(model, dat=datx, algorithm="Fixed_param", iter=1, chains=1, refresh=0)
    return(extract(out, "y")[[1]])
  })
  return(Y)
}

posterior_predictive_paths <- function(model, par_array, base) {
  # Extract par_array here?
  ndim = length(dim(par_array))
  P = apply(par_array, 1:(ndim-1), simplify=TRUE, function(x){
    datx = c(base, as.list(x))
    out = sampling(model, dat=datx, algorithm="Fixed_param", iter=1, chains=1, refresh=0)
    paths = extract(out, "path")[[1]][1,,]
    return(paths)
  })
  P_dim = dim(P)
  if (with(base, N*(K+1) != P_dim[1])) {
    print("MISMATCH OF EXPECTED DIMENSIONS IN PATH")
  }
  P_array = array(as.numeric(P), dim=with(base, c(N, K+1, P_dim[-1])))
  return(P_array)
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

discretized_depth_distribution_2 <- function(paths, depth_breaks, epoch_levels) {
  df_list = list()
  dim_paths = dim(paths)
  n_steps = dim_paths[2]
  n_epochs = dim_paths[length(dim_paths)]
  # depths_first = paths[,1,,,drop=TRUE]
  depths_high = apply(paths, c(1, 3, 4), max)
  depths_final = paths[,n_steps,,,drop=TRUE]
  depths_high_matrix = matrix(as.numeric(depths_high), ncol=n_epochs, byrow=FALSE)
  depths_final_matrix = matrix(as.numeric(depths_final), ncol=n_epochs, byrow=FALSE)
  for (i in 1:n_epochs) {
    df_list[[i]] = data.frame(
      shallowest_cm = depths_high_matrix[,i],
      depth_cm = depths_final_matrix[,i]
    )
    df_list[[i]]$epoch = i
  }
  df = do.call(rbind, df_list)
  df$epoch = as.ordered(df$epoch)
  levels(df$epoch) = epoch_levels
  df$depth_bin = cut(df$depth_cm, breaks=depth_breaks, include.lowest=TRUE, ordered_result=TRUE)
  return(df)
}


discretized_depth_distribution_3 <- function(refined_paths, depth_breaks, epoch_levels) {
  df_list = list()
  paths = refined_paths
  dim_paths = dim(paths)
  n_epochs = dim_paths[3]
  # depths_first = paths[,1,,,drop=TRUE]
  depths_high = apply(paths, c(1, 3), max)
  depths_final = paths[,dim_paths[2],]
  depths_high_matrix = matrix(as.numeric(depths_high), ncol=n_epochs, byrow=FALSE)
  depths_final_matrix = matrix(as.numeric(depths_final), ncol=n_epochs, byrow=FALSE)
  for (i in 1:n_epochs) {
    df_list[[i]] = data.frame(
      shallowest_cm = depths_high_matrix[,i],
      depth_cm = depths_final_matrix[,i]
    )
    df_list[[i]]$epoch = i
  }
  df = do.call(rbind, df_list)
  df$epoch = as.ordered(df$epoch)
  levels(df$epoch) = epoch_levels
  df$depth_bin = cut(df$depth_cm, breaks=depth_breaks, include.lowest=TRUE, ordered_result=TRUE)
  return(df)
}


log_transform <- function(y, y0) {
  return(log(-(y - y0)))
}


inv_log_transform <- function(x, y0) {
  return(y0 - exp(x))
}


refine_paths <- function(paths, mu_dx, n, transform) {
  # Assuming we are on an equally spaced grid

  paths_t = aperm(paths, c(3, 1, 2, 4))
  paths_t_dim = dim(paths_t)
  n_epoch = paths_t_dim[4]
  Kp1 = paths_t_dim[3]

  EP_list = list()
  for (i in 1:n_epoch) {
    P = cbind(
      0.0,
      matrix(
        paths_t[,,,i],
        nrow=prod(paths_t_dim[1:2]),
        ncol=paths_t_dim[3],
        byrow=FALSE
      )
    )
    Y0 = P[,1:Kp1,drop=FALSE]
    DY = P[,2:(Kp1+1),drop=FALSE] - Y0
    one = rep(1, n)
    delta = mu_dx * seq(1, n, length.out=n) / n
    expanded_paths = aperm(outer(DY, delta) / mu_dx + outer(Y0, one), c(1, 3, 2))
    expanded_paths_dim = dim(expanded_paths)
    EP = matrix(expanded_paths, nrow=expanded_paths_dim[1], byrow=FALSE)
    EP_list[[i]] = EP
  }
  
  EP_array = simplify2array(EP_list)

  EP_array = EP_array + transform$offset
  if (transform$log_transform) {
    EP_array = -1. * exp(EP_array)
  }
  
  return(EP_array)
}


sim_plants <- function(mdl, N, p, mu_dm, sig_dm, shape_dm, time_grid, max_trials, r, mu_dx) {

  t_if_needed = function(x) {
    if (length(dim(x)) < 2) {
      return(t(x))
    }
    return(x)
  }

  K = get_K(r, mu_dx)
  ntimes = length(time_grid)
  
  p = t_if_needed(p)
  mu_dm = t_if_needed(mu_dm)
  sig_dm = t_if_needed(sig_dm)
  shape_dm = t_if_needed(shape_dm)

  n1 = nrow(p)
  n2 = nrow(mu_dm)

  # Two possible methods
  if (N == 0) {
    # Cross
    idc1 = rep(1:n1, times=n2)
    idc2 = rep(1:n2, n1)
  } else {
    # Sample with replacement
    idc1 = sample(1:n1, replace=TRUE, size=N)
    idc2 = sample(1:n2, replace=TRUE, size=N)
  }

  N = ifelse(N > 1, N, 1)

  # Now make samples
  df_list = list()
  for (i in 1:N) {
    
    i1 = idc1[i]
    i2 = idc2[i]
    
    num_roots = rbinom(ntimes, p=p[i1,], size=max_trials)
    root_times = rep(1:ntimes, times=num_roots)

    n_roots = length(root_times)

    # Check
    if (!all(table(root_times) == num_roots[num_roots > 0])) {
      print("PROBLEM!!!")
    }

    dat = list(
      mu_dx = mu_dx,
      K = K,
      r = r,
      n_times = ntimes,
      n_observed_times = n_roots,
      time_index = root_times,
      mu_dm = mu_dm[i2,],
      sig_dm = sig_dm[i2,],
      shape_dm = shape_dm[i2,]
    )

    sim_plant_out = sampling(
      sim_plant_mdl,
      data=dat,
      iter=1,
      warmup=0,
      chains=1,
      cores=1,
      algorithm="Fixed_param",
      refresh=FALSE
    )

    root_df = data.frame(
      time_index = root_times,
      plant = i,
      p_id = i1,
      d_id = i2,
      root_id = 1:n_roots,
      angle = runif(n_roots, min=0, max=2 * pi),
      starting_hour = sample(0:23, size=n_roots, replace=TRUE)
    )

    paths_df = as.data.frame(drop(extract(sim_plant_out, "path")[[1]]))
    colnames(paths_df) = paste("depth", 1:(K+1), sep="_")

    df_list[[i]] = cbind(root_df, paths_df)
    
  }

  df = do.call(rbind, df_list)
  return(df)
}
