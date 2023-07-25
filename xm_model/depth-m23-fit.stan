/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
#include "common.stan"
data {
	int<lower=1> n_obs;
	int<lower=1> n_times;
	real time_grid[n_times];
	int<lower=1, upper=n_times> time_group[n_obs];
	real<lower=0.> r;
	real<lower=0.> mu_dx;
	int<lower=1> K;
	real y[n_obs];
	real mu_y_prior_mean;
	real mu_y_prior_std;
	real sig_y_prior_mean;
	real sig_y_prior_std;
	real shape_dm_prior_mean;
	real shape_dm_prior_std;
}
transformed data {
	vector[K+1] x_left = x_left_from_mu_dx(r, mu_dx, K);
	real x_left_l1 = sum(x_left);
	real x_left_l2 = l2_norm(x_left);
	real mu_dm_prior_mean = mu_y_prior_mean / x_left_l1;
	real mu_dm_prior_std = mu_y_prior_std / x_left_l1;
	real sig_dm_prior_mean = sig_y_prior_mean / x_left_l2;
	real sig_dm_prior_std = sig_y_prior_std / x_left_l2;
	real sqrt_2_over_pi = sqrt(2. / pi());
	real tiny = 1e-9; // Try to ensure positive definite
	real sqrt_tiny = sqrt(tiny);
}
parameters {
	/* Gaussian process */
	real<lower=0.> gp_ell[3];
	real<lower=0.> gp_sig[3];
	vector[n_times] eta[3];
	/* Global parameters values */
	real mu_y_gbl_mean;
	real<lower=0.> sig_y_gbl_mean;
	real shape_dm_gbl_mean;
	/* Changes in slope */
	vector[K+1] dm[n_obs];
}
transformed parameters {
	/* Gaussian process */
	vector[n_times] f_gp[3];
	for (i in 1:3) {
		matrix[n_times, n_times] L;
		matrix[n_times, n_times] C = gp_exp_quad_cov(time_grid, gp_sig[i], gp_ell[i]);
		// diagonal elements - this may not be necessary since I
		// control the design, i.e. `time_grid`.
		for (n in 1:n_times) {
			C[n, n] = C[n, n] + tiny;
		}
		L = cholesky_decompose(C);
		f_gp[i] = L * eta[i];
	}
	/* Transformation to location, scale, shape */
	vector[n_times] mu_dm = (f_gp[1] + mu_y_gbl_mean) / x_left_l1;
	vector<lower=0.>[n_times] sig_dm = (exp(f_gp[2]-0.5*square(gp_sig[2])) * sig_y_gbl_mean) / x_left_l2;
	vector[n_times] shape_dm = f_gp[3] + shape_dm_gbl_mean;
	vector[n_times] delta_dm = shape_dm ./ sqrt(1. + square(shape_dm));
	vector[n_times] b_delta_dm = sqrt_2_over_pi * delta_dm;
	vector<lower=0.>[n_times] omega_dm = sig_dm ./ sqrt(1. - square(b_delta_dm));
	vector[n_times] xi_dm = mu_dm - b_delta_dm .* omega_dm;
	/* Transformation to depths */
	real z[n_obs];
	for (i in 1:n_obs) {
		z[i] = dot_product(dm[i], x_left);
	}
}
model {
	gp_ell ~ inv_gamma(5, 5);
	gp_sig ~ std_normal();
	mu_y_gbl_mean ~ normal(mu_y_prior_mean, mu_y_prior_std);
	sig_y_gbl_mean ~ normal(sig_y_prior_mean, sig_y_prior_std) T[0,];
	shape_dm_gbl_mean ~ normal(shape_dm_prior_mean, shape_dm_prior_std);
	for (i in 1:3) {
		eta[i] ~ std_normal();
	}
	for (i in 1:n_obs) {
		int t = time_group[i];
		dm[i] ~ skew_normal(xi_dm[t], omega_dm[t], shape_dm[t]);
		y[i] ~ normal(z[i], 0.01);
	}
}
generated quantities {
	real mu_dm_gbl_mean = mu_y_gbl_mean / x_left_l1;
	real<lower=0.> sig_dm_gbl_mean = sig_y_gbl_mean / x_left_l2;
	vector[n_times] mu_y = mu_dm * x_left_l1;
	vector<lower=0.>[n_times] sig_y = sig_dm * x_left_l2;
}
