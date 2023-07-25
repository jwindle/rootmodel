/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
#include "common.stan"
data {
	real<lower=0.> mu_dx;        /* dx step size */
	int<lower=1> K;              /* number of kinks                    */
	real<lower=0> r;             /* radius at which a root is detected */
	int<lower=1> n_times;
	int<lower=1> n_observed_times;
	int<lower=1, upper=n_times> time_index[n_observed_times];
	vector[n_times] mu_dm;
	vector<lower=0.>[n_times] sig_dm;
	vector[n_times] shape_dm;
}
transformed data {
	real sqrt_2_over_pi = sqrt(2. / pi());
	vector[K+1] x_left = x_left_from_mu_dx(r, mu_dx, K);
	vector[n_times] delta_dm = shape_dm ./ sqrt(1. + square(shape_dm));
	vector[n_times] b_delta_dm = sqrt_2_over_pi * delta_dm;
	vector<lower=0.>[n_times] omega_dm = sig_dm ./ sqrt(1. - square(b_delta_dm));
	vector[n_times] xi_dm = mu_dm - b_delta_dm .* omega_dm;
	matrix[K+1,K+1] PT = path_transform(r, mu_dx, K);
}
generated quantities {
	vector[K+1] dm[n_observed_times];
	/* real z[n_observed_times]; */
	real y[n_observed_times];
	vector[K+1] path[n_observed_times];
	for (i in 1:n_observed_times) {
		int t = time_index[i];
		for (j in 1:(K+1)) {
			dm[i][j] = skew_normal_rng(xi_dm[t], omega_dm[t], shape_dm[t]);
		}
		y[i] = dot_product(dm[i], x_left);
		/* y[i] = normal_rng(z[i], 0.01); */
		path[i] = PT * dm[i];
	}
}
