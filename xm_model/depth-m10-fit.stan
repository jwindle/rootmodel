/* -*- c-basic-offset: 4; -*- */
#include "common.stan"
data {
	int<lower=1> N;              /* number of total observations       */
	int<lower=1> M;              /* number of groups                   */
	real<lower=0> mu_dx;         /* dx step size */
	int<lower=1> K;              /* number of kinks                    */
	real<lower=0> r;             /* radius at which a root is detected */
	int<lower=1, upper=M> group[N];
	real mu_y_prior_mean;
	real mu_y_prior_std;
	real sig_y_prior_mean;
	real sig_y_prior_std;
	real shape_dm_prior_mean;
	real shape_dm_prior_std;
	real y[N];
}
transformed data {
	vector[K+1] x_left = x_left_from_mu_dx(r, mu_dx, K);
	real x_left_l1 = sum(x_left);
	real x_left_l2 = l2_norm(x_left);
	real mu_dm_prior_mean = mu_y_prior_mean / x_left_l1;
	real mu_dm_prior_std = mu_y_prior_std / x_left_l1;
	real sig_dm_prior_mean = sig_y_prior_mean / x_left_l2;
	real sig_dm_prior_std = sig_y_prior_std / x_left_l2;
	/* real sig_dm_alpha = square(sig_dm_prior_mean / sig_dm_prior_std) + 2; */
	/* real sig_dm_beta = sig_dm_prior_mean * (sig_dm_alpha - 1); */
	real sqrt_2_over_pi = sqrt(2. / pi());
}
parameters {
	real mu_dm[M];
	real<lower=0> sig_dm[M];
	real shape_dm[M];
	vector[K+1] dm[N];
}
transformed parameters {
	real z[N];
	real delta[M];
	real b_delta[M];
	real xi_dm[M];
	real eta_dm[M];
	for (j in 1:M) {
		delta[j] = shape_dm[j] / sqrt(1. + square(shape_dm[j]));
		b_delta[j] = sqrt_2_over_pi * delta[j];
		eta_dm[j] = sig_dm[j] / sqrt(1 - square(b_delta[j]));
		xi_dm[j] = mu_dm[j] - eta_dm[j] * b_delta[j];
	}
	for (i in 1:N) {
		z[i] = dot_product(dm[i], x_left);
	}
}
model {
	for (i in 1:M) {
		mu_dm[i] ~ normal(mu_dm_prior_mean, mu_dm_prior_std);
		sig_dm[i] ~ normal(sig_dm_prior_mean, sig_dm_prior_std) T[0,];
		/* sig_dm[i] ~ inv_gamma(sig_dm_alpha, sig_dm_beta); */
		shape_dm[i] ~ normal(shape_dm_prior_mean, shape_dm_prior_std);
	}
	for (i in 1:N) {
		int j = group[i];
		dm[i] ~ skew_normal(xi_dm[j], eta_dm[j], shape_dm[j]);
		y[i] ~ normal(z[i], 0.01);
	}
}
generated quantities {
	real mu_y[M];
	real sig_y[M];
	for (j in 1:M) {
		mu_y[j] = mu_dm[j] * x_left_l1;
		sig_y[j] = sig_dm[j] * x_left_l2;
	}
}
