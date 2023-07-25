/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
data {
	int<lower=1> N;              /* number of total observations       */
	int<lower=1> M;              /* number of groups                   */
	real<lower=0> mu_dx;         /* dx step size */
	int<lower=1> K;              /* number of kinks                    */
	int<lower=1, upper=M> group[N]; /* group membership */
	real<lower=0> r;             /* radius at which a root is detected */
	real mu_dm_prior_mean;
	real mu_dm_prior_std;
	real sig_dm_prior_mean;
	real sig_dm_prior_std;
	real shape_dm_prior_mean;
	real shape_dm_prior_std;
	real y[N];
}
transformed data {
	vector[K] rvec = rep_vector(r, K);
	vector[K] dx = rep_vector(mu_dx, K);
	vector[K] x = cumulative_sum(dx);
	vector[K] x_left = fdim(rvec, x);
	real sig_dm_alpha = square(sig_dm_prior_mean / sig_dm_prior_std) + 2;
	real sig_dm_beta = sig_dm_prior_mean * (sig_dm_alpha - 1);
	real sqrt_2_over_pi = sqrt(2. / pi());
}
parameters {
	real mu_dm[M];
	real<lower=0> sig_dm[M];
	real shape_dm[M];
	/* real<lower=-1., upper=1.> delta[M]; */
	real ep_0[N];
	vector[K] ep[N];
}
transformed parameters {
	real dm_0[N];
	vector[K] dm[N];
	real z[N];
	real<lower=-1., upper=1.> delta[M];
	/* real shape_dm[M]; */
	real ep_mean[M];
	real<lower=0.> ep_std[M];
	for (j in 1:M) {
		delta[j] = shape_dm[j] / sqrt(1. + square(shape_dm[j]));
		/* shape_dm[j] = delta[j] / sqrt(1 - square(delta[j])); */
		ep_mean[j] = sqrt_2_over_pi * delta[j];
		ep_std[j] = sqrt(1. - square(sqrt_2_over_pi * delta[j]));
	}
	for (i in 1:N) {
		int j = group[i];
		dm_0[i] = mu_dm[j] + sig_dm[j] * (ep_0[i] - ep_mean[j]) / ep_std[j];
		dm[i] = mu_dm[j] + sig_dm[j] * (ep[i] - ep_mean[j]) / ep_std[j];
		z[i] = dm_0[i] * r + sum(dm[i] .* x_left);
	}
}
model {
	for (j in 1:M) {
		mu_dm[j] ~ normal(mu_dm_prior_mean, mu_dm_prior_std);
		sig_dm[j] ~ normal(sig_dm_prior_mean, sig_dm_prior_std) T[0,];
		/* sig_dm[j] ~ inv_gamma(sig_dm_alpha, sig_dm_beta); */
		shape_dm[j] ~ normal(shape_dm_prior_mean, shape_dm_prior_std);
		/* delta[j] ~ uniform(-1., 1.); */
	}
	for (i in 1:N) {
		int j = group[i];
		ep_0[i] ~ skew_normal(0.0, 1.0, shape_dm[j]);
		ep[i] ~ skew_normal(0.0, 1.0, shape_dm[j]);
		y[i] ~ normal(z[i], 0.01);
	}
}
