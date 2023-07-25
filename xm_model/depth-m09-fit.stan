/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
#include "common.stan"
data {
	int<lower=1> N;
	int<lower=1> M;
	real<lower=0> r;
	real <upper=0.0> min_depth;
	real <upper=0.0> max_depth;
	real tilde_mu_prior_mean;
	real tilde_mu_prior_std;
	real tilde_sig_prior_mean;
	real tilde_sig_prior_std;
	real <lower=0.0, upper=1.0> frac_var;
	int<lower=1, upper=M> group[N]; /* group membership */
	real y[N];
	/* /\* Below is only necessary for going back to mu_dm and sig_dm *\/ */
	/* real <lower=0> mu_dx; */
	/* int<lower=1> K; /\* ceil(r / mu_dx) - 1; does not work in stan *\/ */
}
transformed data {
	real lb = min_depth;
	real ub = max_depth;
	/* /\* Below is only necessary for going back to mu_dm and sig_dm *\/ */
	/* vector[K+1] x_left = x_left_from_mu_dx(r, mu_dx, K); */
	/* real x_left_one = sum(x_left); */
	/* real x_left_norm = sqrt(sum(square(x_left))); */
}
parameters {
	real tilde_mu[M];
	real<lower=0> tilde_sig[M];
	real z[N];
}
transformed parameters {
	real sig_z[M];
	real sig_y[M];
	for (i in 1:M) {
		sig_z[i] = frac_var * tilde_sig[i];
		sig_y[i] = (1 - frac_var) * tilde_sig[i];
	}
}
model {
	for (i in 1:M) {
		tilde_mu[i] ~ normal(tilde_mu_prior_mean, tilde_mu_prior_std);
		tilde_sig[i] ~ normal(tilde_sig_prior_mean, tilde_sig_prior_std) T[0,];
	}
	for (i in 1:N) {
		int g_i = group[i];
		z[i] ~ normal(tilde_mu[g_i], sig_z[g_i]);
		y[i] ~ normal(z[i], sig_y[g_i]) T[lb, ub];
	}
}
generated quantities {
	/* real mu[M]; */
	/* real sig1[M]; */
	/* real sig2[M]; */
	/* for (i in 1:M) { */
	/* 	mu[i] = tilde_mu[i] / x_left_one; */
	/* 	sig1[i] = sqrt(frac_var) * tilde_sig[i] / x_left_one; */
	/* 	sig2[i] = sqrt(1 - frac_var) * tilde_sig[i] / x_left_norm; */
	/* } */
}
