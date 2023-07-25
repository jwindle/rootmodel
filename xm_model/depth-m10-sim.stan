/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
#include "common.stan"
data {
	int<lower=1> N;
	real<lower=0> r;
	real <lower=0> mu_dx;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	real mu_dm;
	real<lower=0> sig_dm;
	real shape_dm;
}
transformed data {
	vector[K+1] x_left = x_left_from_mu_dx(r, mu_dx, K);
	real delta = shape_dm / sqrt(1. + square(shape_dm));
	real sqrt_2_over_pi = sqrt(2. / pi());
	real ep_mean = sqrt_2_over_pi * delta;
	real ep_std = sqrt(1 - square(sqrt_2_over_pi * delta));
	matrix[K+1,K+1] PT = path_transform(r, mu_dx, K);
}
generated quantities {
	vector[K+1] dm[N];
	vector[K+1] ep[N];
	real z[N];
	real y[N];
	vector[K+1] path[N];
	/* vector[K+1] m[N]; */
	/* vector[K+1] theta[N]; */
	/* vector[K] dl[N]; */
	for (i in 1:N) {
		for (j in 1:(K+1)) {
			ep[i][j] = skew_normal_rng(0.0, 1.0, shape_dm);
		}
		dm[i] = (sig_dm * (ep[i] - ep_mean)) / ep_std + mu_dm;
		z[i] = dot_product(dm[i], x_left);
		y[i] = normal_rng(z[i], 0.01);
		/* m[i][1] = dm_0[i]; */
		/* for (j in 1:K) { */
		/* 	m[i][j+1] = m[i][j] + dm[i][j]; */
		/* } */
		/* theta[i] = atan(m[i]); */
		/* dl[i] = dx ./ cos(theta[i][1:K]); */
		path[i] = PT * dm[i];
	}
}
