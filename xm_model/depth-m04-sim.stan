/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
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
	vector[K] rvec = rep_vector(r, K);
	vector[K] dx = rep_vector(mu_dx, K);
	vector[K] x = cumulative_sum(dx);
	vector[K] x_left = fdim(rvec, x);
	real delta = shape_dm / sqrt(1. + square(shape_dm));
	real sqrt_2_over_pi = sqrt(2. / pi());
	real ep_mean = sqrt_2_over_pi * delta;
	real ep_std = sqrt(1 - square(sqrt_2_over_pi * delta));
}
generated quantities {
	real dm_0[N];
	real ep_0[N];
	vector[K] dm[N];
	vector[K] ep[N];
	real z[N];
	real y[N];
	vector[K+1] m[N];
	vector[K+1] theta[N];
	vector[K] dl[N];
	for (i in 1:N) {
		ep_0[i] = skew_normal_rng(0.0, 1.0, shape_dm);
		dm_0[i] = mu_dm + sig_dm * (ep_0[i] - ep_mean) / ep_std;
		for (j in 1:K) {
			ep[i][j] = skew_normal_rng(0.0, 1.0, shape_dm);
			dm[i][j] = mu_dm + sig_dm * (ep[i][j] - ep_mean) / ep_std;
		}
		z[i] = dm_0[i] * r + sum(dm[i] .* x_left);
		y[i] = normal_rng(z[i], 0.01);
		m[i][1] = dm_0[i];
		for (j in 1:K) {
			m[i][j+1] = m[i][j] + dm[i][j];
		}
		theta[i] = atan(m[i]);
		dl[i] = dx ./ cos(theta[i][1:K]);
	}
}
