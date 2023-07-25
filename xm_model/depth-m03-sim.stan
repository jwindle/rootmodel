/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
data {
	int<lower=1> N;
	real<lower=0> r;
	real <lower=0> mu_dx;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	real xi_dm;
	real<lower=0> omega_dm;
	real shape_dm;
}
transformed data {
	vector[K] rvec = rep_vector(r, K);
	vector[K] dx = rep_vector(mu_dx, K);
	vector[K] x = cumulative_sum(dx);
	vector[K] x_left = fdim(rvec, x);
}
generated quantities {
	real dm_0[N];
	vector[K] dm[N];
	real z[N];
	real y[N];
	vector[K+1] m[N];
	vector[K+1] theta[N];
	vector[K] dl[N];
	for (i in 1:N) {
		dm_0[i] = skew_normal_rng(xi_dm, omega_dm, shape_dm);
		for (j in 1:K) {
			dm[i][j] = skew_normal_rng(xi_dm, omega_dm, shape_dm);
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
