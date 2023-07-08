/* -*- c-basic-offset: 4; -*- */
data {
	int<lower=1> N;
	real<lower=0> r;
	real min_angle;
	real max_angle;
	real theta_start;
	real<lower=0.> alpha;
	real <lower=0> mu_dx;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	real mu_dxi;
	real<lower=0> sig_dxi;
}
transformed data {
	real lb = min_angle;
	real ub = max_angle;
	/* real<lower=0.> alpha = logistic_scale; */
	real xi_start = - log((ub - lb) / (theta_start - lb) - 1) / alpha;
	matrix[K+1, K+1] C = rep_matrix(0., K+1, K+1);
	for (i in 1:(K+1)) {
		C[1:i,i] = rep_vector(1., i);
	}
	vector[K+1] dx = rep_vector(mu_dx, K+1);
	dx[K+1] = r - mu_dx * K;
}
generated quantities {
	matrix[N, K+1] dxi;
	for (i in 1:N) {
		for (j in 1:(K+1)) {
			dxi[i,j] = normal_rng(mu_dxi, sig_dxi);
		}
	}
	matrix[N, K+1] xi = dxi * C + xi_start;
	matrix[N, K+1] theta = (ub - lb) ./ (1 + exp(-alpha * xi)) + lb;
	matrix[N, K+1] m = tan(theta);
	vector[N] y = m * dx;
}
