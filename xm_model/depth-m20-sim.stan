/* -*- c-basic-offset: 4; -*- */
functions {
	real real_logistic(real x, real lb, real ub) {
		real d = ub - lb;
		return d / (1 + exp(-x)) + lb;
	}
	real tan_real_logistic(real x, real lb, real ub) {
		return tan(real_logistic(x, lb, ub));
	}
}
data {
	int<lower=1> N;
	real<lower=0> r;
	real lba;
	real uba;
	real <lower=0> mu_dx;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	/* real theta_start; */
	real xi_start;
	real mu_yp;
	real sig_yp;
}
transformed data {
	real lb = lba;
	real ub = uba;
	/* real<lower=0.> alpha = logistic_scale; */
	/* real xi_start = - log((ub - lb) / (theta_start - lb) - 1); */
	real xi2 = -3.;
	real y1 = tan_real_logistic(xi_start, lb, ub);
	real y2 = tan_real_logistic(xi2, lb, ub);
	real beta_1 = (y1 - y2) / (xi_start - xi2);
	real beta_0 = y2 - beta_1 * xi2;
	real mu_dxi_rescale = 2. / ((K+2) * (r*beta_1));
	real mu_dxi_shift_1 = r * beta_0 * mu_dxi_rescale;
	real mu_dxi_shift_2 = 2 * xi_start / (K+2);
	/* real mu_hat = (mu_yp - r*beta_0) / (r * beta_1) - xi_start; */
	/* real mu_dxi = mu_hat * 2 / (K+2); */
	real mu_dxi = mu_dxi_rescale * mu_yp - mu_dxi_shift_1 - mu_dxi_shift_2;
	real sig_dxi_rescale = 6 * (K+1) / ((K+2) * (2*K+3) * (r*beta_1));
	real<lower=0> sig_dxi = sig_dxi_rescale * sig_yp;
	matrix[K+1, K+1] C = rep_matrix(0., K+1, K+1);
	for (i in 1:(K+1)) {
		C[1:i,i] = rep_vector(1., i);
	}
	vector[K+1] dx = rep_vector(mu_dx, K+1);
	dx[K+1] = r - mu_dx * K;
	matrix [K+1, K+1] P = rep_matrix(0., K+1, K+1);
	for (i in 1:(K+1)) {
		P[1:i,i] = dx[1:i];
	}
}
generated quantities {
	real beta_0_copy = beta_0;
	real beta_1_copy = beta_1;
	real mu_dxi_copy = mu_dxi;
	real sig_dxi_copy = sig_dxi;
	matrix[N, K+1] dxi;
	for (i in 1:N) {
		for (j in 1:(K+1)) {
			dxi[i,j] = normal_rng(mu_dxi, sig_dxi);
		}
	}
	matrix[N, K+1] xi = dxi * C + xi_start;
	matrix[N, K+1] theta = (ub - lb) ./ (1 + exp(-xi)) + lb;
	matrix[N, K+1] m = tan(theta);
	vector[N] y = m * dx;
	matrix[N, K+1] path = m * P;
}
