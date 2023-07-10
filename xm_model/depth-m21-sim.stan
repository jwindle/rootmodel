/* -*- c-basic-offset: 4; -*- */
functions {
	real real_logistic(real x, real lb, real ub) {
		real d = ub - lb;
		return d / (1 + exp(-x)) + lb;
	}
	real inv_real_logistic(real theta, real lb, real ub) {
		real d = ub - lb;
		return -log(d / (theta - lb) - 1.);
	}
	vector vec_logistic(vector v, vector lb, vector ub) {
		int n = num_elements(v);
		vector[n] d = ub - lb;
		return d ./ (1 + exp(-v)) + lb;
	}
	vector inv_vec_logistic(vector theta, vector lb, vector ub) {
		int n = num_elements(theta);
		vector[n] d = ub - lb;
		return -log(d ./ (theta - lb) - 1.);
	}
	real tan_real_logistic(real x, real lb, real ub) {
		return tan(real_logistic(x, lb, ub));
	}
	real inv_real_logistic_atan(real x, real lb, real ub) {
		return inv_real_logistic(atan(x), lb, ub);
	}
}
data {
	int<lower=1> N;
	real<lower=0> r;
	real min_angle;
	real max_angle;
	real frac_start_lower;
	real frac_start_upper;
	real <lower=0.> mu_dx;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	real mu_yp;
	real<lower=0.> sig_yp;
	real lba;
	real uba;
	real<lower=frac_start_lower, upper=frac_start_upper> start_frac;
}
transformed data {
	/* Consants */
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
	/* Prior transforms */
	/* secant approximation */
	real theta0 = (max_angle - min_angle) * 0.5*(frac_start_lower + frac_start_upper) + min_angle;
	real xi1 = inv_real_logistic(theta0, min_angle, max_angle);
	real xi2 = -3.;
	real y1 = tan_real_logistic(xi1, min_angle, max_angle);
	real y2 = tan_real_logistic(xi2, min_angle, max_angle);
	real beta_1 = (y1 - y2) / (xi1 - xi2);
	real beta_0 = y2 - beta_1 * xi2;
	real mu_dxi_rescale = 2. / ((K+2) * (r*beta_1));
	real mu_dxi_shift_1 = r * beta_0 * mu_dxi_rescale;
	real sig_dxi_rescale = 6 * (K+1) / ((K+2) * (2*K+3) * (r*beta_1));
	/* Transformed parameters */
	real theta_start = lba + start_frac * (uba - lba);
	real xi_start = inv_real_logistic(theta_start, lba, uba);
	real mu_dxi_shift_2 = 2 * xi_start / (K+2);
	real mu_dxi = mu_dxi_rescale * mu_yp - mu_dxi_shift_1 - mu_dxi_shift_2;
	real sig_dxi = sig_dxi_rescale * sig_yp;
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
	matrix[N, K+1] theta = (uba - lba) ./ (1 + exp(-xi)) + lba;
	matrix[N, K+1] m = tan(theta);
	vector[N] y = m * dx;
	matrix[N, K+1] path = m * P;
}
