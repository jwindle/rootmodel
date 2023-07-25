/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
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
	int<lower=1> M;
	real min_angle;
	real max_angle;
	real <lower=0> mu_dx;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	int<lower=1, upper=M> group[N]; /* group membership */
	real y[N];
	/* prior */
	real mu_yp_mean;
	real<lower=0.> mu_yp_std;
	real<lower=0.> sig_yp_mean;
	real<lower=0.> sig_yp_std;
	real xi_start_mean;
	real<lower=0.> xi_start_std;
}
transformed data {
	/* real<lower=0.> alpha = logistic_scale; */
	real<lower=0.> sig_yp_alpha = sig_yp_mean^2 / sig_yp_std^2 + 2;
	real<lower=0.> sig_yp_beta = sig_yp_mean * (sig_yp_alpha - 1);
	matrix[K+1, K+1] C = rep_matrix(0., K+1, K+1);
	for (i in 1:(K+1)) {
		C[1:i,i] = rep_vector(1., i);
	}
	vector[K+1] dx = rep_vector(mu_dx, K+1);
	dx[K+1] = r - mu_dx * K;
	matrix[N,M] IM = rep_matrix(0., N, M);
	for (i in 1:N) {
		int j = group[i];
		IM[i,j] = 1;
	}
	row_vector[K+1] oner_kp1 = rep_row_vector(1, K+1);
	/* secant approximation */
	real xi1 = xi_start_mean;
	real xi2 = -3.;
	real y1 = tan_real_logistic(xi1, min_angle, max_angle);
	real y2 = tan_real_logistic(xi2, min_angle, max_angle);
	real beta_1 = (y1 - y2) / (xi1 - xi2);
	real beta_0 = y2 - beta_1 * xi2;
	real mu_dxi_rescale = 2. / ((K+2) * (r*beta_1));
	real mu_dxi_shift_1 = r * beta_0 * mu_dxi_rescale;
	real sig_dxi_rescale = 6 * (K+1) / ((K+2) * (2*K+3) * (r*beta_1));
	vector[M] lba = rep_vector(min_angle, M);
	vector[M] uba = rep_vector(max_angle, M);
}
parameters {
	/* vector<lower=min_angle, upper=min_angle+0.1>[M] lba; */
	/* vector<lower=max_angle-0.1, upper=max_angle>[M] uba; */
	vector[M] mu_yp;
	vector<lower=0.>[M] sig_yp;
	vector[M] xi_start;
	matrix[N, K+1] dxi;
}
transformed parameters {
	vector[M] mu_dxi_shift_2 = 2 * xi_start / (K+2);
	vector[M] mu_dxi = mu_dxi_rescale * mu_yp - mu_dxi_shift_1 - mu_dxi_shift_2;
	vector<lower=0.>[M] sig_dxi = sig_dxi_rescale * sig_yp;
	matrix[N, K+1] xi;
	{
		matrix[N, K+1] xi_start_mat = (IM * xi_start) * oner_kp1;
		xi = dxi * C + xi_start_mat;
	}
	matrix[N, K+1] theta;
	{
		matrix[N, K+1] lb_mat = (IM * lba) * oner_kp1;
		matrix[N, K+1] ub_mat = (IM * uba) * oner_kp1;
		theta = (ub_mat - lb_mat) ./ (1. + exp(-1. * xi)) + lb_mat;
	}
	matrix[N, K+1] m = tan(theta);
	vector[N] z = m * dx;
}
model {
	/* lba ~ uniform(min_angle, min_angle+0.1); */
	/* uba ~ uniform(max_angle-0.1, max_angle); */
	xi_start ~ normal(xi_start_mean, xi_start_std);
	mu_yp ~ normal(mu_yp_mean, mu_yp_std);
	sig_yp ~ inv_gamma(sig_yp_alpha, sig_yp_beta);
	{
		/* matrix[N, K+1] mu_dxi_mat = (IM * mu_dxi) * oner_kp1; */
		/* matrix[N, K+1] sig_dxi_mat = (IM * sig_dxi) * oner_kp1; */
		/* dxi ~ normal(mu_dxi_mat, sig_dxi_mat);*/
		for (i in 1:N) {
			int g_i = group[i];
			dxi[i,] ~ normal(mu_dxi[g_i], sig_dxi[g_i]);
		}
	}
	y ~ normal(z, 0.01);
}
