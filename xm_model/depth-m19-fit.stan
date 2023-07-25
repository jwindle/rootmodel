/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
data {
	int<lower=1> N;
	real<lower=0> r;
	int<lower=1> M;
	real min_angle;
	real max_angle;
	real<lower=0.> alpha;
	real <lower=0> mu_dx;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	int<lower=1, upper=M> group[N]; /* group membership */
	real y[N];
	/* prior */
	real mu_dxi_mean;
	real mu_dxi_std;
	real sig_dxi_mean;
	real sig_dxi_std;
	real xi_start_mean;
	real xi_start_std;
}
transformed data {
	/* real<lower=0.> alpha = logistic_scale; */
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
}
parameters {
	vector<lower=min_angle, upper=min_angle+0.1>[M] lb;
	vector<lower=max_angle-0.1, upper=max_angle>[M] ub;
	vector[M] xi_start;
	vector[M] mu_dxi;
	vector<lower=0.>[M] sig_dxi;
	matrix[N, K+1] dxi;
}
transformed parameters {
	matrix[N, K+1] xi;
	{
		matrix[N, K+1] xi_start_mat = (IM * xi_start) * oner_kp1;
		xi = dxi * C + xi_start_mat;
	}
	matrix[N, K+1] theta;
	{
		matrix[N, K+1] lb_mat = (IM * lb) * oner_kp1;
		matrix[N, K+1] ub_mat = (IM * ub) * oner_kp1;
		theta = (ub_mat - lb_mat) ./ (1. + exp(-1. * alpha * xi)) + lb_mat;
	}
	matrix[N, K+1] m = tan(theta);
	vector[N] z = m * dx;
}
model {
	lb ~ uniform(min_angle, min_angle+0.1);
	ub ~ uniform(max_angle-0.1, max_angle);
	xi_start ~ normal(xi_start_mean, xi_start_std);
	mu_dxi ~ normal(mu_dxi_mean, mu_dxi_std);
	sig_dxi ~ chi_square(1);
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
