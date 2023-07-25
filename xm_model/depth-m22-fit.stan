/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
functions {
	real r_sig(real x) {
		return (1./pi()) * atan(x) + 0.5;
	}
	vector v_sig(vector x) {
		return (1./pi()) * atan(x) + 0.5;
	}
	matrix m_sig(matrix x) {
		return (1./pi()) * atan(x) + 0.5;
	}
	real r_sig_inv(real x) {
		return tan( pi() * (x - 0.5) );
	}
	vector v_sig_inv(vector x) {
		return tan( pi() * (x - 0.5) );
	}
}
data {
	/* Constants */
	int<lower=1> N;
	int<lower=1> M;
	real<lower=0> r;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	real rho_l;
	real rho_u;
	int<lower=1, upper=M> group[N];
	real y[N];
	/* Priors */
	real min_angle;
	real max_angle;
	real mu_yp_prior_mean;
	real<lower=0.> mu_yp_prior_std;
	real<lower=0.> sig_yp_prior_mean;
	real<lower=0.> sig_yp_prior_std;
	real frac_start_lower;
	real<lower=0.> frac_start_upper;
	real<lower=0.> lba_width;
	real<lower=0.> uba_width;

}
transformed data {
	real mu_dx = r / (K+1);
	matrix[K+1, K+1] Cp = rep_matrix(0., K+1, K+1);
	for (i in 1:(K+1)) {
		Cp[1:i,i] = rep_vector(1., i);
	}
	matrix[K+1, K+1] C = Cp';
	vector[K+1] dx = rep_vector(mu_dx, K+1);
	real xi_l = r_sig_inv(rho_l);
	real xi_u = r_sig_inv(rho_u);
	vector[K+1] rev_x = Cp * dx;
	real x_l1 = sum(rev_x);
	real x_l2 = sqrt(sum(square(rev_x)));
}
parameters {
	/* Param */
	vector<lower=min_angle, upper=min_angle+lba_width>[M] lba;
	vector<lower=max_angle-uba_width, upper=max_angle>[M] uba;
	vector<lower=frac_start_lower, upper=frac_start_upper>[M] rho_start;
	vector[M] mu_yp;
	vector<lower=0.>[M] sig_yp;
	matrix[K+1, N] dxi;
}
transformed parameters {
	/* secant transform */
	vector[M] theta_l = rho_l * uba + (1 - rho_l) * lba;
	vector[M] theta_u = rho_u * uba + (1 - rho_u) * lba;
	vector[M] beta_1 = (tan(theta_u) - tan(theta_l)) / (xi_u - xi_l);
	vector[M] beta_0 = tan(theta_u) - beta_1 * xi_u;
	/* dxi transform */
	vector[M] xi_0 = v_sig_inv(rho_start);
	vector[M] mu_dxi_shift_1 = r * beta_0;
	vector[M] mu_dxi_shift_2 = r * beta_1 .* xi_0;
	vector[M] mu_dxi_rescale = beta_1 * x_l1;
	vector[M] mu_dxi = (mu_yp - mu_dxi_shift_1 - mu_dxi_shift_2) ./  mu_dxi_rescale;
	vector[M] sig_dxi_rescale = beta_1 * x_l2;
	vector<lower=0>[M] sig_dxi = sig_yp ./ sig_dxi_rescale;
	/* latent variables */
	matrix[K+1, N] xi;
	matrix[K+1, N] theta;
	for (i in 1:N) {
		int j = group[i];
		xi[,i] = C * dxi[,i] + xi_0[j];
		theta[,i] = (uba[j] - lba[j]) * v_sig(xi[,i]) + lba[j];
	}
	matrix[K+1,N] m = tan(theta);
	vector[N] z = m' * dx;
}
model {
	lba ~ uniform(min_angle, min_angle + lba_width);
	uba ~ uniform(max_angle - uba_width, max_angle);
	rho_start ~ uniform(frac_start_lower, frac_start_upper);
	mu_yp ~ normal(mu_yp_prior_mean, mu_yp_prior_std);
	sig_yp ~ normal(sig_yp_prior_mean, sig_yp_prior_std);
	for (i in 1:N) {
		int j = group[i];
		dxi[,i] ~ normal(mu_dxi[j], sig_dxi[j]);
	}
	y ~ normal(z, 0.01);
}
