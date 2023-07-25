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
	real<lower=0> r;
	/* real <lower=0> mu_dx; */
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	real rho_l;
	real rho_u;
	/* Param */
	real lba;
	real uba;
	real rho_start;
	real mu_yp;
	real sig_yp;
}
transformed data {
	real lb = lba;
	real ub = uba;
	real mu_dx = r / (K+1);
	/* Constants */
	matrix[K+1, K+1] Cp = rep_matrix(0., K+1, K+1);
	for (i in 1:(K+1)) {
		Cp[1:i,i] = rep_vector(1., i);
	}
	vector[K+1] dx = rep_vector(mu_dx, K+1);
	matrix [K+1, K+1] P = rep_matrix(0., K+1, K+1);
	for (i in 1:(K+1)) {
		P[1:i,i] = dx[1:i];
	}
	real xi_0 = r_sig_inv(rho_start);
	/* Secant transform */
	real xi_l = r_sig_inv(rho_l);
	real xi_u = r_sig_inv(rho_u);
	real theta_l = rho_l * ub + (1 - rho_l) * lb;
	real theta_u = rho_u * ub + (1 - rho_u) * lb;
	real beta_1 = (tan(theta_u) - tan(theta_l)) / (xi_u - xi_l);
	real beta_0 = tan(theta_u) - beta_1 * xi_u;
	/* dxi transform */
	vector[K+1] rev_x = Cp * dx;
	real x_l1 = sum(rev_x);
	real x_l2 = sqrt(sum(square(rev_x)));
	real mu_dxi_shift_1 = r * beta_0;
	real mu_dxi_shift_2 = r * beta_1 * xi_0;
	real mu_dxi_rescale = beta_1 * x_l1;
	real mu_dxi = (mu_yp - mu_dxi_shift_1 - mu_dxi_shift_2) / mu_dxi_rescale;
	real sig_dxi_rescale = beta_1 * x_l2;
	real<lower=0> sig_dxi = sig_yp / sig_dxi_rescale;
}
generated quantities {
	real beta_0_copy = beta_0;
	real beta_1_copy = beta_1;
	real xi_0_copy = xi_0;
	real mu_dxi_copy = mu_dxi;
	real sig_dxi_copy = sig_dxi;
	real lb_copy = lb;
	real ub_copy = ub;
	matrix[N, K+1] dxi;
	for (i in 1:N) {
		for (j in 1:(K+1)) {
			dxi[i,j] = normal_rng(mu_dxi, sig_dxi);
		}
	}
	matrix[N, K+1] xi = dxi * Cp + xi_0;
	matrix[N, K+1] theta = (ub - lb) * m_sig(xi) + lb;
	matrix[N, K+1] m = tan(theta);
	vector[N] y = m * dx;
	matrix[N, K+1] path = m * P;
}
