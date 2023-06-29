/* -*- c-basic-offset: 4; -*- */
functions {
	real normal_lub_rng(real mu, real sigma, real lb, real ub) {
		/* From the Stan user guide */
		real p_lb = normal_cdf(lb, mu, sigma);
		real p_ub = normal_cdf(ub, mu, sigma);
		real u = uniform_rng(p_lb, p_ub);
		real y = mu + sigma * inv_Phi(u);
		return y;
	}
	vector x_left_from_mu_dx(real r, real mu_dx, int K) {
		vector[K+1] rvec = rep_vector(r, K+1);
		vector[K+1] dx = rep_vector(mu_dx, K+1);
		dx[1] = 0.0;
		vector[K+1] x = cumulative_sum(dx);
		vector[K+1] x_left = fdim(rvec, x);
		return x_left;
	}
	matrix Q_from_x_left(vector x_left) {
		int n = num_elements(x_left);
		vector[n] one = rep_vector(1., n);
		matrix[n,n] M = diag_matrix(one);
		real x_left_norm = sqrt(sum(square(x_left)));
		vector[n] x_left_til = x_left / x_left_norm;
		M[,1] = x_left_til;
		matrix[n, n] Q = qr_Q(M);
		return Q;
	}
	real l2_norm(vector x) {
		return sqrt(sum(square(x)));
	}
	matrix path_transform(real r, real mu_dx, int K) {
		vector[K+1] dx = rep_vector(mu_dx, K+1);
		vector[K+1] xgrid = cumulative_sum(dx);
		xgrid[K+1] = r;
		dx[1] = 0.0;
		vector[K+1] x = cumulative_sum(dx);
		matrix[K+1,K+1] A;
		for (i in 1:(K+1)) {
			vector[K+1] grid_i_rep = rep_vector(xgrid[i], K+1);
			vector[K+1] x_left = fdim(grid_i_rep, x);
			A[i,] = x_left';
		}
		return A;
	}
}
