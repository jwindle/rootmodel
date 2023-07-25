/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
functions {
	/* From the Stan user guide */
	real normal_lub_rng(real mu, real sigma, real lb, real ub) {
		real p_lb = normal_cdf(lb, mu, sigma);
		real p_ub = normal_cdf(ub, mu, sigma);
		real u = uniform_rng(p_lb, p_ub);
		real y = mu + sigma * inv_Phi(u);
		return y;
	}
}
data {
	int<lower=1> N;
	real<lower=0> r;
	real min_depth;
	real max_depth;
	real <lower=0> mu_dx;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	real mu_dm;
	real<lower=0> sig_dm;
}
transformed data {
	real lb = min_depth;
	real ub = max_depth;
	vector[K+1] rvec = rep_vector(r, K+1);
	vector[K+1] dx = rep_vector(mu_dx, K+1);
	dx[1] = 0.0;
	vector[K+1] x = cumulative_sum(dx);
	vector[K+1] x_left = fdim(rvec, x);
	vector[K+1] one = rep_vector(1., K+1);
	matrix[K+1,K+1] M = diag_matrix(one);
	real x_left_norm = sqrt(sum(square(x_left)));
	vector[K+1] x_left_til = x_left / x_left_norm;
	M[,1] = x_left_til;
	real lb_til = lb / x_left_norm;
	real ub_til = ub / x_left_norm;
	matrix[K+1, K+1] Q = qr_Q(M);
	matrix[K+1, K+1] Qp = (Q');
	/* w = Q' dm */
	vector[K+1] QpOne = Qp * one;
}
generated quantities {
	vector[K+1] w[N];
	vector[K+1] dm[N];
	real z[N];
	real y[N];
	vector[K+1] m[N];
	vector[K+1] theta[N];
	vector[K+1] dl[N];
	for (i in 1:N) {
		w[i][1] = normal_lub_rng(mu_dm * QpOne[1], sig_dm, lb_til, ub_til);
		for (j in 2:(K+1)) {
			w[i][j] = normal_rng(mu_dm * QpOne[j], sig_dm);
		}
		dm[i] = Q * w[i];
		z[i] = dot_product(dm[i], x_left); /* Could simplify to x_til_norm * w[i][1] */
		y[i] = normal_rng(z[i], 0.01);
		m[i][1] = dm[i][1];
		for (j in 1:K) {
			m[i][j+1] = m[i][j] + dm[i][j+1];
		}
		theta[i] = atan(m[i]);
		dl[i] = dx ./ cos(theta[i]);
	}
}
