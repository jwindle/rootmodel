/* -*- c-basic-offset: 4; -*- */
#include "common.stan"
data {
	int<lower=1> N;
	real<lower=0> r;
	real min_depth;
	real max_depth;
	real <lower=0> mu_dx;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	real tilde_mu;
	real tilde_sig;
	real gamma;
}
transformed data {
	real lb = min_depth;
	real ub = max_depth;
	vector[K+1] x_left = x_left_from_mu_dx(r, mu_dx, K);
	/* w = b * Q' dm */
	matrix[K+1, K+1] Q = Q_from_x_left(x_left);
	matrix[K+1, K+1] Qp = (Q');
	real x_left_l1 = sum(x_left);
	real x_left_l2 = l2_norm(x_left);
	real x_left_ratio = x_left_l2 / x_left_l1;
	vector[K+1] one = rep_vector(1.0, K+1);
	vector[K+1] rhoQpOne = Qp * one * x_left_ratio;
	/* real mu_dm = mu_y / x_left_l1; */
	/* real sig_dm = sig_y / x_left_l2; */
	matrix[K+1,K+1] PT = path_transform(r, mu_dx, K);
}
generated quantities {
	vector[K+1] w[N];
	vector[K+1] dm[N];
	real z[N];
	real y[N];
	vector[K+1] path[N];
	for (i in 1:N) {
		z[i] = normal_rng(tilde_mu, sqrt(gamma) * tilde_sig);
		w[i][1] = normal_lub_rng(z[i] * rhoQpOne[1], sqrt(1 - gamma) * tilde_sig, lb, ub);
		for (j in 2:(K+1)) {
			w[i][j] = normal_rng(z[i] * rhoQpOne[j], sqrt(1 - gamma) * tilde_sig);
		}
		dm[i] = Q * w[i] / x_left_l2;
		y[i] = dot_product(dm[i], x_left); /* Could simplify to x_til_norm * w[i][1] */
		path[i] = PT * dm[i];
	}
}

