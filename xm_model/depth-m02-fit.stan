/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
data {
	int<lower=1> N;
	int<lower=1> M;
	real<lower=0> r;
	real <lower=0> mu_dx;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	real <upper=0.0> min_depth;
	real <upper=0.0> max_depth;
	real mu_dm_prior_mean;
	real mu_dm_prior_std;
	real sig_dm_prior_mean;
	real sig_dm_prior_std;
	int<lower=1, upper=M> group[N]; /* group membership */
	real y[N];
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
	matrix[K+1,K+1] A = diag_matrix(one);
	real x_left_norm = sqrt(sum(square(x_left)));
	vector[K+1] x_left_til = x_left / x_left_norm;
	A[,1] = x_left_til;
	real lb_til = lb / x_left_norm;
	real ub_til = ub / x_left_norm;
	matrix[K+1, K+1] Q = qr_Q(A);
	matrix[K+1, K+1] Qp = (Q');
	/* w = Q' dm */
	vector[K+1] QpOne = Qp * one;
}
parameters {
	real mu_dm[M];
	real<lower=0> sig_dm[M];
	real<lower=lb_til, upper=ub_til> w_0[N];
	vector[K] w[N];
}
transformed parameters {
	vector[K+1] dm[N];
	real z[N];
	for(i in 1:N) {
		dm[i] = (Q[:,1] * w_0[i]) + (Q[:,2:(K+1)] * w[i]);
		z[i] = dot_product(dm[i], x_left); /* Could simplify to x_til_norm * w[i][1] */
	}
}
model {
	for (i in 1:M) {
		mu_dm[i] ~ normal(mu_dm_prior_mean, mu_dm_prior_std);
		sig_dm[i] ~ normal(sig_dm_prior_mean, sig_dm_prior_std) T[0,];
		/* sig_dm[i] ~ inv_gamma(sig_dm_alpha, sig_dm_beta); */
	}
	for (i in 1:N) {
		int g_i = group[i];
		w_0[i] ~ normal(mu_dm[g_i] * QpOne[1], sig_dm[g_i]) T[lb_til, ub_til];
		for (j in 2:(K+1)) {
			w[i][j-1] ~ normal(mu_dm[g_i] * QpOne[j], sig_dm[g_i]);
		}
		y[i] ~ normal(z[i], 0.01);
	}
}
generated quantities {
	vector[K+1] m[N];
	vector[K+1] theta[N];
	vector[K+1] dl[N];
	for (i in 1:N) {
		m[i][1] = dm[i][1];
		for (j in 1:K) {
			m[i][j+1] = m[i][j] + dm[i][j+1];
		}
		theta[i] = atan(m[i]);
		dl[i] = dx ./ cos(theta[i]);
	}
}
