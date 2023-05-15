/* -*- c-basic-offset: 4; -*- */
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
	real lb;
	real ub;
	real <lower=0> mu_dx;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	real mu_dm;
	real<lower=0> sig_dm;
}
transformed data {
	vector[K+1] rvec = rep_vector(r, K+1);
	vector[K+1] dx = rep_vector(mu_dx, K+1);
	dx[1+1] = 0.0;
	vector[K+1] x = cumulative_sum(dx);
	vector[K+1] x_left = fdim(rvec, x);
	real alpha = sum(x_left);
	real beta = sqrt(sum(square(x_left)));
}
generated quantities {
	real z[N];
	real y[N];
	for (i in 1:N) {
		/* z[i] = normal_rng(mu_dm * alpha, sig_dm * beta); */
		z[i] = normal_lub_rng(mu_dm * alpha, sig_dm * beta, lb, ub);
		y[i] = normal_rng(z[i], 0.01);
	}
}
