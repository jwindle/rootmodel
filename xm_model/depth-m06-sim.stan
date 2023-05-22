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
	real<lower=0.> r;
	real<upper=0.> min_depth;
	real<upper=0.> max_depth;
	real <lower=0> mu_dx;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	real mu_dm;
	real<lower=0> sig_dm;
}
transformed data {
	vector[K+1] rvec = rep_vector(r, K+1);
	vector[K+1] dx = rep_vector(mu_dx, K+1);
	dx[1] = 0.0;
	vector[K+1] x = cumulative_sum(dx);
	vector[K+1] x_left = fdim(rvec, x);
	real ub = max_depth / sum(x_left);
	real lb = min_depth / sum(x_left);
}
generated quantities {
	vector[K+1] dm[N];
	real z[N];
	real y[N];
	vector[K+1] m[N];
	vector[K+1] theta[N];
	vector[K+1] dl[N];
	for (i in 1:N) {
		for (j in 1:(K+1)) {
			dm[i][j] = normal_lub_rng(mu_dm, sig_dm, lb, ub);
		}
		z[i] = dot_product(dm[i], x_left);
		y[i] = normal_rng(z[i], 0.01);
		m[i][1] = dm[i][1];
		for (j in 1:K) {
			m[i][j+1] = m[i][j] + dm[i][j+1];
		}
		theta[i] = atan(m[i]);
		dl[i] = dx ./ cos(theta[i]);
	}
}
