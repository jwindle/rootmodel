/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
functions {
	real normal_lub_rng(real mu, real sigma, real lb, real ub) {
		/* From the Stan user guide */
		real p_lb = normal_cdf(lb, mu, sigma);
		real p_ub = normal_cdf(ub, mu, sigma);
		if (p_lb == 1) {
			return mu;
		}
		if (p_ub == 0) {
			return mu;
		}
		real u = uniform_rng(p_lb, p_ub);
		real y = mu + sigma * inv_Phi(u);
		return y;
	}
}
data {
	int<lower=1> N;
	int<lower=1> K;
	real<lower=0> r;
	real mu_dm;
	real<lower=0> sig_dm;
	real mu_dx;
	real<lower=0> sig_dx;
}
transformed data {
	vector[K] rvec = rep_vector(r, K);
}
generated quantities {
	real dm_0[N];
	vector[K] dm[N];
	vector<lower=0.0, upper=r/(1.5)>[K] dx[N];
	vector[K] x[N];
	vector[K] x_left[N];
	real z[N];
	real y[N];
	vector[K+1] m[N];
	vector[K+1] theta[N];
	vector[K+1] delta_l[N];
	for (i in 1:N) {
		dm_0[i] = normal_rng(mu_dm, sig_dm);
		for (j in 1:K) {
			dm[i][j] = normal_rng(mu_dm, sig_dm);
			dx[i][j] = normal_lub_rng(mu_dx, sig_dx, 0.0, r/(1.5));
		}
		x[i] = cumulative_sum(dx[i]);
		x_left[i] = fdim(rvec, x[i]);
		z[i] = dm_0[i] * r + sum(dm[i] .* x_left[i]);
		y[i] = normal_rng(z[i], 0.01);
		m[i][1] = dm_0[i];
		for (j in 1:K) {
			m[i][j+1] = m[i][j] + dm[i][j];
		}
		theta[i] = atan(m[i]);
		delta_l[i] = dx[i] ./ cos(theta[i][1:K]);
	}
}
