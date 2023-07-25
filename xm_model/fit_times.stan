/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
data {
	int<lower=1> n_obs;
	int<lower=1> n_times;
	int<lower=1> max_trials;
	real time_grid[n_times];
	int<lower=1, upper=n_times> time_group[n_obs];
	int<lower=0, upper=max_trials> y[n_obs];
	int<lower=1, upper=max_trials> n_trials[n_obs];
	real<lower=0> per_low;
	real<lower=0> per_high;
}
transformed data {
	real delta = 1e-9; // Try to ensure positive definite
	real sqrt_delta = sqrt(delta);
	matrix[n_obs, n_times] Z1 = rep_matrix(0., n_obs, n_times);
	for (i in 1:n_obs) {
		Z1[i,time_group[i]] = 1;
	}
}
parameters {
	real<lower=0.> rho1;
	real<lower=0.> sigma1;
	real<lower=0.> sigma_plant;
	real<lower=sqrt_delta> sigma_nug;
	real<lower=per_low, upper=per_high> per1;
	real mu;
	vector[n_times] eta1;
	// vector[n_times] ep_time;
	vector[n_obs] ep_plant;
}
transformed parameters {
	vector[n_times] f1;
	{
		matrix[n_times, n_times] L1;
		// matrix[n_times, n_times] L2;
		matrix[n_times, n_times] K1 = gp_periodic_cov(time_grid, sigma1, rho1, per1);
		// matrix[n_times, n_times] K2 = gp_exp_quad_cov(x, sigma2, rho2);
		
		// diagonal elements - this may not be necessary since I
		// control the design, i.e. `time_grid`.
		for (n in 1:n_times) {
			// K1[n, n] = K1[n, n] + delta;
			K1[n, n] = K1[n, n] + square(sigma_nug);
		}
		
		L1 = cholesky_decompose(K1);
		f1 = L1 * eta1;
	}
	// vector[n_times] lods1 = f1 + mu;
	// vector[n_times] lods2 = lods1 + ep_time;
	vector[n_times] lods2 = f1 + mu;
	vector[n_obs] z = Z1 * lods2 + ep_plant;
}
model {
	rho1 ~ inv_gamma(5, 5);
	sigma1 ~ std_normal();
	per1 ~ uniform(per_low, per_high);
	eta1 ~ std_normal();
	sigma_plant ~ std_normal();
	sigma_nug ~ std_normal() T[sqrt_delta,];
	mu ~ normal(-3.0, 1.5);
	// ep_time ~ normal(0.0, sigma_nug);
	ep_plant ~ normal(0.0, sigma_plant);
	y ~ binomial_logit(n_trials, z);
}
generated quantities {
	// vector[n_times] p1 = inv_logit(lods1);
	vector[n_times] p2 = inv_logit(lods2);
}
