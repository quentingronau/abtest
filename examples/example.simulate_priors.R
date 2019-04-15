# prior parameters
prior_par <- list(mu_psi = 0, sigma_psi = 1,
                  mu_beta = 0, sigma_beta = 1)

# obtain prior samples
samples <- simulate_priors(nsamples = 1000, prior_par = prior_par)

# plot, e.g., prior samples for absolute risk
hist(samples$arisk)
