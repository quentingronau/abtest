# synthetic data
data <- list(y1 = 10, n1 = 28, y2 = 14, n2 = 26)

# Bayesian A/B test with default settings
ab <- ab_test(data = data)
print(ab)

# different prior parameter settings
prior_par <- list(mu_psi = 0.2, sigma_psi = 0.8,
                  mu_beta = 0, sigma_beta = 0.7)
ab2 <- ab_test(data = data, prior_par = prior_par)
print(ab2)

# different prior probabilities
prior_prob <- c(.1, .3, .2, .4)
names(prior_prob) <- c("H1", "H+", "H-", "H0")
ab3 <- ab_test(data = data, prior_prob = prior_prob)
print(ab3)

# also possible to obtain posterior samples
ab4 <- ab_test(data = data, posterior = TRUE)

# plot parameter posterior
plot_posterior(x = ab4, what = "logor")
