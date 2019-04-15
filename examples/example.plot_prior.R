# prior parameters
prior_par <- list(mu_psi = 0, sigma_psi = 1,
                  mu_beta = 0, sigma_beta = 1)

# plot prior
plot_prior(prior_par = prior_par, what = "logor")
plot_prior(prior_par = prior_par, what = "or")
plot_prior(prior_par = prior_par, what = "p1p2")
plot_prior(prior_par = prior_par, what = "p1")
plot_prior(prior_par = prior_par, what = "p2")
plot_prior(prior_par = prior_par, what = "p2givenp1", p1 = 0.3)
plot_prior(prior_par = prior_par, what = "rrisk")
plot_prior(prior_par = prior_par, what = "arisk")
