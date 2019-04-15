# prior parameters
prior_par <- list(mu_psi = 0, sigma_psi = 1,
                  mu_beta = 0, sigma_beta = 1)

# evaluate prior CDF
pprior(q = 0.1, prior_par = prior_par, what = "logor")
pprior(q = 1.1, prior_par = prior_par, what = "or")
pprior(q = 1.05, prior_par = prior_par, what = "rrisk")
pprior(q = 0.02, prior_par = prior_par, what = "arisk")

# also works for vectors
pprior(q = c(-0.1, 0, 0.1, 0.2), prior_par = prior_par, what = "logor")
