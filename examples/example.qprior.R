# prior parameters
prior_par <- list(mu_psi = 0, sigma_psi = 1,
                  mu_beta = 0, sigma_beta = 1)

# evaluate prior quantile function
qprior(p = .1, prior_par = prior_par, what = "logor")
qprior(p = .7, prior_par = prior_par, what = "or")
qprior(p = .9, prior_par = prior_par, what = "rrisk")
qprior(p = .7, prior_par = prior_par, what = "arisk")

# also works for vectors
qprior(p = c(.1, .2, .5, .7, .9), prior_par = prior_par, what = "logor")
