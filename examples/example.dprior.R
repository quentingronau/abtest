# prior parameters
prior_par <- list(mu_psi = 0, sigma_psi = 1,
                  mu_beta = 0, sigma_beta = 1)

# prior density
dprior(x1 = 0.1, prior_par = prior_par, what = "logor")
dprior(x1 = 1.1, prior_par = prior_par, what = "or")
dprior(x1 = 0.49, x2 = 0.51, prior_par = prior_par, what = "p1p2")
dprior(x1 = 0.45, prior_par = prior_par, what = "p1")
dprior(x1 = 0.45, prior_par = prior_par, what = "p2")
dprior(x1 = 0.49, x2 = 0.51, prior_par = prior_par, what = "p2givenp1")
dprior(x1 = 1.05, prior_par = prior_par, what = "rrisk")
dprior(x1 = 0.02, prior_par = prior_par, what = "arisk")

# also works for vectors
dprior(x1 = c(-0.1, 0, 0.1, 0.2), prior_par = prior_par, what = "logor")
