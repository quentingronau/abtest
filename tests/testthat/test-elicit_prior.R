context('elicit_prior')

test_that("elicit_prior works", {

  set.seed(1)

  prior_par <- list(mu_psi = -0.4, sigma_psi = 1.2,
                    mu_beta = 0.2, sigma_beta = 0.6)
  probs <- c(.25, .5, .75)

  whats <- c("logor", "or", "rrisk", "arisk")
  hypotheses <- c("H1", "H+", "H-")

  for (what in whats) {
    for (hypothesis in hypotheses) {

      sp <- simulate_priors(nsamples = 1e5, prior_par = prior_par,
                            hypothesis = hypothesis)

      q <- quantile(sp[[what]], probs = probs)
      ep <- elicit_prior(q = q, p = probs, what = what,
                         hypothesis = hypothesis,
                         mu_beta = prior_par[["mu_beta"]],
                         sigma_beta = prior_par[["sigma_beta"]])

      expect_equal(ep[["mu_psi"]], expected = prior_par[["mu_psi"]],
                   tolerance = 0.1)
      expect_equal(ep[["sigma_psi"]], expected = prior_par[["sigma_psi"]],
                   tolerance = 0.1)

    }
  }

})
