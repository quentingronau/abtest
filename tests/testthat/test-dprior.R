context('dprior')

test_that("prior density integrates to one", {

  #------------------------------------------------------------------------
  # univariate densities
  #------------------------------------------------------------------------

  whats <- c("logor", "or", "p1", "p2", "p2givenp1", "rrisk", "arisk")
  hypotheses <- c("H1", "H+", "H-")
  prior_par <- list(mu_psi = -.2, sigma_psi = 1.2,
                    mu_beta = 0.3, sigma_beta = 0.6)

  for (what in whats) {
    for (hypothesis in hypotheses) {

      if (what == "p2givenp1") {
        x2 <- 0.3
      } else {
        x2 <- NULL
      }

      if (what == "logor") {

        bounds <- switch(hypothesis,
                         "H1" =  c(-Inf, Inf),
                         "H+" =  c(0, Inf),
                         "H-" =  c(-Inf, 0))

      } else if (what == "or") {

        bounds <- switch(hypothesis,
                         "H1" =  c(0, Inf),
                         "H+" =  c(exp(0), Inf),
                         "H-" =  c(0, exp(0)))

      } else if (what == "rrisk") {

        bounds <- switch(hypothesis,
                         "H1" =  c(0, Inf),
                         "H+" =  c(1, Inf),
                         "H-" =  c(0, 1))

      } else if (what == "arisk") {

        bounds <- switch(hypothesis,
                         "H1" =  c(-1, 1),
                         "H+" =  c(0, 1),
                         "H-" =  c(-1, 0))

      } else if (what == "p2givenp1") {

        bounds <- switch(hypothesis,
                         "H1" =  c(0, 1),
                         "H+" =  c(x2, 1),
                         "H-" =  c(0, x2))
      }


      int <- integrate(dprior, lower = bounds[1], upper = bounds[2],
                       x2 = x2, prior_par = prior_par, what = what,
                       hypothesis = hypothesis)$value

      expect_equal(int, expected = 1, tolerance = 1e-5)

    }
  }

  #------------------------------------------------------------------------
  # bivariate density
  #------------------------------------------------------------------------

  dp1 <- function(p1, prior_par, hypothesis) {

    vapply(p1, FUN = function(y, prior_par, hypothesis) {

      bounds <- switch(hypothesis,
                       "H1" =  c(0, 1),
                       "H+" =  c(y, 1),
                       "H-" =  c(0, y))

      integrate(function(p2, y, prior_par, hypothesis) {
        dprior(x1 = y, x2 = p2, prior_par = prior_par, what = "p1p2",
               hypothesis = hypothesis)
      },
      lower = 0, upper = 1,
      y = y, prior_par = prior_par,
      hypothesis = hypothesis)$value

    }, FUN.VALUE = 0, prior_par = prior_par,
    hypothesis = hypothesis)

  }

  for (hypothesis in hypotheses) {

    int <- integrate(dp1, lower = 0, upper = 1,
                     prior_par = prior_par,
                     hypothesis = hypothesis)$value

    expect_equal(int, expected = 1, tolerance = 1e-5)

  }

})
