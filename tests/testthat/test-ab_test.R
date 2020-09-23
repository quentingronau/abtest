context('ab_test')

test_that("basic ab_test works", {

  library(truncnorm)

  #------------------------------------------------------------------------
  # functions for numerical integration
  #------------------------------------------------------------------------

  H0int <- function(beta, data, prior_par) {

    p <- plogis(beta)
    p^(data$y1 + data$y2) *
      (1 - p)^(data$n1 + data$n2 - data$y1 - data$y2) *
      dnorm(beta, prior_par$mu_beta, prior_par$sigma_beta)

  }

  inner_integrand <- function(beta, psi, data, prior_par) {

    p1 <- plogis(beta - psi / 2)
    p2 <- plogis(beta + psi / 2)

    p1^data$y1 * (1 - p1)^(data$n1 - data$y1) *
      p2^data$y2 * (1 - p2)^(data$n2 - data$y2) *
      dnorm(beta, prior_par$mu_beta, prior_par$sigma_beta)

  }

  Haltint <- function(psi, data, prior_par, hypothesis) {

    bounds <- switch(hypothesis,
                     "H1" =  c(-Inf, Inf),
                     "H+" =  c(0, Inf),
                     "H-" =  c(-Inf, 0))

    vapply(psi, FUN = function(x) {

      integrate(inner_integrand, lower = -Inf, upper = Inf,
                psi = x, data = data, prior_par = prior_par)$value

    }, FUN.VALUE = 0) * dtruncnorm(psi, a = bounds[1], b = bounds[2],
                                   mean = prior_par$mu_psi,
                                   sd = prior_par$sigma_psi)

  }

  mean_psi_int <- function(psi, data, prior_par, hypothesis) {

    bounds <- switch(hypothesis,
                     "H1" =  c(-Inf, Inf),
                     "H+" =  c(0, Inf),
                     "H-" =  c(-Inf, 0))

    vapply(psi, FUN = function(x) {

      integrate(inner_integrand, lower = -Inf, upper = Inf,
                psi = x, data = data, prior_par = prior_par)$value

    }, FUN.VALUE = 0) * psi * dtruncnorm(psi, a = bounds[1], b = bounds[2],
                                         mean = prior_par$mu_psi,
                                         sd = prior_par$sigma_psi)

  }

  var_psi_int <- function(psi, m, data, prior_par, hypothesis) {

    bounds <- switch(hypothesis,
                     "H1" =  c(-Inf, Inf),
                     "H+" =  c(0, Inf),
                     "H-" =  c(-Inf, 0))

    vapply(psi, FUN = function(x) {

      integrate(inner_integrand, lower = -Inf, upper = Inf,
                psi = x, data = data, prior_par = prior_par)$value

    }, FUN.VALUE = 0) * (psi - m)^2 * dtruncnorm(psi, a = bounds[1],
                                                 b = bounds[2],
                                                 mean = prior_par$mu_psi,
                                                 sd = prior_par$sigma_psi)

  }

  mean_beta_int <- function(beta, data, prior_par, hypothesis) {

    bounds <- switch(hypothesis,
                     "H1" =  c(-Inf, Inf),
                     "H+" =  c(0, Inf),
                     "H-" =  c(-Inf, 0))

    vapply(beta, FUN = function(x) {

      integrate(function(psi, beta, data, prior_par, hypothesis, bounds) {

        p1 <- plogis(beta - psi / 2)
        p2 <- plogis(beta + psi / 2)

        p1^data$y1 * (1 - p1)^(data$n1 - data$y1) *
          p2^data$y2 * (1 - p2)^(data$n2 - data$y2) *
          dtruncnorm(psi, a = bounds[1], b = bounds[2],
                     mean = prior_par$mu_psi,
                     sd = prior_par$sigma_psi)

      }, lower = -Inf, upper = Inf,
      beta = x, data = data, prior_par = prior_par,
      hypothesis = hypothesis, bounds = bounds)$value

    }, FUN.VALUE = 0) * beta *
      dnorm(beta, prior_par$mu_beta, prior_par$sigma_beta)

  }

  var_beta_int <- function(beta, m, data, prior_par, hypothesis) {

    bounds <- switch(hypothesis,
                     "H1" =  c(-Inf, Inf),
                     "H+" =  c(0, Inf),
                     "H-" =  c(-Inf, 0))

    vapply(beta, FUN = function(x) {

      integrate(function(psi, beta, data, prior_par, hypothesis, bounds) {

        p1 <- plogis(beta - psi / 2)
        p2 <- plogis(beta + psi / 2)

        p1^data$y1 * (1 - p1)^(data$n1 - data$y1) *
          p2^data$y2 * (1 - p2)^(data$n2 - data$y2) *
          dtruncnorm(psi, a = bounds[1], b = bounds[2],
                     mean = prior_par$mu_psi,
                     sd = prior_par$sigma_psi)

      }, lower = -Inf, upper = Inf,
      beta = x, data = data, prior_par = prior_par,
      hypothesis = hypothesis, bounds = bounds)$value

    }, FUN.VALUE = 0) * (beta - m)^2 *
      dnorm(beta, prior_par$mu_beta, prior_par$sigma_beta)

  }

  #------------------------------------------------------------------------
  # example 1
  #------------------------------------------------------------------------

  set.seed(1)

  data <- list(y1 = 11, n1 = 15, y2 = 5, n2 = 13)
  prior_par <- list(mu_psi = -0.4, sigma_psi = 1.2,
                    mu_beta = 0.2, sigma_beta = 0.6)
  ab <- ab_test(data = data, prior_par = prior_par)

  # compare to numerical integration results
  mlH0 <- integrate(H0int, lower = -Inf, upper = Inf, data = data,
                    prior_par = prior_par)$value
  mlH1 <- integrate(Haltint, lower = -Inf, upper = Inf, data = data,
                    prior_par = prior_par, hypothesis = "H1")$value
  mlHplus <- integrate(Haltint, lower = 0, upper = Inf, data = data,
                       prior_par = prior_par, hypothesis = "H+")$value
  mlHminus <- integrate(Haltint, lower = -Inf, upper = 0, data = data,
                        prior_par = prior_par, hypothesis = "H-")$value

  BF10 <- mlH1 / mlH0
  BFplus0 <- mlHplus / mlH0
  BFminus0 <- mlHminus / mlH0

  # test that the Bayes factors match
  expect_equal(ab$bf$bf10, expected = BF10, tolerance = 0.1)
  expect_equal(ab$bf$bfplus0, expected = BFplus0, tolerance = 0.1)
  expect_equal(ab$bf$bfminus0, expected = BFminus0, tolerance = 0.1)

  # test that posterior probabilities match
  prior_prob <- ab$input$prior_prob
  mls <- c(mlH1, mlHplus, mlHminus, mlH0)
  post_prob <- mls * prior_prob / sum(mls * prior_prob)
  expect_equal(ab$post_prob, expected = post_prob, tolerance = 0.05)

  prior_prob2 <- c(0.1, 0.3, 0.2, 0.4)
  names(prior_prob2) <- names(prior_prob)
  ab2 <- ab_test(data = data, prior_par = prior_par,
                 prior_prob = prior_prob2)
  post_prob2 <- mls * prior_prob2 / sum(mls * prior_prob2)
  expect_equal(ab2$post_prob, expected = post_prob2, tolerance = 0.05)

  #------------------------------------------------------------------------
  # test that posterior samples look ok
  #------------------------------------------------------------------------

  ab3 <- ab_test(data = data, prior_par = prior_par, posterior = TRUE)
  ml <- list("H1" = mlH1, "H+" = mlHplus, "H-" = mlHminus)

  for (hypothesis in c("H1", "H+", "H-")) {

    bounds <- switch(hypothesis,
                     "H1" =  c(-Inf, Inf),
                     "H+" =  c(0, Inf),
                     "H-" =  c(-Inf, 0))

    hypothesis2 <- switch(hypothesis,
                          "H1" =  "H1",
                          "H+" =  "Hplus",
                          "H-" =  "Hminus")

    # mean psi
    mnum <- integrate(mean_psi_int, lower = bounds[1], upper = bounds[2],
                      data = data, prior_par = prior_par,
                      hypothesis = hypothesis)$value / ml[[hypothesis]]
    expect_equal(mean(ab3$post[[hypothesis2]]$psi), expected = mnum,
                 tolerance = 0.05)

    # variance psi
    vnum <- integrate(var_psi_int, lower = bounds[1], upper = bounds[2],
                      m = mnum, data = data, prior_par = prior_par,
                      hypothesis = hypothesis)$value / ml[[hypothesis]]
    expect_equal(var(ab3$post[[hypothesis2]]$psi), expected = vnum,
                 tolerance = 0.05)

    # mean beta
    mnum <- integrate(mean_beta_int, lower = -Inf, upper = Inf,
                      data = data, prior_par = prior_par,
                      hypothesis = hypothesis)$value / ml[[hypothesis]]
    expect_equal(mean(ab3$post[[hypothesis2]]$beta), expected = mnum,
                 tolerance = 0.05)

    # variance beta
    vnum <- integrate(var_beta_int, lower = -Inf, upper = Inf, m = mnum,
                      data = data, prior_par = prior_par,
                      hypothesis = hypothesis)$value / ml[[hypothesis]]
    expect_equal(var(ab3$post[[hypothesis2]]$beta), expected = vnum,
                 tolerance = 0.05)

  }

})

test_that("different sequential data formats yield same result", {

  set.seed(1)

  data <- list(y1 = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4),
               n1 = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10),
               y2 = c(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9),
               n2 = c(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10))
  data_seq <- data.frame(outcome = c(1, 1, 0, 1, 0, 1, 0, 1, 0, 1,
                                     0, 1, 0, 1, 1, 1, 1, 1, 1, 0),
                         group = rep(c(1, 2), 10))
  data_seq_mat <- as.matrix(data_seq)

  ab <- ab_test(data)
  ab_seq <- ab_test(data_seq)
  ab_seq_mat <- ab_test(data_seq_mat)

  # test that the Bayes factors match
  expect_equal(ab_seq$bf$bf10, expected = ab$bf$bf10)
  expect_equal(ab_seq_mat$bf$bf10, expected = ab$bf$bf10)

})

test_that("y and n alternative data specification works", {

  set.seed(1)

  data <- list(y1 = 11, n1 = 15, y2 = 5, n2 = 13)

  ab <- ab_test(data)
  ab2 <- ab_test(y = c(11, 5), n = c(15, 13))

  # test that the Bayes factors match
  expect_equal(ab2$bf$bf10, expected = ab$bf$bf10)

})
