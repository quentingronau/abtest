context('plot_robustness')

test_that("Different cores settings yield same results", {

  testthat::skip_on_cran()

  # synthetic data
  data <- list(y1 = 10, n1 = 28, y2 = 14, n2 = 26)

  # Bayesian A/B test with default settings
  ab <- ab_test(data = data)

  # BF10
  p1 <- plot_robustness(ab, mu_steps = 2, sigma_steps = 2, cores = 1)
  p2 <- plot_robustness(ab, mu_steps = 2, sigma_steps = 2, cores = 2)

  expect_equal(p1$bf, p2$bf)

  # BF0+
  p3 <- plot_robustness(ab, bftype = "BF0+",
                        mu_steps = 2,
                        sigma_steps = 2,
                        cores = 1)
  p4 <- plot_robustness(ab, bftype = "BF0+",
                        mu_steps = 2,
                        sigma_steps = 2,
                        cores = 1)
  expect_equal(p3$bf, p4$bf, tol = 0.1)

})
