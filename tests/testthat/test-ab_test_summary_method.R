context('ab_test summary method')

test_that("raw TRUE and FALSE yield similar results", {

  ab <- ab_test(data = list(y1 = 12000, n1 = 24000, y2 = 15000, n2 = 21000))
  s1 <- summary(ab, raw = FALSE)
  s2 <- summary(ab, raw = TRUE)

  expect_equal(s1$post$post_summary, s2$post$post_summary, tol = 0.001)

})
