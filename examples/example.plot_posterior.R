# synthetic data
data <- list(y1 = 10, n1 = 28, y2 = 14, n2 = 26)

# Bayesian A/B test with default settings
ab <- ab_test(data = data, posterior = TRUE)

# plot parameter posterior
plot_posterior(x = ab, what = "logor")
plot_posterior(x = ab, what = "or")
plot_posterior(x = ab, what = "p1p2")
plot_posterior(x = ab, what = "rrisk")
plot_posterior(x = ab, what = "arisk")

\donttest{
# example of good width and height values for saving to file
cairo_pdf(file.path(tempdir(), "test_plot.pdf"),
          width = 530 / 72, height = 400 / 72)
plot_posterior(ab, what = "p1p2")
dev.off()
}
