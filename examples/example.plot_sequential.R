### 1.

# synthetic sequential data (observations alternate between the groups)
# note that the cumulative number of successes and trials need to be provided
data <- list(y1 = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4),
             n1 = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10),
             y2 = c(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9),
             n2 = c(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10))

# conduct Bayesian A/B test with default settings
ab <- ab_test(data = data)
print(ab)

\donttest{
# produce sequential plot of posterior probabilities of the hypotheses
# (using recommended width and height values for saving to file)
cairo_pdf(file.path(tempdir(), "test_plot.pdf"),
          width = 530 / 72, height = 400 / 72)
plot_sequential(ab)
dev.off()
}
\dontrun{
### 2.
data(seqdata)

# conduct Bayesian A/B test with default settings
ab2 <- ab_test(data = seqdata)
print(ab2)

# produce sequential plot of posterior probabilities of the hypotheses
# (using recommended width and height values for saving to file)
cairo_pdf(file.path(tempdir(), "test_plot2.pdf"),
          width = 530 / 72, height = 400 / 72)
plot_sequential(ab2, thin = 4)
dev.off()
}
