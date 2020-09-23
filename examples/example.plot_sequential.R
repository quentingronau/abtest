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

### 2.

# synthetic sequential data (observations alternate between the groups)
# this time provided in the alternative format
data2 <- data.frame(outcome = c(1, 1, 0, 1, 0, 1, 0, 1, 0, 1,
                                0, 1, 0, 1, 1, 1, 1, 1, 1, 0),
                    group = rep(c(1, 2), 10))

# conduct Bayesian A/B test with default settings
ab2 <- ab_test(data = data2)
print(ab2)

\donttest{
# produce sequential plot of posterior probabilities of the hypotheses
# (using recommended width and height values for saving to file)
cairo_pdf(file.path(tempdir(), "test_plot2.pdf"),
          width = 530 / 72, height = 400 / 72)
plot_sequential(ab2)
dev.off()
}

\dontrun{
### 3.
data(seqdata)

# conduct Bayesian A/B test with default settings
ab3 <- ab_test(data = seqdata)
print(ab3)

# produce sequential plot of posterior probabilities of the hypotheses
# (using recommended width and height values for saving to file)
cairo_pdf(file.path(tempdir(), "test_plot3.pdf"),
          width = 530 / 72, height = 400 / 72)
plot_sequential(ab3, thin = 4)
dev.off()
}
