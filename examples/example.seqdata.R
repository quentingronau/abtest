\dontrun{

data(seqdata)

# conduct Bayesian A/B test with default settings
ab <- ab_test(data = seqdata)
print(ab)

# produce sequential plot of posterior probabilities of the hypotheses
plot_sequential(ab, thin = 4)

# example of good width and height values for saving to file
cairo_pdf(file.path(tempdir(), "test_plot.pdf"),
          width = 530 / 72, height = 400 / 72)
plot_sequential(ab)
dev.off()

}
