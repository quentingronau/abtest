\dontrun{
# synthetic data
data <- list(y1 = 10, n1 = 28, y2 = 14, n2 = 26)

# Bayesian A/B test with default settings
ab <- ab_test(data = data)

# plot robustness check (i.e., prior sensitivity analysis)
p <- plot_robustness(ab)

# returned object contains the Bayes factors for the different prior settings
head(p)
}
