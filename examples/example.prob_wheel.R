# synthetic data
data <- list(y1 = 10, n1 = 28, y2 = 14, n2 = 26)

# Bayesian A/B test with default settings
ab <- ab_test(data = data)
print(ab)

# visualize prior probabilities of the hypotheses
prob_wheel(ab, type = "prior")

# visualize posterior probabilities of the hypotheses
prob_wheel(ab, type = "posterior")
