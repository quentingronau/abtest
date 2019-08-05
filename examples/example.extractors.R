# synthetic data
data <- list(y1 = 10, n1 = 28, y2 = 14, n2 = 26)

# Bayesian A/B test with default settings
ab <- ab_test(data = data, posterior = TRUE)

# extract Bayes factors
get_bf(ab)

# extract prior probabilities
get_prior_prob(ab)

# extract posterior probabilities
get_post_prob(ab)

# extract posterior samples for H1
s <- get_post_samples(ab, hypothesis = "H1")
