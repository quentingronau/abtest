# elicit prior
prior_par <- elicit_prior(q = c(0.1, 0.3, 0.5),
                          prob = c(.025, .5, .975),
                          what = "arisk")
print(prior_par)

# plot elicited prior (absolute risk)
plot_prior(prior_par = prior_par, what = "arisk")

# plot corresponding normal prior on log odds ratio
plot_prior(prior_par = prior_par, what = "logor")
