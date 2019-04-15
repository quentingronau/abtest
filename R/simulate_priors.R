
#' Function for simulating from the parameter prior distributions.
#' @title Simulate from Parameter Priors
#' @param nsamples number of samples.
#' @param prior_par list with prior parameters. This list needs to contain the
#'   following elements: \code{mu_psi} (prior mean for the normal prior on the
#'   test-relevant log odds ratio), \code{sigma_psi} (prior standard deviation
#'   for the normal prior on the test-relevant log odds ratio), \code{mu_beta}
#'   (prior mean for the normal prior on the grand mean of the log odds),
#'   \code{sigma_beta} (prior standard deviation for the normal prior on the
#'   grand mean of the log odds). Each of the elements needs to be a real number
#'   (the standard deviations need to be positive). The default are standard
#'   normal priors for both the log odds ratio parameter and the grand mean of
#'   the log odds parameter.
#' @param hypothesis character specifying whether to sample from a two-sided
#'   prior (i.e., "H1"), a one-sided prior with lower truncation point (i.e.,
#'   "H+"), or a one-sided prior with upper truncation point (i.e., "H-").
#'
#' @return a data frame with prior samples for the following quantities (see
#'   \code{?ab_test} for a description of the underlying model): \itemize{ \item
#'   \code{beta}: prior samples for the grand mean of the log odds. \item
#'   \code{psi}: prior samples for the log odds ratio. \item \code{p1}: prior
#'   samples for the latent "success" probability in the control group. \item
#'   \code{p2}: prior samples for the latent "success" probability in the
#'   experimental group. \item \code{logor}: prior samples for the log odds
#'   ratio (identical to \code{psi}, only included for easier reference). \item
#'   \code{or}: prior samples for the odds ratio. \item \code{rrisk}: prior
#'   samples for the relative risk (i.e., the ratio of the "success" probability
#'   in the experimental and the control condition). \item \code{arisk}: prior
#'   samples for the absolute risk (i.e., the difference of the "success"
#'   probability in the experimental and control condition)}.
#' @author Quentin F. Gronau
#' @example examples/example.simulate_priors.R
#' @importFrom stats rnorm
#' @importFrom truncnorm rtruncnorm
#' @export
simulate_priors <- function(nsamples,
                            prior_par = list(mu_psi = 0, sigma_psi = 1,
                                             mu_beta = 0, sigma_beta = 1),
                            hypothesis = "H1") {


  ### check arguments ###

  # check nsamples
  if (length(nsamples) > 1 || !is.numeric(nsamples) || nsamples <= 0) {
    stop('nsamples needs to be a single number larger than 0',
         call. = FALSE)
  }

  # check prior_par
  if ( ! is.list(prior_par) ||
       ! all(c("mu_psi", "sigma_psi", "mu_beta", "sigma_beta") %in%
          names(prior_par))) {
    stop('prior_par needs to be a named list with elements
         "mu_psi", "sigma_psi", "mu_beta", and "sigma_beta',
         call. = FALSE)
  }

  if (prior_par$sigma_psi <= 0 || prior_par$sigma_beta <= 0) {
    stop('sigma_psi and sigma_beta need to be larger than 0',
         call. = FALSE)
  }

  # check hypothesis
  if ( ! (hypothesis %in% c("H1", "H+", "H-"))) {
    stop('hypothesis needs to be either "H1", "H+", or "H-"',
         call. = FALSE)
  }

  bounds <- switch(hypothesis,
                   "H1" =  c(-Inf, Inf),
                   "H+" =  c(0, Inf),
                   "H-" =  c(-Inf, 0))

  beta <- rnorm(nsamples, mean = prior_par[["mu_beta"]],
                sd = prior_par[["sigma_beta"]])
  psi <- rtruncnorm(nsamples, a = bounds[1], b = bounds[2],
                    mean = prior_par[["mu_psi"]],
                    sd = prior_par[["sigma_psi"]])
  p1 <- inv_logit(beta - psi / 2)
  p2 <- inv_logit(beta + psi / 2)

  or <- p2 / (1 - p2) / (p1 / (1 - p1))
  logor <- log(or)

  rrisk <- p2 / p1
  arisk <- p2 - p1

  out <- data.frame(beta = beta, psi = psi, p1 = p1, p2 = p2,
                    logor = logor, or = or, rrisk = rrisk,
                    arisk = arisk)

  return(out)

}
