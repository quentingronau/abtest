
#' Function for evaluating the prior quantile function.
#' @title Prior Quantile Function
#' @param p numeric vector with probabilities.
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
#' @param what character specifying for which quantity the prior quantile
#'   function should be evaluated. Either \code{"logor"} (i.e., log odds ratio)
#'   , \code{"or"} (i.e., odds ratio), \code{"rrisk"} (i.e., relative risk, the
#'   ratio of the "success" probability in the experimental and the control
#'   condition), or \code{"arisk"} (i.e., absolute risk, the difference of the
#'   "success" probability in the experimental and control condition).
#' @param hypothesis character specifying whether to evaluate the quantile
#'   function for a two-sided prior (i.e., "H1"), a one-sided prior with lower
#'   truncation point (i.e., "H+"), or a one-sided prior with upper truncation
#'   point (i.e., "H-").
#'
#' @return numeric vector with the values of the prior quantile function.
#' @author Quentin F. Gronau
#' @example examples/example.qprior.R
#' @export
qprior <- function(p,
                   prior_par = list(mu_psi = 0, sigma_psi = 1,
                                    mu_beta = 0, sigma_beta = 1),
                   what = "logor",
                   hypothesis = "H1") {


  ### check arguments ###

  # check p
  if (any(p < 0) || any(p > 1)) {
    stop('p needs to be between 0 and one', call. = FALSE)
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

  # check what
  if ( ! (what %in% c("logor", "or", "rrisk", "arisk"))) {
    stop('what needs to be either "logor", "or", "rrisk", or "arisk"',
         call. = FALSE)
  }

  # check hypothesis
  if ( ! (hypothesis %in% c("H1", "H+", "H-"))) {
    stop('hypothesis needs to be either "H1", "H+", or "H-"',
         call. = FALSE)
  }

  if (what == "logor") {

    bounds <- switch(hypothesis,
                     "H1" =  c(-Inf, Inf),
                     "H+" =  c(0, Inf),
                     "H-" =  c(-Inf, 0))
    start_value <- switch(hypothesis,
                          "H1" =  0,
                          "H+" =  .5,
                          "H-" =  -.5)

  } else if (what == "or") {

    bounds <- switch(hypothesis,
                     "H1" =  c(0, Inf),
                     "H+" =  c(exp(0), Inf),
                     "H-" =  c(0, exp(0)))
    start_value <- switch(hypothesis,
                          "H1" =  1,
                          "H+" =  1.5,
                          "H-" =  .5)

  } else if (what == "rrisk") {

    bounds <- switch(hypothesis,
                     "H1" =  c(0, Inf),
                     "H+" =  c(1, Inf),
                     "H-" =  c(0, 1))
    start_value <- switch(hypothesis,
                          "H1" =  1,
                          "H+" =  1.5,
                          "H-" =  .5)

  } else if (what == "arisk") {

    bounds <- switch(hypothesis,
                     "H1" =  c(-1, 1),
                     "H+" =  c(0, 1),
                     "H-" =  c(-1, 0))
    start_value <- switch(hypothesis,
                          "H1" =  0,
                          "H+" =  .5,
                          "H-" =  -.5)

  }

  if (what == "logor") {

   out <- qtruncnorm(p, a = bounds[1], b = bounds[2],
                     mean = prior_par$mu_psi,
                     sd = prior_par$sigma_psi)

  } else {

    out <- vapply(p, FUN = function(x,
                                    start_value,
                                    bounds,
                                    prior_par,
                                    what,
                                    hypothesis) {

      nlminb(start = start_value, objective = function(q,
                                                       p,
                                                       prior_par,
                                                       what,
                                                       hypothesis) {

        (pprior(q, prior_par = prior_par, what = what,
                hypothesis = hypothesis) - p)^2

      }, lower = bounds[1], upper = bounds[2], p = x, prior_par = prior_par,
      what = what, hypothesis = hypothesis)$par

    }, FUN.VALUE = 0, start_value = start_value, bounds = bounds,
    prior_par = prior_par, what = what, hypothesis = hypothesis)

  }

  return(out)

}
