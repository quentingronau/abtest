
#' Function for evaluating the prior cumulative distribution function (CDF).
#' @title Prior Cumulative Distribution Function (CDF)
#' @param q numeric vector with quantiles.
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
#' @param what character specifying for which quantity the prior CDF should be
#'   evaluated. Either \code{"logor"} (i.e., log odds ratio) , \code{"or"}
#'   (i.e., odds ratio), \code{"rrisk"} (i.e., relative risk, the ratio of the
#'   "success" probability in the experimental and the control condition), or
#'   \code{"arisk"} (i.e., absolute risk, the difference of the "success"
#'   probability in the experimental and control condition).
#' @param hypothesis character specifying whether to evaluate the CDF for a
#'   two-sided prior (i.e., "H1"), a one-sided prior with lower truncation point
#'   (i.e., "H+"), or a one-sided prior with upper truncation point (i.e.,
#'   "H-").
#'
#' @return numeric vector with the values of the prior CDF.
#' @note Internally, the test-relevant prior is always a normal prior on the log
#'   odds ratio, consequently, if \code{what} is not \code{"logor"}, the
#'   implied prior CDF for the quantity is returned.
#' @author Quentin F. Gronau
#' @example examples/example.pprior.R
#' @export
pprior <- function(q,
                   prior_par = list(mu_psi = 0, sigma_psi = 1,
                                    mu_beta = 0, sigma_beta = 1),
                   what = "logor",
                   hypothesis = "H1") {


  ### check arguments ###

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

  if (what %in% c("logor", "or")) {

    out <- do.call(what = paste0("p", what),
                   args = list(q = q,
                               mu_psi = prior_par[["mu_psi"]],
                               sigma_psi = prior_par[["sigma_psi"]],
                               hypothesis = hypothesis))

  } else if (what %in% c("rrisk", "arisk")) {

    out <- do.call(what = paste0("p", what),
                   args = list(q = q,
                               mu_psi = prior_par[["mu_psi"]],
                               sigma_psi = prior_par[["sigma_psi"]],
                               mu_beta = prior_par[["mu_beta"]],
                               sigma_beta = prior_par[["sigma_beta"]],
                               hypothesis = hypothesis))

  }

  return(out)

}
