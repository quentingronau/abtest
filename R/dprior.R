
#' Function for evaluating the prior density.
#' @title Prior Density
#' @param x1 numeric vector with values at which the prior density should be
#'   evaluated.
#' @param x2 if \code{what = "p1p2"}, value of p2 (i.e., the latent "success"
#'   probability in the experimental condition) at which the joint prior density
#'   should be evaluated. If \code{what = "p2givenp1"}, the given value of p1
#'   (i.e., the latent "success" probability in the control condition).
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
#' @param what character specifying for which quantity the prior density should
#'   be evaluated. Either \code{"logor"} (i.e., log odds ratio) , \code{"or"}
#'   (i.e., odds ratio), \code{"p1p2"} (i.e., the joint density of the latent
#'   "success" probability in the experimental and control condition),
#'   \code{"p1"} (i.e., latent "success" probability in the control condition),
#'   \code{"p2"} (i.e., latent "success" probability in the experimental
#'   condition), \code{"p2givenp1"} (i.e., conditional distribution of the
#'   latent "success" probability in the experimental condition given a
#'   "success" probability of \code{p1} in the control condition),
#'   \code{"rrisk"} (i.e., relative risk, the ratio of the "success" probability
#'   in the experimental and the control condition), or \code{"arisk"} (i.e.,
#'   absolute risk, the difference of the "success" probability in the
#'   experimental and control condition).
#' @param hypothesis character specifying whether to evaluate the two-sided
#'   prior density (i.e., "H1"), the one-sided prior density with lower
#'   truncation point (i.e., "H+"), or the one-sided prior density with upper
#'   truncation point (i.e., "H-").
#'
#' @return numeric vector with the values of the prior density.
#'
#' @note Internally, the test-relevant prior is always a normal prior on the log
#'   odds ratio, consequently, if \code{what} is not \code{"logor"}, the
#'   implied prior density for the quantity is returned.
#'
#' @author Quentin F. Gronau
#' @example examples/example.dprior.R
#' @export
dprior <- function(x1,
                   x2 = NULL,
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
  if ( ! (what %in%
          c("logor", "or", "p1p2", "p1", "p2",
            "p2givenp1", "rrisk", "arisk"))) {
    stop('what needs to be either "logor", "or", "p1p2", "p1", "p2",
         "p2givenp1", "rrisk", or "arisk"', call. = FALSE)
  }

  # check hypothesis
  if ( ! (hypothesis %in% c("H1", "H+", "H-"))) {
    stop('hypothesis needs to be either "H1", "H+", or "H-"',
         call. = FALSE)
  }

  if (what %in% c("logor", "or")) {

    out <- do.call(what = paste0("d", what),
                   args = list(x = x1,
                               mu_psi = prior_par[["mu_psi"]],
                               sigma_psi = prior_par[["sigma_psi"]],
                               hypothesis = hypothesis))

  } else if (what == "p1p2") {

    out <- do.call(what = paste0("d", what),
                   args = list(p1 = x1,
                               p2 = x2,
                               mu_psi = prior_par[["mu_psi"]],
                               sigma_psi = prior_par[["sigma_psi"]],
                               mu_beta = prior_par[["mu_beta"]],
                               sigma_beta = prior_par[["sigma_beta"]],
                               hypothesis = hypothesis))

  } else if (what == "p2givenp1") {

    out <- do.call(what = paste0("d", what),
                   args = list(p2 = x1,
                               p1 = x2,
                               mu_psi = prior_par[["mu_psi"]],
                               sigma_psi = prior_par[["sigma_psi"]],
                               mu_beta = prior_par[["mu_beta"]],
                               sigma_beta = prior_par[["sigma_beta"]],
                               hypothesis = hypothesis))

  } else if (what %in% c("rrisk", "arisk", "p1", "p2")) {

    out <- do.call(what = paste0("d", what),
                   args = list(x = x1,
                               mu_psi = prior_par[["mu_psi"]],
                               sigma_psi = prior_par[["sigma_psi"]],
                               mu_beta = prior_par[["mu_beta"]],
                               sigma_beta = prior_par[["sigma_beta"]],
                               hypothesis = hypothesis))

  }

  return(out)

}
