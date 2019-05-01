
#' Function for eliciting a prior distribution.
#' @title Elicit Prior
#' @param q vector with quantiles for the quantity of interest.
#' @param prob vector with probabilities corresponding to the quantiles (e.g., for
#'   the median the corresponding element of \code{prob} would need to be .5).
#' @param what character specifying for which quantity a prior should be
#'   elicited. Either \code{"logor"} (i.e., log odds ratio) , \code{"or"} (i.e.,
#'   odds ratio), \code{"rrisk"} (i.e., relative risk, the ratio of the
#'   "success" probability in the experimental and the control condition), or
#'   \code{"arisk"} (i.e., absolute risk, the difference of the "success"
#'   probability in the experimental and control condition).
#' @param hypothesis character specifying whether the provided quantiles
#'   correspond to a two-sided prior (i.e., "H1"), a one-sided prior with lower
#'   truncation point (i.e., "H+"), or a one-sided prior with upper truncation
#'   point (i.e., "H-").
#' @param mu_beta prior mean of the nuisance parameter \eqn{\beta} (i.e., the
#'   grand mean of the log odds). The default is 0.
#' @param sigma_beta prior standard deviation of the nuisance parameter
#'   \eqn{\beta} (i.e., the grand mean of the log odds). The default is 1.
#'
#' @details It is assumed that the prior on the grand mean of the log odds
#'   (i.e., \eqn{\beta}) is not the primary target of prior elicitation and is
#'   fixed (e.g., to a standard normal prior). The reason is that the grand mean
#'   nuisance parameter \eqn{\beta} is not the primary target of inference and
#'   changes in the prior on this nuisance parameter do not affect the results
#'   much in most cases (see Kass & Vaidyanathan, 1992). Nevertheless, it should
#'   be emphasized that the implemented approach allows users to set the prior
#'   parameters \code{mu_beta} and \code{sigma_beta} flexibly; the only
#'   constraint is that this takes place before the prior on the test-relevant
#'   log odds ratio parameter \eqn{\psi} is elicited. The \code{elicit_prior}
#'   function allows the user to elicit a prior not only in terms of the log
#'   odds ratio parameter \eqn{\psi}, but also in terms of the odds ratio, the
#'   relative risk (i.e., the ratio of the "success" probability in the
#'   experimental and the control condition), or the absolute risk (i.e., the
#'   difference of the "success" probability in the experimental and control
#'   condition). In case the prior is not elicited for the log odds ratio
#'   directly, the elicited prior is always translated to the closest
#'   corresponding normal prior on the log odds ratio. The prior parameters
#'   \code{mu_psi} and \code{sigma_psi} are obtained using least squares
#'   minimization.
#'
#' @return list with the elicited prior parameters. Specifically, this list
#'   consists of: \itemize{ \item\code{mu_psi} (prior mean for the normal prior
#'   on the test-relevant log odds ratio). \item \code{sigma_psi} (prior
#'   standard deviation for the normal prior on the test-relevant log odds
#'   ratio), \item \code{mu_beta} (prior mean for the normal prior on the grand
#'   mean of the log odds), \item \code{sigma_beta} (prior standard deviation
#'   for the normal prior on the grand mean of the log odds).} Note that the
#'   prior on the grand mean of the log odds is not part of the elicitation and
#'   is assumed to be fixed by the user (using the arguments \code{mu_beta} and
#'   \code{sigma_beta}). Consequently, the returned values for \code{mu_beta}
#'   and \code{sigma_beta} simply correspond to the input values.
#' @author Quentin F. Gronau
#' @references Kass, R. E., & Vaidyanathan, S. K. (1992). Approximate Bayes
#'   factors and orthogonal parameters, with application to testing equality of
#'   two binomial proportions. \emph{Journal of the Royal Statistical Society,
#'   Series B, 54}, 129-144. \url{https://doi.org/10.1111/j.2517-6161.1992.tb01868.x}
#' @example examples/example.elicit_prior.R
#'
#' @seealso The \code{\link{plot_prior}} function allows the user to visualize
#'   the elicited prior distribution.
#'
#' @export
elicit_prior <- function(q,
                         prob,
                         what = "logor",
                         hypothesis = "H1",
                         mu_beta = 0,
                         sigma_beta = 1) {

  ### check arguments ###

  # check q
  if (what == "logor") {

    qbounds <- switch(hypothesis,
                      "H1" =  c(-Inf, Inf),
                      "H+" =  c(0, Inf),
                      "H-" =  c(-Inf, 0))

  } else if (what == "or") {

    qbounds <- switch(hypothesis,
                      "H1" =  c(0, Inf),
                      "H+" =  c(exp(0), Inf),
                      "H-" =  c(0, exp(0)))

  } else if (what == "rrisk") {

    qbounds <- switch(hypothesis,
                      "H1" =  c(0, Inf),
                      "H+" =  c(1, Inf),
                      "H-" =  c(0, 1))

  } else if (what == "arisk") {

    qbounds <- switch(hypothesis,
                      "H1" =  c(-1, 1),
                      "H+" =  c(0, 1),
                      "H-" =  c(-1, 0))

  }

  if (any(q < qbounds[1]) || any(q > qbounds[2])) {
    stop(paste0('All specified quantiles need to be within ',
                qbounds[1], ' and ', qbounds[2]),
         call. = FALSE)
  }

  # check prob
  if (any(prob <= 0) || any(prob >= 1)) {
    stop('All specified probabilities need to be within 0 and 1',
         call. = FALSE)
  }

  # check hypothesis
  if ( ! (hypothesis %in% c("H1", "H+", "H-"))) {
    stop('hypothesis needs to be either "H1", "H+", or "H-"',
         call. = FALSE)
  }

  # check what
  if ( ! (what %in%
          c("logor", "or", "rrisk", "arisk"))) {
    stop('what needs to be either "logor", "or", "rrisk", or "arisk"',
         call. = FALSE)
  }

  # check sigma_beta
  if (sigma_beta <= 0) {
    stop('sigma_beta needs to be > 0',
         call. = FALSE)
  }

  start_value <- switch(hypothesis,
                        "H1" =  0,
                        "H+" =  .5,
                        "H-" =  -.5)

  o <- nlminb(start = c(start_value, 1), objective = .objective_elicit,
             lower = c(-Inf, 0),
             upper = c(Inf, Inf),
             q = q, prob = prob,
             mu_beta = mu_beta, sigma_beta = sigma_beta,
             what = what, hypothesis = hypothesis)
  par <- o$par
  names(par) <- c("mu_psi", "sigma_psi")

  out <- list(mu_psi = par[["mu_psi"]], sigma_psi = par[["sigma_psi"]],
              mu_beta = mu_beta, sigma_beta = sigma_beta)

  return(out)

}

.objective_elicit <- function(par, q, prob, mu_beta, sigma_beta,
                              what, hypothesis) {

  mu_psi <- par[1]
  sigma_psi <- par[2]
  prior_par <- list(mu_psi = mu_psi, sigma_psi = sigma_psi,
                    mu_beta = mu_beta, sigma_beta = sigma_beta)

  objvalue <- tryCatch({

    sum((pprior(q = q, prior_par = prior_par, what = what,
                hypothesis = hypothesis) - prob)^2)

  }, error = function(e) 1e5)

  return(objvalue)

}
