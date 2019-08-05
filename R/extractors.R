#' Extraction functions for ab objects
#'
#' Extraction functions for objects returned from the \code{\link{ab_test}}
#' function.
#'
#' @param x object of class \code{"ab"} as returned from \code{\link{ab_test}}.
#'
#' @return \code{get_bf} returns the Bayes factors in favor of "H1", "H+", and
#'   "H-" (compared to H0). \code{get_prior_prob} returns the prior
#'   probabilities of the hypotheses. \code{get_post_prob} returns the posterior
#'   probabilities of the hypotheses. \code{get_post_samples} returns posterior
#'   samples for the specified hypothesis.
#'
#' @example examples/example.extractors.R
#'
#' @name extractors
NULL


# extract bayes factors

#' @rdname extractors
#' @param log determines whether the log Bayes factors are returned.
#' @export
get_bf <- function(x, log = FALSE) {

  # make sure that object is of class ab
  if ( ! inherits(x, "ab")) {
    stop("x needs to be of class 'ab'", call. = FALSE)
  }

  if ( log ) {
    out <- x$logbf
  } else {
    out <- x$bf
  }

  return(out)

}

# extract prior probabilities

#' @rdname extractors
#' @export
get_prior_prob <- function(x) {

  # make sure that object is of class ab
  if ( ! inherits(x, "ab")) {
    stop("x needs to be of class 'ab'", call. = FALSE)
  }

  out <- x$input$prior_prob

  return(out)

}

# extract posterior probabilities

#' @rdname extractors
#' @export
get_post_prob <- function(x) {

  # make sure that object is of class ab
  if ( ! inherits(x, "ab")) {
    stop("x needs to be of class 'ab'", call. = FALSE)
  }

  out <- x$post_prob

  return(out)

}

# extract posterior samples

#' @rdname extractors
#' @param hypothesis determines for which hypothesis posterior samples are
#'   returned. Needs to be either "H1", "H+", or "H-" (the default is "H1").
#' @export
get_post_samples <- function(x, hypothesis = "H1") {

  # make sure that object is of class ab
  if ( ! inherits(x, "ab")) {
    stop("x needs to be of class 'ab'", call. = FALSE)
  }

  if (! hypothesis %in% c("H1", "H+", "H-")) {
    stop("hypothesis needs to either 'H1', 'H+', or 'H-'", call. = FALSE)
  }

  if ( ! x$input$posterior) {
    x <- ab_test(data = x$input$data,
                 prior_par = x$input$prior_par,
                 prior_prob = x$input$prior_prob,
                 nsamples = x$input$nsamples,
                 is_df = x$input$is_df,
                 posterior = TRUE)
  }

  name <- switch(hypothesis,
                 "H1" = "H1",
                 "H+" = "Hplus",
                 "H-" = "Hminus")
  out <- x$post[[name]]

  return(out)

}
