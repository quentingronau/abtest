#' Methods for ab objects
#'
#' Methods defined for objects returned from the \code{\link{ab_test}} function.
#'
#' @param object,x object of class \code{ab} as returned from
#'   \code{\link{ab_test}}.
#' @param digits number of digits to print for the summary.
#' @param raw if \code{TRUE}, the raw posterior samples are used to estimate the
#'   mean, sd, and quantiles for the summary of the posterior. If \code{FALSE},
#'   parametric fits to the marginal posteriors are used to obtain the mean, sd,
#'   and quantiles. Specifically, a normal distribution is fitted for \code{psi
#'   (logor)} and \code{beta}; a log-normal distribution is fitted for \code{or}
#'   and \code{rrisk}; beta distributions are fitted for \code{p1} and
#'   \code{p2}; a scaled beta distribution is fitted for \code{arisk}. These
#'   distributional fits are also used in \code{\link{plot_posterior}}.
#' @param ... further arguments, currently ignored.
#'
#' @return The \code{print} methods prints the Bayes factors, prior
#'   probabilities of the hypotheses, and posterior probabilities of the
#'   hypotheses (and returns nothing).
#'
#'   The \code{plot} method visualizes the prior probabilities of the hypotheses
#'   and posterior probabilities of the hypotheses (the next plots is obtained
#'   by hitting Return) using the \code{\link{prob_wheel}} function.
#'
#'   The \code{summary} methods returns the \code{ab} object that is guaranteed
#'   to contain posterior samples (i.e., it adds posterior samples if they were
#'   not included already). Additionally, it adds to the object a posterior
#'   summary matrix (i.e., \code{ab$post$post_summary}) for the posterior under
#'   H1 and the arguments \code{digits} (used for printing) and \code{raw}
#'   (added to \code{ab$input}).
#'
#' @name ab-methods
NULL


# summary method

#' @rdname ab-methods
#' @method summary ab
#' @export
summary.ab <- function(object, digits = 3, raw = FALSE, ...) {

  # make sure that object contains posterior samples
  if ( ! object$input$posterior) {
    object <- ab_test(data = object$input$data,
                 prior_par = object$input$prior_par,
                 prior_prob = object$input$prior_prob,
                 nsamples = object$input$nsamples,
                 is_df = object$input$is_df,
                 posterior = TRUE)
  }

  # obtain posterior summary under H1
  pars <- c("beta", "psi (logor)", "or", "arisk", "rrisk", "p1", "p2")
  post_summary <- matrix(nrow = length(pars), ncol = 7)
  rownames(post_summary) <- pars
  colnames(post_summary) <- c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%")
  p <- c(0.025, 0.25, 0.5, 0.75, 0.975)

  if ( ! raw) {
    fit <- vector("list", 7)
    names(fit) <- pars
    for (i in pars) {
      whati <- ifelse(i == "psi (logor)",  "logor", i)
      fit[[i]] <- fitdist(post_samples = object$post$H1, what = whati)
      post_summary[i, "mean"] <- mean_posterior(fit = fit[[i]])
      post_summary[i, "sd"] <- sd_posterior(fit = fit[[i]])
      post_summary[i,3:7] <- qposterior(p = p, what = whati,
                                       fit = fit[[i]],
                                       hypothesis = "H1")
    }
  } else {
    for (i in pars) {
      if (i == "psi (logor)") {
        s <- object$post$H1[["logor"]]
      } else {
        s <- object$post$H1[[i]]
      }
      post_summary[i, "mean"] <- mean(s)
      post_summary[i, "sd"] <- sd(s)
      post_summary[i,3:7] <- quantile(s, probs = p)
    }
  }

  out <- object
  out$post[["post_summary"]] <- post_summary
  out$input$digits <- digits
  out$input$raw <- raw
  class(out) <- "summary.ab"
  return(out)

}

# print summary method

#' @rdname ab-methods
#' @method print summary.ab
#' @export
print.summary.ab <- function(x, ...) {

  index <- x$input$prior_prob != 0
  hypotheses <- names(x$input$prior_prob)

  cat("Bayesian A/B Test Summary:",
      "\n\n Bayes Factors:",
      "\n\n BF10: ", round(x$bf$bf10, x$input$digits),
      "\n BF+0: ", round(x$bf$bfplus0, x$input$digits),
      "\n BF-0: ", round(x$bf$bfminus0, x$input$digits),
      "\n\n",
      " Prior Probabilities Hypotheses:",
      "\n\n ",
      paste(hypotheses[index], round(x$input$prior_prob[index], x$input$digits),
            sep = ": ", collapse = "\n "),
      "\n\n",
      " Posterior Probabilities Hypotheses:",
      "\n\n ",
      paste(hypotheses[index], round(x$post_prob[index], x$input$digits),
            sep = ": ", collapse = "\n "),
      "\n\n ",
      paste0("Posterior Summary under H1 (based on ",
             x$input$nsamples, " samples):"),
      "\n\n", sep = "")
      print(round(x$post$post_summary, x$input$digits), ...)

}

# print method

#' @rdname ab-methods
#' @method print ab
#' @export
print.ab <- function(x, ...) {

  index <- x$input$prior_prob != 0
  hypotheses <- names(x$input$prior_prob)
  cat("Bayesian A/B Test Results:",
      "\n\n Bayes Factors:",
      "\n\n BF10: ", x$bf$bf10,
      "\n BF+0: ", x$bf$bfplus0,
      "\n BF-0: ", x$bf$bfminus0,
      "\n\n",
      " Prior Probabilities Hypotheses:",
      "\n\n ",
      paste(hypotheses[index], round(x$input$prior_prob[index], 4),
            sep = ": ", collapse = "\n "),
      "\n\n",
      " Posterior Probabilities Hypotheses:",
      "\n\n ",
      paste(hypotheses[index], round(x$post_prob[index], 4),
            sep = ": ", collapse = "\n "), "\n",
      sep = "")

}

#' @rdname ab-methods
#' @method plot ab
#' @export
plot.ab <- function(x, ...) {

  if (x$input$posterior) {
    post_samples <- x$post$H1
  } else {
    x <- ab_test(data = x$input$data,
                 prior_par = x$input$prior_par,
                 prior_prob = x$input$prior_prob,
                 nsamples = x$input$nsamples,
                 is_df = x$input$is_df,
                 posterior = TRUE)
  }

  userask <- grDevices::devAskNewPage()
  grDevices::devAskNewPage(ask = TRUE)

  prob_wheel(x, type = "prior")
  prob_wheel(x, type = "posterior")

  grDevices::devAskNewPage(ask = userask)

}
