#' Methods for ab objects
#'
#' Methods defined for objects returned from the \code{\link{ab_test}} function.
#'
#' @param object,x object of class \code{ab} as returned from
#'   \code{\link{ab_test}}.
#' @param digits number of digits to print for the summary.
#' @param ... further arguments, currently ignored.
#'
#' @return The \code{summary} methods returns the \code{ab} object that is
#'   guaranteed to contain posterior samples (i.e., it adds posterior samples if
#'   they were not included already). Additionally, the list output contains the
#'   number of digits to print for the posterior summary.
#'
#'   The \code{print} methods simply print and return nothing.
#'
#'
#' @name ab-methods
NULL


# summary method

#' @rdname ab-methods
#' @method summary ab
#' @export
summary.ab <- function(object, digits = 3, ...) {

  # make sure that object contains posterior samples
  if ( ! object$input$posterior) {
    object <- ab_test(data = object$input$data,
                 prior_par = object$input$prior_par,
                 prior_prob = object$input$prior_prob,
                 nsamples = object$input$nsamples,
                 is_df = object$input$is_df,
                 posterior = TRUE)
  }

  out <- object
  out$digits <- digits

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

  # obtain summary posterior under H1
  pars <- c("psi", "beta", "logor", "or", "arisk", "rrisk", "p1", "p2")
  post_matrix <- matrix(nrow = 8, ncol = 7)
  rownames(post_matrix) <- pars
  colnames(post_matrix) <- c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%")
  p <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  fit <- vector("list", 6)
  names(fit) <- pars[1:6]
  for (i in pars[! pars %in% c("beta", "psi")]) {
    fit[[i]] <- fitdist(post_samples = x$post$H1, what = i)
    post_matrix[i, "mean"] <- mean(x$post$H1[[i]])
    post_matrix[i, "sd"] <- sd(x$post$H1[[i]])
    post_matrix[i,3:7] <- qposterior(p = p, what = i,
                                           fit = fit[[i]],
                                           hypothesis = "H1")
  }
  post_matrix["psi",] <- post_matrix["logor",]
  post_matrix["beta","mean"] <- mean(x$post$H1$beta)
  post_matrix["beta","sd"] <- sd(x$post$H1$beta)
  post_matrix["beta",3:7] <- quantile(x$post$H1$beta, probs = p)

  cat("Bayesian A/B Test Summary:",
      "\n\n Bayes Factors:",
      "\n\n BF10: ", round(x$bf$bf10, x$digits),
      "\n BF+0: ", round(x$bf$bfplus0, x$digits),
      "\n BF-0: ", round(x$bf$bfminus0, x$digits),
      "\n\n",
      " Prior Probabilities Hypotheses:",
      "\n\n ",
      paste(hypotheses[index], round(x$input$prior_prob[index], x$digits),
            sep = ": ", collapse = "\n "),
      "\n\n",
      " Posterior Probabilities Hypotheses:",
      "\n\n ",
      paste(hypotheses[index], round(x$post_prob[index], x$digits),
            sep = ": ", collapse = "\n "),
      "\n\n ",
      paste0("Posterior Summary for H1 (based on ",
             x$input$nsamples, " samples):"),
      "\n\n", sep = "")
      print(round(post_matrix, x$digits), ...)

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
