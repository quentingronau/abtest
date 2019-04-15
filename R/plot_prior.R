
#' Function for plotting parameter prior distributions.
#' @title Plot Prior
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
#' @param what character specifying for which quantity the prior should be
#'   plotted. Either \code{"logor"} (i.e., log odds ratio) , \code{"or"} (i.e.,
#'   odds ratio), \code{"p1p2"} (i.e., plots the joint distribution of the
#'   latent "success" probability in the experimental and control condition),
#'   \code{"p1"} (i.e., latent "success" probability in the control condition),
#'   \code{"p2"} (i.e., latent "success" probability in the experimental
#'   condition), \code{"p2givenp1"} (i.e., plots the conditional distribution of
#'   the latent "success" probability in the experimental condition given a
#'   "success" probability of \code{p1} in the control condition),
#'   \code{"rrisk"} (i.e., relative risk, the ratio of the "success" probability
#'   in the experimental and the control condition), or \code{"arisk"} (i.e.,
#'   absolute risk, the difference of the "success" probability in the
#'   experimental and control condition).
#' @param hypothesis character specifying whether to plot a two-sided prior
#'   (i.e., "H1"), a one-sided prior with lower truncation point (i.e., "H+"),
#'   or a one-sided prior with upper truncation point (i.e., "H-").
#' @param p1 value of the "success" probability in the control condtion. Only
#'   used when \code{what = "p2givenp1"}.
#' @param ... further arguments.
#'
#' @note Internally, the test-relevant prior is always a normal prior on the log
#'   odds ratio, however, the \code{plot_prior} function also allows one to plot
#'   the implied prior on different quantities.
#'
#' @author Quentin F. Gronau
#' @example examples/example.plot_prior.R
#' @export
#' @importFrom graphics axis contour image lines mtext plot
plot_prior <- function(prior_par = list(mu_psi = 0, sigma_psi = 1,
                                        mu_beta = 0, sigma_beta = 1),
                       what = "logor",
                       hypothesis = "H1",
                       p1 = 0.5,
                       ...) {

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

  # check p1
  if (p1 < 0 || p1 > 1) {
    stop('p1 needs to be between 0 and 1',
         call. = FALSE)
  }

  # plotting settings (maybe put in arguments at some point)
  cexAxis <- 1.2
  cexYlab <- 1.5
  cexXlab <- 1.5

  if (what == "logor") {
    xlab <- "Log Odds Ratio"
  } else if (what == "or") {
    xlab <- "Odds Ratio"
  } else if (what == "rrisk") {
    xlab <- "Relative Risk"
  } else if (what == "arisk") {
    xlab <- "Absolute Risk"
  } else if (what == "p2givenp1") {

    bounds <- switch(hypothesis,
                     "H1" =  c(0, 1),
                     "H+" =  c(p1, 1),
                     "H-" =  c(0, p1))

    xticks <- pretty(bounds)
    xlim <- range(xticks)
    xlab <- paste0("p2 given p1 = ", round(p1, 3))

  } else if (what == "p1") {

    bounds <- c(0, 1)

    xticks <- pretty(bounds)
    xlim <- range(xticks)
    xlab <- "p1"

  } else if (what == "p2") {

    bounds <- c(0, 1)

    xticks <- pretty(bounds)
    xlim <- range(xticks)
    xlab <- "p2"

  } else if (what == "p1p2") {

    xlab <- "p1"
    ylab <- "p2"
    xlim <- ylim <- c(0, 1)

  }

  if (what %in% c("logor", "or", "rrisk", "arisk")) {

    plo <- switch(what,
                  "logor" = 0.001,
                  "or" = 0.02,
                  "rrisk" = 0.02,
                  "arisk" = 0.001)
    pup <- 1 - plo
    olo <- qprior(p = plo, prior_par = prior_par,
                  what = what, hypothesis = hypothesis)
    oup <- qprior(p = pup, prior_par = prior_par,
                  what = what, hypothesis = hypothesis)
    xlim <- c(olo, oup)
    xticks <- pretty(xlim)
    xlim <- range(xticks)

  }

  # transparent color
  col_trans <- rgb(0, 0, 0, alpha = 0.15)

  if (what != "p1p2") {

    xx <- seq(xlim[1], xlim[2], length.out = 500)
    yy <- dprior(x1 = xx, x2 = p1, prior_par = prior_par,
                 hypothesis = hypothesis, what = what)
    yticks <- pretty(yy)
    ylim <- range(yticks)

    plot(1, axes = FALSE, type = "n", xlim = xlim, ylim = ylim,
         xlab = "", ylab = "", ...)
    lines(xx, yy, lwd = 2, ...)
    graphics::polygon(c(xx, xx[length(xx)], xx[1]),
                      y = c(yy, rep(ylim[1], 2)),
                      col = col_trans)
    axis(1, at = xticks, cex.axis = cexAxis,  ...)
    axis(2, at = yticks, cex.axis = cexAxis, las = 1, ...)
    mtext(text = xlab, side = 1, line = 2.5, cex = cexXlab, ...)
    mtext(text = "Density", side = 2, line = 2.8, cex = cexYlab, ...)

  } else if (what == "p1p2") {

    xx <- seq(xlim[1], xlim[2], length.out = 100)
    yy <- xx
    zz <- outer(xx, yy, dprior, prior_par = prior_par,
                hypothesis = hypothesis, what = what)
    image(xx, yy, zz, xlab = "", ylab = "", las = 1,
          col = grDevices::heat.colors(50),
          cex.axis = cexAxis, ...)
    mtext(text = "p1", side = 1, line = 2.5, cex = cexXlab, ...)
    mtext(text = "p2", side = 2, line = 2.6, cex = cexYlab, las = 1, ...)
    contour(xx, yy, zz, add = TRUE, col = "darkblue", ...)

  }

}
