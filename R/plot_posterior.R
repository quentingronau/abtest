
#' Function for plotting the posterior distribution.
#' @title Plot Posterior
#' @param x object of class \code{"ab"}.
#' @param what character specifying for which quantity the posterior should be
#'   plotted. Either \code{"logor"} (i.e., log odds ratio) , \code{"or"} (i.e.,
#'   odds ratio), \code{"p1p2"} (i.e., the marginal posteriors of the latent
#'   "success" probabilities in the experimental and control condition),
#'   \code{"rrisk"} (i.e., relative risk, the ratio of the "success" probability
#'   in the experimental and the control condition), or \code{"arisk"} (i.e.,
#'   absolute risk, the difference of the "success" probability in the
#'   experimental and control condition).
#' @param hypothesis character specifying whether to plot the two-sided
#'   posterior distribution (i.e., "H1"), the one-sided posterior distribution
#'   with lower truncation point (i.e., "H+"), or the one-sided posterior
#'   distribution with upper truncation point (i.e., "H-").
#' @param ci numeric value specifying the \code{ci}\% central credible interval.
#'   The default is 0.95 which yields a 95\% central credible interval.
#' @param p1lab determines p1 x-axis label. Only relevant for \code{what =
#'   "p1p2"}.
#' @param p2lab determines p2 x-axis label. Only relevant for \code{what =
#'   "p1p2"}.
#' @param p1adj determines p1 x-axis label adjustment. Only relevant for
#'   \code{what = "p1p2"}.
#' @param p2adj determines p2 x-axis label adjustment. Only relevant for
#'   \code{what = "p1p2"}.
#' @param ... further arguments
#'
#' @details The resulting plot displays the posterior density for the quantitiy
#'   of interest and also displays the corresponding prior density. The values
#'   of the posterior median and a \code{ci}\% central credible interval are
#'   displayed on top of the plot.
#'
#' @author Quentin F. Gronau
#' @example examples/example.plot_posterior.R
#' @export
#' @importFrom graphics axis contour image lines mtext plot arrows
#' @importFrom stats dbeta pbeta var qbeta qlnorm qnorm
plot_posterior <- function(x,
                           what = "logor",
                           hypothesis = "H1",
                           ci = .95,
                           p1lab = "p1",
                           p2lab = "p2",
                           p1adj = 0.44,
                           p2adj = 0.56,
                           ...) {

  ### check arguments ###

  # check what
  if ( ! (what %in%
          c("logor", "or", "rrisk", "arisk", "p1p2"))) {
    stop('what needs to be either "logor", "or", "rrisk", "arisk",
         or "p1p2',
         call. = FALSE)
  }

  # check hypothesis
  if ( ! (hypothesis %in% c("H1", "H+", "H-"))) {
    stop('hypothesis needs to be either "H1", "H+", or "H-"',
         call. = FALSE)  }


  # plotting settings (maybe put in arguments at some point)
  lwd <- 2
  cexPoints <- 1.5
  cexAxis <- 1.2
  cexYlab <- 1.5
  cexXlab <- 1.5
  cexTextBF <- 1.4
  cexCI <- 1.1
  cexLegend <- 1.2
  lwdAxis <- 1.2

  # extract posterior samples (2-sided)
  if (x$input$posterior) {
    post_samples <- x$post$H1
  } else {
    x <- ab_test(data = x$input$data, prior_par = x$input$prior_par,
                  prior_prob = x$input$prior_prob,
                  nsamples = x$input$nsamples, is_df = x$input$is_df,
                  posterior = TRUE)
  }

  # set limits plot
  if (what == "p1p2") {

    post_samples <- switch(hypothesis,
      "H1" = x$post$H1,
      "H+" = x$post$Hplus,
      "H-" = x$post$Hminus
    )

    # fit distributions
    fit1 <- fitdist(post_samples = post_samples, what = "p1")
    fit2 <- fitdist(post_samples = post_samples, what = "p2")

    # compute quantiles
    pci <- c((1 - ci) / 2, 1 - (1 - ci) / 2)
    CIlow1 <- qposterior(p = pci[1], what = "p1", hypothesis = hypothesis,
                         fit = fit1)
    CIhigh1 <- qposterior(p = pci[2], what = "p1", hypothesis = hypothesis,
                          fit = fit1)
    medianPosterior1 <- qposterior(p = .5, what = "p1",
                                   hypothesis = hypothesis,
                                   fit = fit1)
    CIlow2 <- qposterior(p = pci[1], what = "p2", hypothesis = hypothesis,
                         fit = fit2)
    CIhigh2 <- qposterior(p = pci[2], what = "p2", hypothesis = hypothesis,
                          fit = fit2)
    medianPosterior2 <- qposterior(p = .5, what = "p2",
                                   hypothesis = hypothesis,
                                   fit = fit2)

    if (abs(CIhigh1 - CIlow1) < .01 && abs(CIhigh2 - CIlow2) < .01) {
      qlow1 <- qposterior(p = .001, what = "p1", hypothesis = hypothesis,
                          fit = fit1)
      qhigh1 <- qposterior(p = .999, what = "p1", hypothesis = hypothesis,
                           fit = fit1)
      qlow2 <- qposterior(p = .001, what = "p2", hypothesis = hypothesis,
                          fit = fit2)
      qhigh2 <- qposterior(p = .999, what = "p2", hypothesis = hypothesis,
                           fit = fit2)
      xlim <- c(min(CIlow1, CIlow2, qlow1, qlow2),
                max(CIhigh1, CIhigh2, qhigh1, qhigh2))
    } else {
      xlim <- c(0, 1)
    }

    xticks <- pretty(xlim, n = 5)
    xlim <- range(xticks)

    # prior and posterior density
    stretch <- ifelse(hypothesis == "H1", 1.2, 1.32)
    xx <- seq(min(xticks), max(xticks), length.out = 1e3)
    priorLine1 <- dprior(x1 = xx, x2 = NULL, prior_par = x$input$prior_par,
                         what = "p1", hypothesis = hypothesis)
    priorLine2 <- dprior(x1 = xx, x2 = NULL, prior_par = x$input$prior_par,
                         what = "p2", hypothesis = hypothesis)
    postLine1 <- dposterior(x = xx, what = "p1", hypothesis = hypothesis,
                            fit = fit1)
    postLine2 <- dposterior(x = xx, what = "p2", hypothesis = hypothesis,
                            fit = fit2)

    dpriormax1 <- max(priorLine1)
    dpriormax2 <- max(priorLine2)
    dpostmax1 <- max(postLine1)
    dpostmax2 <- max(postLine2)

    ylim <- c(0, stretch * max(c(dpriormax1, dpriormax2,
                                 dpostmax1, dpostmax2)))
    yticks <- pretty(ylim)
    yticks <- c(yticks,
                yticks[length(yticks)] +
                  diff(yticks[(length(yticks)-1):length(yticks)]),
                yticks[length(yticks)] +
                  2*diff(yticks[(length(yticks)-1):length(yticks)]))
    ylim <- range(yticks)

  } else {

    post_samples <- x$post$H1

    # fit distribution
    fit <- fitdist(post_samples = post_samples, what = what)

    prior_range <- qprior(p = c(.05, .95), prior_par = x$input$prior_par,
                          what = what, hypothesis = "H1")
    post_range <- qposterior(p = c(.001, .999), what = what,
                             hypothesis = "H1",
                             fit = fit)
    xticks <- pretty(c(prior_range, post_range))
    xlim <- range(xticks)
    stretch <- ifelse(hypothesis == "H1", 1.2, 1.32)

    xx <- seq(min(xticks), max(xticks), length.out = 1e3)
    priorLine <- dprior(x1 = xx, x2 = NULL, prior_par = x$input$prior_par,
                        what = what, hypothesis = hypothesis)
    postLine <- dposterior(x = xx, what = what, hypothesis = hypothesis,
                           fit = fit)

    dpriormax <- max(priorLine)
    dpostmax <- max(postLine)

    # compute 95% credible interval & median
    pci <- c((1 - ci) / 2, 1 - (1 - ci) / 2)
    CIlow <- qposterior(p = pci[1], what = what, hypothesis = hypothesis,
                        fit = fit)
    CIhigh <- qposterior(p = pci[2], what = what, hypothesis = hypothesis,
                         fit = fit)
    ylim <- c(0, max(stretch * c(dpriormax, dpostmax)))
    yticks <- pretty(ylim)
    xlim <- c(min(CIlow, range(xticks)[1]), max(range(xticks)[2], CIhigh))
    plot(1, 1, xlim = xlim, ylim = ylim, ylab = "", xlab = "",
         type = "n", axes = FALSE)
    yCI <- grconvertY(max(dpriormax, dpostmax), "user", "ndc") + 0.04
    yCI <- grconvertY(yCI, "ndc", "user")
    medianPosterior <- qposterior(p = .5, what = what,
                                  hypothesis = hypothesis,
                                  fit = fit)

    if ((yCI > yticks[length(yticks) - 2] &&
        diff(yticks[(length(yticks)-1):length(yticks)]) /
        diff(range(yticks)) < 1 / 6) || (yCI > yticks[length(yticks) - 1])) {
      yticks <- c(yticks,
                  yticks[length(yticks)] +
                    diff(yticks[(length(yticks)-1):length(yticks)]))
      ylim <- range(yticks)
    }

  }

  yticks <- pretty(ylim)
  ylim <- range(yticks)
  ylabels <- as.character(yticks)

  xlab <- switch(what,
                 "logor" = "Log Odds Ratio",
                 "or" = "Odds Ratio",
                 "rrisk" = "Relative Risk",
                 "arisk" = "Absolute Risk")

  op <- par(mar = c(5.6, 5, 7, 4) + 0.1, las = 1, xpd = TRUE)
  plot(1, 1, xlim = xlim, ylim = ylim, ylab = "", xlab = "",
       type = "n", axes = FALSE)

  if (what == "p1p2") {

    # transparent colors
    col_trans1 <- rgb(0, 0, 0, alpha = 0.5)
    col_trans2 <- rgb(0, 0, 0, alpha = 0.8)

    lines(xx, postLine1, lwd = lwd, col = col_trans1)
    lines(xx, postLine2, lwd = lwd, col = col_trans2)
    lines(xx, priorLine1, lwd = lwd, lty = 3, col = col_trans1)
    lines(xx, priorLine2, lwd = lwd, lty = 3, col = col_trans2)

    # credible interval
    yCI1 <- grconvertY(max(dpriormax1, dpriormax2, dpostmax1, dpostmax2),
                       "user", "ndc") + 0.03
    yCI1 <- grconvertY(yCI1, "ndc", "user")
    yCI2 <- grconvertY(max(dpriormax1, dpriormax2, dpostmax1, dpostmax2),
                       "user", "ndc") + 0.06
    yCI2 <- grconvertY(yCI2, "ndc", "user")
    arrows(CIlow1, yCI1 , CIhigh1, yCI1, angle = 90, code = 3,
           length = 0.1, lwd = lwd, col = col_trans1)
    arrows(CIlow2, yCI2 , CIhigh2, yCI2, angle = 90, code = 3,
           length = 0.1, lwd = lwd, col = col_trans2)

    medianText1 <- formatC(medianPosterior1, digits = 3, format = "f")
    medianText2 <- formatC(medianPosterior2, digits = 3, format = "f")

    offsetTopPart <- 0.06
    yy <- grconvertY(0.756 + offsetTopPart, "ndc", "user")
    yy2 <- grconvertY(0.812 + offsetTopPart, "ndc", "user")

    textci <- as.character(ci)
    percentage <- strsplit(textci, split = "0.")[[1]][2]
    CIText1 <- paste0(percentage, "% CI: [",
                     bquote(.(formatC(CIlow1, 3, format = "f"))), ", ",
                     bquote(.(formatC(CIhigh1, 3, format = "f"))), "]")
    CIText2 <- paste0(percentage, "% CI: [",
                      bquote(.(formatC(CIlow2, 3, format = "f"))), ", ",
                      bquote(.(formatC(CIhigh2, 3, format = "f"))), "]")
    medianLegendText1 <- paste("median =", medianText1)
    medianLegendText2 <- paste("median =", medianText2)

    text(min(xticks), yy2, medianLegendText1, cex = 1.1,
         pos = 4, col = col_trans1)
    text(min(xticks), yy, CIText1, cex = 1.1, pos = 4, col = col_trans1)
    text(max(xticks), yy2, medianLegendText2, cex = 1.1,
         pos = 2, col = col_trans2)
    text(max(xticks), yy, CIText2, cex = 1.1, pos = 2, col = col_trans2)
    index <- which.min(c(medianPosterior1, medianPosterior2,
                         1 - medianPosterior1, 1 - medianPosterior2))
    legendPosition <- switch(index,
                             "1" = max(xticks),
                             "2" = max(xticks),
                             "3" = min(xticks),
                             "4" = min(xticks))
    legend(legendPosition, max(yticks),
           legend = c("Posterior", "Prior"), lty = c(1, 3), bty = "n",
           lwd = rep(lwd, 2), cex = cexLegend,
           xjust = ifelse(all.equal(legendPosition, max(xticks)), 1, 0),
           yjust = 1, x.intersp = .6, seg.len = 1.2)
    mtext(c(p1lab, "&", p2lab), col = c(col_trans1, "black", col_trans2),
          adj = c(p1adj, 0.5, p2adj), side = 1, cex = cexXlab, line = 2.5)

  } else {

    lines(xx, postLine, lwd = lwd)
    lines(xx, priorLine, lwd = lwd, lty = 3)

    # credible interval
    arrows(CIlow, yCI, CIhigh, yCI, angle = 90, code = 3,
           length = 0.1, lwd = lwd)

    medianText <- formatC(medianPosterior, digits = 3, format = "f")

    offsetTopPart <- 0.06
    yy <- grconvertY(0.756 + offsetTopPart, "ndc", "user")
    yy2 <- grconvertY(0.812 + offsetTopPart, "ndc", "user")

    textci <- as.character(ci)
    percentage <- strsplit(textci, split = "0.")[[1]][2]
    CIText <- paste0(percentage, "% CI: [",
                     bquote(.(formatC(CIlow, 3, format = "f"))), ", ",
                     bquote(.(formatC(CIhigh, 3, format = "f"))), "]")
    medianLegendText <- paste("median =", medianText)

    text(max(xticks), yy2, medianLegendText, cex = 1.1, pos = 2)
    text(max(xticks), yy, CIText, cex = 1.1, pos = 2)

    if (medianPosterior <= mean(xlim) ) {
      legendPosition <- max(xticks)
      legend(legendPosition, max(yticks), legend = c("Posterior", "Prior"),
             lty = c(1, 3), bty = "n", lwd = rep(lwd, 2), cex = cexLegend,
             xjust = 1, yjust = 1, x.intersp = .6, seg.len = 1.2)
    } else {
      legendPosition <- min(xticks)
      legend(legendPosition, max(yticks), legend = c("Posterior", "Prior"),
             lty = c(1, 3), bty = "n", lwd = rep(lwd, 2), cex = cexLegend,
             xjust = 0, yjust = 1, x.intersp = .6, seg.len = 1.2)
    }

    mtext(xlab, side = 1, cex = cexXlab, line = 2.5)

  }

  axis(1, at = xticks, cex.axis = cexAxis, lwd = lwdAxis)
  axis(2, at = yticks, labels = ylabels, cex.axis = cexAxis, lwd = lwdAxis)

  if (nchar(ylabels[length(ylabels)]) > 4) {
    mtext(text = "Density", side = 2, las = 0, cex = cexYlab, line = 4)
  } else if (nchar(ylabels[length(ylabels)]) == 4) {
    mtext(text = "Density", side = 2, las = 0, cex = cexYlab, line = 3.25)
  } else if (nchar(ylabels[length(ylabels)]) < 4) {
    mtext(text = "Density", side = 2, las = 0, cex = cexYlab, line = 2.85)
  }

  par(op)

}

fitdist <- function(post_samples, what) {

  # assumes that samples come from H1

  samples <- post_samples[[what]]

  if (what == "arisk") {
    # scale back to [0,1] interval for fitting beta distribution
    samples <- (samples + 1) / 2
  }

  # compute mean and variance for method of moments
  m <- mean(samples)
  v <- var(samples)

  if (what %in% c("logor", "beta")) {
    family <- "normal"
    pars <- list(mean = m, sd = sqrt(v))
  } else if (what %in% c("or", "rrisk")) {
    family <- "log-normal"
    pars <- list(meanlog = log(m / sqrt(1 + v / m^2)),
                 sdlog = sqrt(log(1 + v / m^2)))
  } else if (what %in% c("arisk", "p1", "p2")) {
    family <- switch(what,
                     "arisk" = "scaled-beta",
                     "p1" = "beta",
                     "p2" = "beta")
    m <- mean(samples)
    v <- var(samples)
    pars <- list(alpha = m * (m * (1 - m) / v - 1),
                 beta = (1 - m) * (m * (1 - m) / v - 1))
  }

  out <- list(what = what, family = family, pars = pars)
  return(out)

}

mean_posterior <- function(fit) {

  what <- fit$what
  pars <- fit$pars

  if (what %in% c("logor", "beta")) {
    m <- pars[["mean"]]
  } else if (what %in% c("or", "rrisk")) {
    m <- exp(pars[["meanlog"]] + pars[["sdlog"]]^2 / 2)
  } else if (what == "arisk") {
    mtmp <- pars[["alpha"]] / (pars[["alpha"]] + pars[["beta"]])
    m <- 2 * mtmp - 1
  } else if (what %in% c("p1", "p2")) {
    m <- pars[["alpha"]] / (pars[["alpha"]] + pars[["beta"]])
  }

  return(m)

}

sd_posterior <- function(fit) {

  what <- fit$what
  pars <- fit$pars

  if (what %in% c("logor", "beta")) {
    s <- pars[["sd"]]
  } else if (what %in% c("or", "rrisk")) {
    s2 <- (exp(pars[["sdlog"]]^2) - 1) *
      exp(2 * pars[["meanlog"]] + pars[["sdlog"]]^2)
    s <- sqrt(s2)
  } else if (what == "arisk") {
    s2tmp <- pars[["alpha"]] * pars[["beta"]]  /
      ((pars[["alpha"]] + pars[["beta"]])^2 *
         (pars[["alpha"]] + pars[["beta"]] + 1))
    s <- 2 * sqrt(s2tmp)
  } else if (what %in% c("p1", "p2")) {
    s2 <- pars[["alpha"]] * pars[["beta"]]  /
      ((pars[["alpha"]] + pars[["beta"]])^2 *
         (pars[["alpha"]] + pars[["beta"]] + 1))
    s <- sqrt(s2)
  }

  return(s)

}

dposterior <- function(x, what, hypothesis, fit, log = FALSE) {

  if (what == "logor") {

    out <- dnorm(x, mean = fit$pars$mean, sd = fit$pars$sd, log = TRUE)
    if (hypothesis == "H+") {
      out[x < 0] <- -Inf
    } else if (hypothesis == "H-") {
      out[x > 0] <- -Inf
    }
    logconst <- switch(hypothesis,
                       "H1" = 0,
                       "H+" = pnorm(0, mean = fit$pars$mean,
                                    sd = fit$pars$sd,
                                    lower.tail = FALSE,
                                    log.p = TRUE),
                       "H-" = pnorm(0, mean = fit$pars$mean,
                                    sd = fit$pars$sd,
                                    lower.tail = TRUE,
                                    log.p = TRUE))

  } else if (what %in% c("or", "rrisk")) {

    out <- dlnorm(x, meanlog = fit$pars$meanlog,
                  sdlog = fit$pars$sdlog, log = TRUE)
    if (hypothesis == "H+") {
      out[x < 1] <- -Inf
    } else if (hypothesis == "H-") {
      out[x > 1] <- -Inf
    }
    logconst <- switch(hypothesis,
                       "H1" = 0,
                       "H+" = plnorm(1, meanlog = fit$pars$meanlog,
                                     sdlog = fit$pars$sdlog,
                                     lower.tail = FALSE,
                                     log.p = TRUE),
                       "H-" = plnorm(1, meanlog = fit$pars$meanlog,
                                     sdlog = fit$pars$sdlog,
                                     lower.tail = TRUE,
                                     log.p = TRUE))

  } else if (what == "arisk") {

    out <- log(.5) + dbeta( (x + 1) / 2, shape1 = fit$pars$alpha,
                           shape2 = fit$pars$beta, log = TRUE)
    if (hypothesis == "H+") {
      out[x < 0] <- -Inf
    } else if (hypothesis == "H-") {
      out[x > 0] <- -Inf
    }
    logconst <- switch(hypothesis,
                       "H1" = 0,
                       "H+" = pbeta(.5, shape1 = fit$pars$alpha,
                                    shape2 = fit$pars$beta,
                                    lower.tail = FALSE,
                                    log.p = TRUE),
                       "H-" = pbeta(.5, shape1 = fit$pars$alpha,
                                    shape2 = fit$pars$beta,
                                    lower.tail = TRUE,
                                    log.p = TRUE))

  } else if (what %in% c("p1", "p2")) {

    out <- dbeta(x, shape1 = fit$pars$alpha, shape2 = fit$pars$beta,
                 log = TRUE)
    logconst <- 0

  }

  out <- out - logconst

  if ( ! log) {
    out <- exp(out)
  }

  return(out)

}

qposterior <- function(p, what, hypothesis, fit) {

  if (what == "logor") {

    if (hypothesis == "H1") {
      out <- qnorm(p, mean = fit$pars$mean, sd = fit$pars$sd)
    } else if (hypothesis == "H+") {
      leftprob <- pnorm(0, mean = fit$pars$mean, sd = fit$pars$sd)
      out <- qnorm(leftprob + p * (1 - leftprob),
                   mean = fit$pars$mean, sd = fit$pars$sd)
    } else if (hypothesis == "H-") {
      leftprob <- pnorm(0, mean = fit$pars$mean, sd = fit$pars$sd)
      out <- qnorm(p * leftprob,
                   mean = fit$pars$mean, sd = fit$pars$sd)
    }

  } else if (what == "beta") {

    out <- qnorm(p, mean = fit$pars$mean, sd = fit$pars$sd)

  } else if (what %in% c("or", "rrisk")) {

    if (hypothesis == "H1") {
      out <- qlnorm(p, meanlog = fit$pars$meanlog,
                    sdlog = fit$pars$sdlog)
    } else if (hypothesis == "H+") {
      leftprob <- plnorm(1, meanlog = fit$pars$meanlog,
                         sdlog = fit$pars$sdlog)
      out <- qlnorm(leftprob + p * (1 - leftprob),
                    meanlog = fit$pars$meanlog,
                    sdlog = fit$pars$sdlog)
    } else if (hypothesis == "H-") {
      leftprob <- plnorm(1, meanlog = fit$pars$meanlog,
                         sdlog = fit$pars$sdlog)
      out <- qlnorm(p * leftprob,
                    meanlog = fit$pars$meanlog,
                    sdlog = fit$pars$sdlog)
    }

  } else if (what == "arisk") {

    if (hypothesis == "H1") {
      out <- qbeta(p, shape1 = fit$pars$alpha,
                   shape2 = fit$pars$beta)
    } else if (hypothesis == "H+") {
      leftprob <- pbeta(.5, shape1 = fit$pars$alpha,
                        shape2 = fit$pars$beta)
      out <- qbeta(leftprob + p * (1 - leftprob),
                   shape1 = fit$pars$alpha,
                   shape2 = fit$pars$beta)
    } else if (hypothesis == "H-") {
      leftprob <- pbeta(.5, shape1 = fit$pars$alpha,
                        shape2 = fit$pars$beta)
      out <- qbeta(p * leftprob,
                   shape1 = fit$pars$alpha,
                   shape2 = fit$pars$beta)
    }

    # scale back to [-1, 1] interval
    out <- 2 * out - 1

  } else if (what %in% c("p1", "p2")) {

    out <- qbeta(p,shape1 = fit$pars$alpha, shape2 = fit$pars$beta)

  }

  return(out)

}
