
#' Function for plotting the posterior probabilities of the hypotheses
#' sequentially.
#' @title Plot Sequential Analysis
#' @param x object of class \code{"ab"}. Note that the \code{"ab"} object needs
#'   to contain sequential data, that is, the user needs to have provided
#'   cumulative sequences of "successes"/trials.
#' @param thin allows the user to skip every \eqn{k}th data point for plotting,
#'   where the number \eqn{k} is specified via \code{thin}. For instance, in
#'   case \code{thin = 2}, only every second element of the data is displayed.
#' @param cores number of cores used for the computations.
#' @param ... further arguments
#' @details The plot shows the posterior probabilities of the hypotheses as a
#'   function of the total number of observations across the experimental and
#'   control group. On top of the plot, probability wheels (see also
#'   \code{\link{prob_wheel}}) visualize the prior probabilities of the
#'   hypotheses and the posterior probabilities of the hypotheses after taking
#'   into account all available data.
#'
#'   \strong{N.B.: This plot has been designed to look good in the following
#'   format: In inches, 530 / 72 (width) by 400 / 72 (height); in pixels, 530
#'   (width) by 400 (height).}
#'
#' @author Quentin F. Gronau
#' @example examples/example.plot_sequential.R
#' @importFrom graphics grconvertX grconvertY points text
#' @importFrom grDevices col2rgb rgb
#' @export
plot_sequential <- function(x,
                            thin = 1,
                            cores = 1,
                            ...) {

  # make sure that object is of class ab
  if ( ! inherits(x, "ab")) {
    stop("x needs to be of class 'ab'", call. = FALSE)
  }

  # plotting settings (maybe put in arguments at some point)
  lwd <- 2
  cexPoints <- 1.4
  cexAxis <- 1.2
  cexYlab <- 1.5
  cexXlab <- 1.6
  cexTextBF <- 1.4
  cexText <- 1.2
  cexLegend <- 1.2
  cexEvidence <- 1.6
  lwdAxis <- 1.2

  data <- x$input$data
  ntotal <- unique(vapply(data, length, 0))
  index <- seq(thin, ntotal, thin)

  if (index[length(index)] != ntotal) {
    warning("thinning hides the last observation",
            call. = FALSE)
  }

  nsteps <- length(index)
  n <- vapply(index, function(i) data$n1[i] + data$n2[i], 0)

  if (cores == 1) {

    out <- lapply(index,
                  FUN = compute_ab_seq,
                  ab = x,
                  data = data,
                  ntotal = ntotal)

  } else if (cores > 1) {

    if (.Platform$OS.type == "unix") {

      out <- parallel::mclapply(index, compute_ab_seq, mc.cores = cores,
                                ab = x, data = data, ntotal = ntotal)

    } else {

      cl <- parallel::makeCluster(cores)
      parallel::clusterExport(cl, c("x", "index", "data", "ntotal", "ab_test",
                                    "compute_ab_seq"))
      out <- parallel::parLapply(cl = cl, X = index, fun = compute_ab_seq,
                                 ab = x, data = data, ntotal = ntotal)
      parallel::stopCluster(cl)

    }

  }

  xticks <- pretty(c(0, data$n1[length(data$n1)] +
                       data$n2[length(data$n2)]))
  xlim <- range(xticks)
  ylim <- c(0, 1)
  yticks <- pretty(ylim)

  hyp_index <- x$input$prior_prob != 0
  #col <- c("firebrick", "grey70", "dodgerblue", "darkgoldenrod1")
  col <- c(RColorBrewer::brewer.pal(8, name = "Dark2")[2],
           RColorBrewer::brewer.pal(8, name = "Dark2")[8],
           RColorBrewer::brewer.pal(8, name = "Dark2")[1],
           RColorBrewer::brewer.pal(8, name = "Dark2")[6])
  names(col) <- c("H-", "H0", "H+", "H1")

  op <- par(mar = c(5.6, 6, 7, 7) + 0.1, las = 1, xpd = TRUE)
  plot(1, 1, xlim = xlim, ylim = ylim, ylab = "", xlab = "",
       type = "n", axes = FALSE, asp = (diff(xlim) / diff(ylim)) * 0.001 * 620,
       ...)

  axis(1, at = xticks, cex.axis = cexAxis, lwd = lwdAxis)
  axis(2, at = yticks, cex.axis = cexAxis, lwd = lwdAxis)

  mtext(text = "Posterior Probability", side = 2, las = 0,
        cex = cexYlab, line = 3.1)
  mtext("n", side = 1, cex = cexXlab, line = 2.5)

  # posterior probabilities
  for (hyp in names(hyp_index[hyp_index])) {

    values <- c(x$input$prior_prob[hyp],
                vapply(out, FUN = function(y) y$post_prob[hyp],
                       FUN.VALUE = 0))

    if (length(out) <= 60) {
      lines(c(0, n), values, col = col[hyp], lwd = 2, type = "c")
      points(c(0, n), values, pch = 21, bg = col[hyp], cex = cexPoints,
             lwd = 1.3)
    } else {
      lines(c(0, n), values, col = col[hyp], lwd = 2.7)
    }

  }

  # prior probability wheel
  offsetTopPart <- 0.03
  xx1 <- grconvertX(0.215 + 0.001 * 6, "ndc", "user")
  yy <- grconvertY(0.788 + offsetTopPart, "ndc", "user")
  radius <- grconvertX(0.2, "ndc", "user") -
    grconvertX(0.16, "ndc", "user")

  col_rgb <- col2rgb(col)
  col_trans <- rgb(t(col_rgb), alpha = max(col_rgb) / 3,
                   maxColorValue = max(col_rgb))
  names(col_trans) <- names(col)

  p_prior <- x$input$prior_prob
  p_prior <- p_prior[c("H-", "H0", "H+", "H1")]
  p_prior[p_prior == 0] <- 1e-100 # to avoid not plotting the correct colors for extreme cases
  floating.pie(xx1, yy, p_prior,
               radius = radius, col = col_trans, lwd = 2,
               lwd = 2, startpos = pi/2)
  legend_order <- c("H1", "H+", "H-", "H0")
  hypnice <- names(x$input$prior_prob)
  hypnice[hypnice == "H-"] <- "H\u2212" # to ensure equal spacing
  legend_lab <- paste0("P(", hypnice, ") = ",
                       sprintf("%.2f", round(p_prior[legend_order], 2)))
  names(legend_lab) <- names(p_prior[legend_order])

  legend_col <- col_trans[legend_order]
  legend_index <- hyp_index[legend_order]
  legend_lab <- legend_lab[names(legend_index)]
  legend(xx1 + grconvertX(0.228, "ndc", "user"), yy,
         legend = legend_lab[legend_index],
         pt.bg = legend_col[legend_index],
         pt.cex = rep(1.1, sum(legend_index)),
         pt.lwd = rep(2, sum(legend_index)),
         col = rep("black", sum(legend_index)),
         pch = rep(21, sum(legend_index)),
         bty = "n", cex = 1, xjust = 0,
         yjust = 0.5, text.col = "grey15")

  text(xx1 + grconvertX(0.135, "ndc", "user"),
       yy + grconvertY(0.315, "ndc", "user"),
       labels = "Prior Probabilities",
       cex = 1.2, pos = 4, col = "grey15")

  # posterior probability wheel
  xx2 <- grconvertX(0.495 + 0.001 * 6, "ndc", "user")
  p_post <- x$post_prob
  p_post <- p_post[c("H-", "H0", "H+", "H1")]
  p_post[p_post == 0] <- 1e-100 # to avoid not plotting the correct colors for extreme cases
  floating.pie(xx2, yy, p_post,
               radius = radius, col = col, lwd = 2,
               lwd = 2, startpos = pi/2)
  legend_lab <- paste0("P(", hypnice, " | data) = ",
                       sprintf("%.3f", round(p_post[legend_order], 3)))
  names(legend_lab) <- names(p_post[legend_order])

  legend_lab <- legend_lab[names(legend_index)]
  legend_col <- col[legend_order]
  legend(xx2 + grconvertX(0.228, "ndc", "user"), yy,
         legend = legend_lab[legend_index],
         pt.bg = legend_col[legend_index],
         pt.cex = rep(1.1, sum(legend_index)),
         pt.lwd = rep(2, sum(legend_index)),
         col = rep("black", sum(legend_index)),
         pch = rep(21, sum(legend_index)),
         bty = "n", cex = 1, xjust = 0,
         yjust = 0.5)

  text(xx2 + grconvertX(0.135, "ndc", "user"),
       yy + grconvertY(0.315, "ndc", "user"),
       labels = "Posterior Probabilities",
       cex = 1.2, pos = 4)

  # arrow
  graphics::arrows(x0 = xx2 - grconvertX(0.29, "ndc", "user"),
                   y0 = yy + grconvertY(0.315, "ndc", "user"),
                   x1 = xx2 - grconvertX(0.236, "ndc", "user"),
                   y1 = yy + grconvertY(0.315, "ndc", "user"),
                   code = 2, length = 0.1, lwd = 2)

  par(op)

}

compute_ab_seq <- function(i, ab, data, ntotal) {

  data_tmp <- list(y1 = data$y1[i],
                   n1 = data$n1[i],
                   y2 = data$y2[i],
                   n2 = data$n2[i])


  if (data_tmp$n1 == 0 || data_tmp$n2 == 0) {
    out <- list(post_prob = ab$input$prior_prob)
  } else if (i == ntotal) {
    out <- ab
  } else {
    out <- ab_test(data = data_tmp,
                  prior_par = ab$input$prior_par,
                  prior_prob = ab$input$prior_prob,
                  nsamples = ab$input$nsamples,
                  is_df = ab$input$is_df)
  }

  return(out)

}
