
#' Function for plotting Bayes factor robustness check results (i.e., prior
#' sensitivity analysis).
#' @title Plot Bayes Factor Robustness Check
#' @param x object of class \code{"ab"}.
#' @param bftype character that specifies which Bayes factor is plotted. Either
#'   \code{"BF10"}, \code{"BF01"}, \code{"BF+0"}, \code{"BF0+"}, \code{"BF-0"},
#'   or \code{"BF0-"}.
#' @param log Boolean that specifies whether the log Bayes factor is plotted.
#' @param mu_range numeric vector of length two that specifies the range of
#'   \code{mu_psi} values to consider.
#' @param sigma_range numeric vector of length two that specifies the range of
#'   \code{sigma_psi} values to consider.
#' @param mu_steps numeric value that specifies in how many discrete steps the
#'   interval \code{mu_range} is partitioned.
#' @param sigma_steps numeric value that specifies in how many discrete steps
#'   the interval \code{sigma_range} is partitioned.
#' @param cores number of cores used for the computations.
#' @param ... further arguments passed to \code{filled.contour}.
#' @details The plot shows how the Bayes factor changes as a function of the
#'   normal prior location parameter \code{mu_psi} and the normal prior scale
#'   parameter \code{sigma_psi} (i.e., a prior sensitivity analysis with respect
#'   to the normal prior on the test-relevant log odds ratio).
#'
#' @return Returns a \code{data.frame} with the \code{mu_psi} values,
#'   \code{sigma_psi} values, and corresponding (log) Bayes factors.
#'
#' @example examples/example.plot_robustness.R
#'
#' @author Quentin F. Gronau
#' @export
plot_robustness <- function(x,
                            bftype = "BF10",
                            log = FALSE,
                            mu_range = c(0, 0.3),
                            sigma_range = c(0.25, 1),
                            mu_steps = 40,
                            sigma_steps = 40,
                            cores = 1,
                            ...) {

  # make sure that object is of class ab
  if ( ! inherits(x, "ab")) {
    stop("x needs to be of class 'ab'", call. = FALSE)
  }

  # check bftype
  if ( ! bftype %in% c("BF10", "BF01", "BF+0", "BF0+", "BF-0", "BF0-")) {
    stop("bftype needs to be either 'BF10', 'BF01', 'BF+0', 'BF0+', 'BF-0', or 'BF0-'",
         call. = FALSE)
  }

  # check that sigma_range is positive
  if (any(sigma_range <= 0)) {
    stop("sigma_range may not contain values smaller or equal to 0",
         call. = FALSE)
  }

  # conduct robustness check
  mu <- seq(mu_range[1], mu_range[2], length.out = mu_steps)
  sigma <- seq(sigma_range[1], sigma_range[2], length.out = sigma_steps)

  prior_par_matrix <- expand.grid(mu, sigma)
  colnames(prior_par_matrix) <- c("mu_psi", "sigma_psi")

  if (cores == 1) {

    logbfs <- apply(prior_par_matrix, 1, compute_logbf, ab = x)

  } else if (cores > 1) {

    l <- lapply(seq_len(nrow(prior_par_matrix)),
                function(i) prior_par_matrix[i,])

    if (.Platform$OS.type == "unix") {

      logbfs <- parallel::mclapply(l, compute_logbf, mc.cores = cores, ab = x)

    } else {

      cl <- parallel::makeCluster(cores)
      parallel::clusterExport(cl, c("x", "l", "ab_test", "compute_logbf"))
      logbfs <- parallel::parLapply(cl = cl, X = l, fun = compute_logbf, ab = x)
      parallel::stopCluster(cl)

    }

  }

  if (bftype %in% c("BF10", "BF+0", "BF-0")) {

    bfname <- switch(bftype,
                     "BF10" = "bf10",
                     "BF+0" = "bfplus0",
                     "BF-0" = "bfminus0")
    bf <- vapply(logbfs, function(y) y[[bfname]], 1)

  } else if (bftype %in% c("BF01", "BF0+", "BF0-")) {

    bfname <- switch(bftype,
                     "BF01" = "bf10",
                     "BF0+" = "bfplus0",
                     "BF0-" = "bfminus0")
    bf <- 1 / exp(vapply(logbfs, function(y) y[[bfname]], 1))
    bf <- log(bf)

  }

  subscripts <- strsplit(bftype, split = "")[[1]][3:4]
  if (log) {
    key.title <- bquote("Log(" ~BF[.(subscripts[1])][.(subscripts[2])]~")")
  } else {
    bf <- exp(bf)
    key.title <- bquote(BF[.(subscripts[1])][.(subscripts[2])])
  }

  bf_matrix <- matrix(bf, nrow = mu_steps, ncol = sigma_steps, byrow = FALSE)

  graphics::filled.contour(
    mu,
    sigma,
    bf_matrix,
    key.title = graphics::title(main = key.title, xpd = NA),
    las = 1,
    xlab = expression(mu[psi]),
    cex.axis = 1.3, cex.lab = 1.4,
    ...)
  mtext(expression(sigma[psi]),
        side = 2,
        line = 2.5,
        cex = 1.4,
        las = 1)

  out <- cbind(prior_par_matrix, bf)
  colnames(out) <- c("mu_psi", "sigma_psi", ifelse(log, "logbf", "bf"))
  out <- as.data.frame(out)

  invisible(out)

}

compute_logbf <- function(y, ab) {

  prior_par_i <- list(mu_psi = y[["mu_psi"]],
                      sigma_psi = y[["sigma_psi"]],
                      mu_beta = ab$input$prior_par$mu_beta,
                      sigma_beta = ab$input$prior_par$sigma_beta)
  ab_test(data = ab$input$data,
          prior_par = prior_par_i,
          prior_prob = ab$input$prior_prob,
          nsamples = ab$input$nsamples,
          is_df = ab$input$is_df,
          posterior = FALSE)$logbf

}
