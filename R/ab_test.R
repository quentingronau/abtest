
#--------------------------------------------------------------------------
# Models (for details, see Kass & Vaidyanathan, 1992)
#--------------------------------------------------------------------------

# H0:
#
# y1 ~ Binomial(n1, p)
# y2 ~ Binomial(n2, p)
# log(p / (1 - p)) = beta
# beta ~ N(mu_beta, sigma_beta^2)

# H1:
#
# y1 ~ Binomial(n1, p1)
# y2 ~ Binomial(n2, p2)
# log(p1 / (1 - p1)) = beta - psi / 2
# log(p2 / (1 - p2)) = beta + psi / 2
# beta ~ N(mu_beta, sigma_beta^2)
# psi ~ N(mu_psi, sigma_psi^2)

# H+:
#
# y1 ~ Binomial(n1, p1)
# y2 ~ Binomial(n2, p2)
# log(p1 / (1 - p1)) = beta - psi / 2
# log(p2 / (1 - p2)) = beta + psi / 2
# beta ~ N(mu_beta, sigma_beta^2)
# psi ~ N(mu_psi, sigma_psi^2)T(0,)

# H-:
#
# y1 ~ Binomial(n1, p1)
# y2 ~ Binomial(n2, p2)
# log(p1 / (1 - p1)) = beta - psi / 2
# log(p2 / (1 - p2)) = beta + psi / 2
# beta ~ N(mu_beta, sigma_beta^2)
# psi ~ N(mu_psi, sigma_psi^2)T(,0)


#' Function for conducting a Bayesian A/B test (i.e., test between two
#' proportions).
#' @useDynLib abtest
#' @importFrom Rcpp sourceCpp
#' @title Bayesian A/B Test
#' @param data list or data frame with the data. This list (data frame) needs to
#'   contain the following elements: \code{y1} (number of "successes" in the
#'   control condition), \code{n1} (number of trials in the control condition),
#'   \code{y2} (number of "successes" in the experimental condition), \code{n2}
#'   (number of trials in the experimental condition). Each of these elements
#'   needs to be an integer. Alternatively, the user can provide for each of the
#'   elements a vector with a cumulative sequence of "successes"/trials. This
#'   allows the user to produce a sequential plot of the posterior probabilities
#'   for each hypothesis by passing the result object of class \code{"ab"} to
#'   the \code{\link{plot_sequential}} function.
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
#' @param prior_prob named vector with prior probabilities for the four
#'   hypotheses \code{"H1"}, \code{"H+"}, \code{"H-"}, and \code{"H0"}.
#'   \code{"H1"} states that the "success" probability differs between the
#'   control and the experimental condition but does not specify which one is
#'   higher. \code{"H+"} states that the "success" proability in the
#'   experimental condition is higher than in the control condition, \code{"H-"}
#'   states that the "success" probability in the experimental condition is
#'   lower than in the control condition. \code{"H0"} states that the "success"
#'   probability is identical (i.e., there is no effect). The one-sided
#'   hypotheses \code{"H+"} and \code{"H-"} are obtained by truncating the
#'   normal prior on the log odds ratio so that it assigns prior mass only to
#'   the allowed log odds ratio values (e.g., for \code{"H+"} a normal prior
#'   that is truncated from below at 0). If \code{NULL} (default) the prior
#'   probabilities are set to \code{c(0, 1/4, 1/4, 1/2)}. That is, the default
#'   assigns prior probability .5 to the hypothesis that there is no effect
#'   (i.e., \code{"H0"}). The remaining prior probability (i.e., also .5) is
#'   split evenly across the hypothesis that there is a positive effect (i.e.,
#'   \code{"H+"}) and the hypothesis that there is a negative effect (i.e.,
#'   \code{"H-"}).
#' @param nsamples determines the number of importance samples for obtaining the
#'   log marginal likelihood for \code{"H+"} and \code{"H-"} and the number of
#'   posterior samples in case \code{posterior = TRUE}. The default is
#'   \code{10000}.
#' @param is_df degrees of freedom of the multivariate t importance sampling
#'   proposal density. The default is \code{5}.
#' @param posterior Boolean which indicates whether posterior samples should be
#'   returned. The default is \code{FALSE}.
#' @details The implemented Bayesian A/B test is based on the following model by
#'   Kass and Vaidyanathan (1992, section 3): \deqn{log(p1/(1 - p1)) = \beta -
#'   \psi/2} \deqn{log(p2/(1 - p2)) = \beta + \psi/2} \deqn{y1 ~ Binomial(n1,
#'   p1)} \deqn{y2 ~ Binomial(n2, p2).} \code{"H0"} states that \eqn{\psi = 0},
#'   \code{"H1"} states that \eqn{\psi != 0}, \code{"H+"} states that \eqn{\psi
#'   > 0}, and \code{"H-"} states that \eqn{\psi < 0}. Normal priors are
#'   assigned to the two parameters \eqn{\psi} (i.e., the test-relevant log odds
#'   ratio) and \eqn{\beta} (i.e., the grand mean of the log odds which is a
#'   nuisance parameter). Log marginal likelihoods for \code{"H0"} and
#'   \code{"H1"} are obtained via Laplace approximations (see Kass &
#'   Vaidyanathan, 1992) which work well even for very small sample sizes. For
#'   the one-sided hypotheses \code{"H+"} and \code{"H-"} the log marginal
#'   likelihoods are obtained based on importance sampling which uses as a
#'   proposal a multivariate t distribution with location and scale matrix
#'   obtained via a Laplace approximation to the (log-transformed) posterior. If
#'   \code{posterior = TRUE}, posterior samples are obtained using importance
#'   sampling.
#' @return returns an object of class \code{"ab"} with components: \itemize{
#'   \item \code{input}: a list with the input arguments. \item \code{post}: a
#'   list with parameter posterior samples for the three hypotheses \code{"H1"},
#'   \code{"H+"} (in the output called \code{"Hplus"}), and \code{"H-"} (in the
#'   output called \code{"Hminus"}). Only contains samples if \code{posterior =
#'   TRUE}. \item \code{laplace}: a list with the approximate parameter
#'   posterior mode and variance/covariance matrix for each hypothesis obtained
#'   via a Laplace approximation. \item \code{method}: character that indicates
#'   the method that has been used to obtain the results. The default is
#'   \code{"log-is"} (importance sampling with multivariate t proposal based on
#'   a Laplace approximation to the log transformed posterior). If this method
#'   fails (for the one-sided hypotheses), method \code{"is-sn"} is used (i.e.,
#'   importance sampling is used to obtain unconstrained samples, then a
#'   skew-normal distribution is fitted to the samples to obtain the results for
#'   the one-sided hypotheses). If \code{method = "is-sn"}, posterior samples
#'   can only be obtained for \code{"H1"}. \item \code{logml}: a list with the
#'   estimated log marginal likelihoods for the hypotheses \code{"H0"} (i.e.,
#'   \code{"logml0"}), \code{"H1"} (i.e., \code{"logml1"}), \code{"H+"} (i.e.,
#'   \code{"logmlplus"}), and \code{"H-"} (i.e., \code{"logmlminus"}). \item
#'   \code{post_prob}: a named vector with the posterior probabilities of the
#'   four hypotheses \code{"H1"}, \code{"H+"}, \code{"H-"}, and \code{"H0"}.
#'   \item \code{logbf}: a list with the log Bayes factor in favor of
#'   \code{"H1"} over \code{"H0"}, the log Bayes factor in favor of \code{"H+"}
#'   over \code{"H0"}, and the log Bayes factor in favor of \code{"H-"} over
#'   \code{"H0"}. \item \code{bf}: a list with the Bayes factor in favor of
#'   \code{"H1"} over \code{"H0"} (i.e., \code{"bf10"}), the Bayes factor in
#'   favor of \code{"H+"} over \code{"H0"} (i.e., \code{"bfplus0"}), and the
#'   Bayes factor in favor of \code{"H-"} over \code{"H0"} (i.e.,
#'   \code{"bfminus0"}).}
#' @author Quentin F. Gronau
#' @references Kass, R. E., & Vaidyanathan, S. K. (1992). Approximate Bayes
#'   factors and orthogonal parameters, with application to testing equality of
#'   two binomial proportions. \emph{Journal of the Royal Statistical Society,
#'   Series B, 54}, 129-144. \url{https://doi.org/10.1111/j.2517-6161.1992.tb01868.x}
#' @example examples/example.ab_test.R
#'
#' @seealso \code{\link{elicit_prior}} allows the user to elicit a prior based
#'   on providing quantiles for either the log odds ratio, the odds ratio, the
#'   relative risk, or the absolute risk. The resulting prior is always
#'   translated to the corresponding normal prior on the log odds ratio. The
#'   \code{\link{plot_prior}} function allows the user to visualize the prior
#'   distribution. The \code{\link{simulate_priors}} function produces samples
#'   from the prior distribution. The prior and posterior probabilities of the
#'   hypotheses can be visualized using the \code{\link{prob_wheel}} function.
#'   Parameter posteriors can be visualized using the
#'   \code{\link{plot_posterior}} function. The \code{\link{plot_sequential}}
#'   function allows the user to sequentially plot the posterior probabilities
#'   of the hypotheses (only possible if the \code{data} object contains vectors
#'   with the cumulative "successes"/trials).
#'
#' @importFrom stats dnorm median nlminb pnorm qlogis
#' @importFrom VGAM log1pexp
#' @importFrom Matrix nearPD
#' @export
ab_test <- function(data,
                    prior_par = list(mu_psi = 0, sigma_psi = 1,
                                     mu_beta = 0, sigma_beta = 1),
                    prior_prob = NULL,
                    nsamples = 1e4,
                    is_df = 5,
                    posterior = FALSE) {


  ### argument checking ###

  # check data
  if ( ! is.list(data) ||
       ! all(c("y1", "y2", "n1", "n2") %in% names(data)) ||
       ! is.numeric(unlist(data))) {
    stop('data needs to be a named list with numeric elements
         "y1", "y2", "n1", and "n2',
         call. = FALSE)
  }

  if ( ! (all(vapply(data, length, 0) -
              mean(vapply(data, length, 0)) == 0))) {
    # elements are not of same length
    stop('data needs to be a named list with numeric elements
         "y1", "y2", "n1", and "n2 of same length',
         call. = FALSE)
  }

  # if data elements are vectors, only take last elements for ab test
  # the other ones will be used only for the sequential analysis plot
  data_raw <- data
  data <- list(y1 = data$y1[length(data$y1)],
               n1 = data$n1[length(data$n1)],
               y2 = data$y2[length(data$y2)],
               n2 = data$n2[length(data$n2)])

  # continue checks
  if (data$y1 < 0 || data$y2 < 0 || data$y1 > data$n1 ||
      data$y2 > data$n2) {
    stop('y1 needs to be a number between 0 and n1 and y2 needs to be
         a number between 0 and n2',
         call. = FALSE)
  }

  if (data$n1 <= 0 || data$n2 <= 0) {
    stop('n1 and n2 need to be larger than 0',
         call. = FALSE)
  }

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

  # check prior_prob
  if (is.null(prior_prob)) {
    prior_prob <- c(0, 1 / 4, 1 / 4, 1/ 2)
    names(prior_prob) <- c("H1", "H+", "H-", "H0")
  } else {

    if ( ! is.numeric(prior_prob) || length(prior_prob) != 4 ||
        ! all(c("H1", "H+", "H-", "H0") %in% names(prior_prob))) {
      stop('prior_prob needs to be a named numeric vector of length four
           with elements for "H1", "H+", "H-", "H0"',
           call. = FALSE)
    }

    # make sure that prior_prob has correct order
    prior_prob <- prior_prob[c("H1", "H+", "H-", "H0")]
  }

  if (sum(prior_prob) != 1) {
    stop('prior_prob needs to sum to one',
         call. = FALSE)
  }

  if (sum(prior_prob != 0) == 1) {
    stop('prior_prob needs to be non-zero for at least two hypotheses',
         call. = FALSE)
  }

  # check nsamples
  if (length(nsamples) > 1 || !is.numeric(nsamples) || nsamples <= 0) {
    stop('nsamples needs to be a single number larger than 0',
         call. = FALSE)
  }

  # check is_df
  if (length(is_df) > 1 || !is.numeric(is_df) || is_df <= 1) {
    stop('is_df needs to be a single number larger than or equal to 1',
         call. = FALSE)
  }

  # check posterior
  if ( ! is.logical(posterior)) {
    stop('posterior needs to be a Boolean',
         call. = FALSE)
  }


  # extract prior parameters
  mu_psi <- prior_par[["mu_psi"]]
  sigma_psi <- prior_par[["sigma_psi"]]
  mu_beta <- prior_par[["mu_beta"]]
  sigma_beta <- prior_par[["sigma_beta"]]

  ### H1 and H0 ###

  # starting values
  start_beta <- qlogis((data$y1 + data$y2) / (data$n1 + data$n2))
  eta1 <- qlogis(data$y1 / data$n1)
  eta2 <- qlogis(data$y2 / data$n2)

  if (is.infinite(start_beta)) {
    start_beta <- 0 # a finite number
  }
  if (is.infinite(eta1)) {
    eta1 <- 0 # a finite number
  }
  if (is.infinite(eta2)) {
    eta2 <- 0 # a finite number
  }

  start_psi <- eta2 - eta1

  # optimize
  o0 <- nlminb(start = start_beta, objective = minl0, gradient = gradient_minl0,
               data = data, mu_beta = mu_beta, sigma_beta = sigma_beta)
  o1 <- nlminb(start = c(start_beta, start_psi), objective = minl,
               gradient = gradient_minl, data = data,
               mu_beta = mu_beta, sigma_beta = sigma_beta,
               mu_psi = mu_psi, sigma_psi = sigma_psi)

  # approximate posterior modes
  mode_H0 <- o0$par
  names(mode_H0) <- "beta"
  mode_H1 <- o1$par
  parnames <-  c("beta", "psi")
  names(mode_H1) <- parnames

  # approximate posterior variance/covariance
  var_H0 <- inverse_hessian_minl0(par = o0$par, data = data,
                                  sigma_beta = sigma_beta)
  names(var_H0) <- "beta"
  cov_H1 <- inverse_hessian_minl(par = o1$par, data = data,
                                 sigma_beta = sigma_beta,
                                 sigma_psi = sigma_psi)
  cov_H1 <- as.matrix(nearPD(cov_H1)$mat) # to avoid numerical issues
  colnames(cov_H1) <- parnames
  rownames(cov_H1) <- parnames

  ### compute BF10 using Laplace approximation ###

  logml0 <- unname(.5 * log(2 * pi) +  .5 * log(var_H0) -
                     minl0(par = mode_H0, data = data, mu_beta = mu_beta,
                           sigma_beta = sigma_beta))
  logml1 <- unname(log(2 * pi) -
                     .5 * log(det_hessian_minl(par = mode_H1, data = data,
                                               sigma_beta = sigma_beta,
                                               sigma_psi = sigma_psi)) -
                     minl(par = mode_H1, data = data, mu_beta = mu_beta,
                          sigma_beta = sigma_beta, mu_psi = mu_psi,
                          sigma_psi = sigma_psi))
  logbf10 <- logml1 - logml0

  ### compute one-sided BFs ###

  method <- "log-is" # importance sampling with t-proposal based on
                     # Laplace approximation to log transformed posterior

  r <- try({

    if (start_psi > 0) {
      start_xi_plus <- log(start_psi)
      start_xi_minus <- log(.1)
    } else if (start_psi < 0) {
      start_xi_plus <- log(.1)
      start_xi_minus <- log( - start_psi)
    } else {
      start_xi_plus <- 0
      start_xi_minus <- 0
    }

    # H+
    oplus <- nlminb(start = c(start_beta, start_xi_plus), objective = minlplus,
                    gradient = gradient_minlplus, data = data,
                    mu_beta = mu_beta, sigma_beta = sigma_beta,
                    mu_psi = mu_psi, sigma_psi = sigma_psi)
    mode_Hplus <- oplus$par
    cov_Hplus <- inverse_hessian_minlplus(par = mode_Hplus,
                                          data = data,
                                          sigma_beta = sigma_beta,
                                          mu_psi = mu_psi,
                                          sigma_psi = sigma_psi)
    cov_Hplus <- as.matrix(nearPD(cov_Hplus)$mat) # to avoid numerical issues

    # proposal samples
    prop_samples_Hplus <- mvtnorm::rmvt(nsamples, delta = mode_Hplus,
                                        sigma = cov_Hplus, df = is_df)

    # importance weights
    weights_Hplus <- - apply_minlplus_cpp(prop_samples_Hplus, data$y1,
                                          data$y2, data$n1, data$n2,
                                          mu_beta, sigma_beta, mu_psi,
                                          sigma_psi) -
      mvtnorm::dmvt(prop_samples_Hplus, delta = mode_Hplus,
                    sigma = cov_Hplus, df = is_df, log = TRUE)
    logconst_Hplus <- median(weights_Hplus)
    logmlplus <- log(mean(exp(weights_Hplus - logconst_Hplus))) +
      logconst_Hplus

    # compute BF+0
    logbfplus0 <- logmlplus - logml0

    # H-
    ominus <- nlminb(start = c(start_beta, start_xi_minus),
                     objective = minlminus, gradient = gradient_minlminus,
                     data = data, mu_beta = mu_beta, sigma_beta = sigma_beta,
                     mu_psi = mu_psi, sigma_psi = sigma_psi)
    mode_Hminus <- ominus$par
    cov_Hminus <- inverse_hessian_minlminus(par = mode_Hminus,
                                            data = data,
                                            sigma_beta = sigma_beta,
                                            mu_psi = mu_psi,
                                            sigma_psi = sigma_psi)
    cov_Hminus <- as.matrix(nearPD(cov_Hminus)$mat) # to avoid numerical issues

    # proposal samples
    prop_samples_Hminus <- mvtnorm::rmvt(nsamples, delta = mode_Hminus,
                                         sigma = cov_Hminus, df = is_df)

    # importance weights
    weights_Hminus <- - apply_minlminus_cpp(prop_samples_Hminus,
                                            data$y1, data$y2,
                                            data$n1, data$n2, mu_beta,
                                            sigma_beta, mu_psi,
                                            sigma_psi) -
      mvtnorm::dmvt(prop_samples_Hminus, delta = mode_Hminus,
                    sigma = cov_Hminus, df = is_df, log = TRUE)
    logconst_Hminus <- median(weights_Hminus)
    logmlminus <- log(mean(exp(weights_Hminus - logconst_Hminus))) +
      logconst_Hminus

    # compute BF-0
    logbfminus0 <- logmlminus - logml0

  }, silent = TRUE)

  if (inherits(r, "try-error") || any(is.na(c(logbfplus0, logbfminus0)))) {

    method <- "is-sn" # importance sampling to obtain unconstrained samples
                      # then fit skew-normal and compute areas of interest

    mode_Hplus <- cov_Hplus <- NULL
    mode_Hminus <- cov_Hminus <- NULL

    # compute log of prior area larger than zero for psi
    logPriorPlus <- pnorm(0, mean = mu_psi, sd = sigma_psi,
                          lower.tail = FALSE, log.p = TRUE)
    # compute log of prior area smaller than zero for psi
    logPriorMinus <- pnorm(0, mean = mu_psi, sd = sigma_psi, log.p = TRUE)

    # obtain posterior samples via importance sampling

    # proposal samples
    prop_samples_H1 <- mvtnorm::rmvt(nsamples, delta = mode_H1,
                                     sigma = cov_H1, df = is_df)

    # importance weights
    weights_H1 <- - apply_minl_cpp(prop_samples_H1, data$y1, data$y2,
                                   data$n1, data$n2, mu_beta,
                                   sigma_beta, mu_psi, sigma_psi) -
      mvtnorm::dmvt(prop_samples_H1, delta = mode_H1, sigma = cov_H1,
                    df = is_df, log = TRUE)

    logconst_H1 <- median(weights_H1)

    # normalized importance weights
    normalized_weights_H1 <- exp(weights_H1 - logconst_H1) /
      sum(exp(weights_H1 - logconst_H1))

    # resample according to normalize weights to get posterior samples
    index <- sample(1:nrow(prop_samples_H1), size = nsamples,
                    replace = TRUE, prob = normalized_weights_H1)
    post_samples_H1 <- prop_samples_H1[index,]
    colnames(post_samples_H1) <- parnames

    # fit skew-normal distribution and compute area smaller/larger zero
    psi <- post_samples_H1[,"psi"]
    fit <- sn::selm(formula = psi ~ 1, family = "SN",
                    data = data.frame(psi = psi))
    logPostMinus <- log(sn::psn(x = 0, xi = fit@param$dp[["xi"]],
                                omega = fit@param$dp[["omega"]],
                                alpha = fit@param$dp[["alpha"]]))
    logPostPlus <- log(1 - exp(logPostMinus))

    logbfplus0 <- logbf10 + logPostPlus - logPriorPlus
    logbfminus0 <- logbf10 + logPostMinus - logPriorMinus

  }


  # compute posterior probabilities of hypotheses
  logbfs <- c(logbf10, logbfplus0, logbfminus0, 0)
  maxlogbf <- max(logbfs)
  post_prob <- exp(logbfs - maxlogbf) * prior_prob /
    sum(exp(logbfs - maxlogbf) * prior_prob)
  names(post_prob) <- names(prior_prob)


  post_samples_Hplus <- NULL
  post_samples_Hminus <- NULL

  if (posterior) {

    if (method == "is-sn") {

      # posterior samples for H1 have already been generated
      # generated and generating for Hplus and Hminus is not possible
      # hence, only generate implied posteriors

      beta <- post_samples_H1[,"beta"]
      psi <- post_samples_H1[,"psi"]
      p1 <- inv_logit(beta - psi / 2)
      p2 <- inv_logit(beta + psi / 2)
      or <- p2 / (1 - p2) / (p1 / (1 - p1))
      logor <- log(or)
      rrisk <- p2 / p1
      arisk <- p2 - p1

      post_samples_H1 <- data.frame(beta = beta, psi = psi, p1 = p1,
                                    p2 = p2, logor = logor, or = or,
                                    rrisk = rrisk, arisk = arisk)

    } else if (method == "log-is") {

      ### H1 ###

      # proposal samples
      prop_samples_H1 <- mvtnorm::rmvt(nsamples, delta = mode_H1,
                                       sigma = cov_H1, df = is_df)

      # importance weights
      weights_H1 <- - apply_minl_cpp(prop_samples_H1, data$y1, data$y2,
                                     data$n1, data$n2, mu_beta,
                                     sigma_beta, mu_psi, sigma_psi) -
        mvtnorm::dmvt(prop_samples_H1, delta = mode_H1, sigma = cov_H1,
                      df = is_df, log = TRUE)

      logconst_H1 <- median(weights_H1)

      # normalized importance weights
      normalized_weights_H1 <- exp(weights_H1 - logconst_H1) /
        sum(exp(weights_H1 - logconst_H1))

      # resample according to normalize weights to get posterior samples
      index <- sample(1:nrow(prop_samples_H1), size = nsamples,
                      replace = TRUE, prob = normalized_weights_H1)
      post_samples_H1 <- prop_samples_H1[index,]
      colnames(post_samples_H1) <- parnames
      beta <- post_samples_H1[,"beta"]
      psi <- post_samples_H1[,"psi"]
      p1 <- inv_logit(beta - psi / 2)
      p2 <- inv_logit(beta + psi / 2)
      or <- p2 / (1 - p2) / (p1 / (1 - p1))
      logor <- log(or)
      rrisk <- p2 / p1
      arisk <- p2 - p1

      post_samples_H1 <- data.frame(beta = beta, psi = psi, p1 = p1,
                                    p2 = p2, logor = logor, or = or,
                                    rrisk = rrisk, arisk = arisk)

      ### Hplus ###

      # normalized importance weights
      normalized_weights_Hplus <- exp(weights_Hplus - logconst_Hplus) /
        sum(exp(weights_Hplus - logconst_Hplus))

      # resample according to normalize weights to get posterior samples
      index <- sample(1:nrow(prop_samples_Hplus), size = nsamples,
                      replace = TRUE, prob = normalized_weights_Hplus)
      post_samples_Hplus <- prop_samples_Hplus[index,]
      colnames(post_samples_Hplus) <- c("beta", "log(psi)")
      beta <- post_samples_Hplus[,"beta"]
      psi <- exp(post_samples_Hplus[,"log(psi)"])
      p1 <- inv_logit(beta - psi / 2)
      p2 <- inv_logit(beta + psi / 2)
      or <- p2 / (1 - p2) / (p1 / (1 - p1))
      logor <- log(or)
      rrisk <- p2 / p1
      arisk <- p2 - p1

      post_samples_Hplus <- data.frame(beta = beta, psi = psi, p1 = p1,
                                       p2 = p2, logor = logor, or = or,
                                       rrisk = rrisk, arisk = arisk)

      ### Hminus ###

      # normalized importance weights
      normalized_weights_Hminus <- exp(weights_Hminus - logconst_Hminus) /
        sum(exp(weights_Hminus - logconst_Hminus))

      # resample according to normalize weights to get posterior samples
      index <- sample(1:nrow(prop_samples_Hminus), size = nsamples,
                      replace = TRUE, prob = normalized_weights_Hminus)
      post_samples_Hminus <- prop_samples_Hminus[index,]
      colnames(post_samples_Hminus) <- c("beta", "log(-psi)")
      beta <- post_samples_Hminus[,"beta"]
      psi <- -exp(post_samples_Hminus[,"log(-psi)"])
      p1 <- inv_logit(beta - psi / 2)
      p2 <- inv_logit(beta + psi / 2)
      or <- p2 / (1 - p2) / (p1 / (1 - p1))
      logor <- log(or)
      rrisk <- p2 / p1
      arisk <- p2 - p1

      post_samples_Hminus <- data.frame(beta = beta, psi = psi, p1 = p1,
                                        p2 = p2, logor = logor, or = or,
                                        rrisk = rrisk, arisk = arisk)

    }

  } else {
    post_samples_H1 <- NULL
  }

  out <- list(input = list(data = data_raw,
                           prior_par = prior_par,
                           prior_prob = prior_prob,
                           nsamples = nsamples,
                           is_df = is_df,
                           posterior = posterior),
              post = list(H1 = post_samples_H1,
                          Hplus = post_samples_Hplus,
                          Hminus = post_samples_Hminus),
              laplace = list(mode_H0 = mode_H0,
                             var_H0 = var_H0,
                             mode_H1 = mode_H1,
                             cov_H1 = cov_H1,
                             mode_Hplus = mode_Hplus,
                             cov_Hplus = cov_Hplus,
                             mode_Hminus = mode_Hminus,
                             cov_Hminus = cov_Hminus),
              method = method,
              logml = list(logml0 = logml0,
                           logml1 = logml1,
                           logmlplus = logmlplus,
                           logmlminus = logmlminus),
              post_prob = post_prob,
              logbf = list(bf10 = logbf10,
                           bfplus0 = logbfplus0,
                           bfminus0 = logbfminus0),
              bf = list(bf10 = exp(logbf10),
                        bfplus0 = exp(logbfplus0),
                        bfminus0 = exp(logbfminus0)))

  class(out) <- "ab"

  return(out)

}

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
