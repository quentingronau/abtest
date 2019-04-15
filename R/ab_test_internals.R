#--------------------------------------------------------------------------
# Functions for Laplace Approximation (Kass & Vaidyanathan, 1992)
#--------------------------------------------------------------------------

# minus l0(beta)
minl0 <- function(par, data, mu_beta = 0, sigma_beta = 1) {

  # parameter
  beta <- par[1]

  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2

  # log probabilities
  logp <- beta - log1pexp(beta)
  log1minp <- - log1pexp(beta)

  out <- - (y1 + y2) * logp - (n1 + n2 - y1 - y2) * log1minp -
    dnorm(beta, mu_beta, sigma_beta, log = TRUE)

  return(out)

}

# gradient of minus l0(beta)
gradient_minl0 <- function(par, data, mu_beta = 0,
                           sigma_beta = 1) {

  # parameter
  beta <- par[1]

  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2

  out <- - (y1 + y2 - (n1 + n2 - y1 - y2) * exp(beta)) /
    (1 + exp(beta)) + (beta - mu_beta) / sigma_beta^2

  return(out)

}


# Hessian for minus l0(beta)
hessian_minl0 <- function(par, data, sigma_beta = 1) {

  # parameter
  beta <- par[1]

  # data
  n1 <- data$n1
  n2 <- data$n2

  hessian <- (n1 + n2) * exp(beta) / (1 + exp(beta))^2 + 1 / sigma_beta^2
  return(hessian)

}

# inverse of Hessian for minus l0(beta)
inverse_hessian_minl0 <- function(par, data, sigma_beta = 1) {

  hessian <- hessian_minl0(par = par, data = data, sigma_beta = sigma_beta)
  return(1 / hessian)

}

# minus l(beta, psi)
minl <- function(par, data, mu_beta = 0, sigma_beta = 1,
                 mu_psi = 0, sigma_psi = 1) {

  # parameters
  beta <- par[1]
  psi <- par[2]

  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2

  # log probabilities
  logp1 <- beta - .5 * psi - log1pexp(beta - .5 * psi)
  log1minp1 <- - log1pexp(beta - .5 * psi)
  logp2 <- beta + .5 * psi - log1pexp(beta + .5 * psi)
  log1minp2 <- - log1pexp(beta + .5 * psi)

  out <- y1 * logp1 + (n1 - y1) * log1minp1 +
    y2 * logp2 + (n2 - y2) * log1minp2 +
    dnorm(beta, mu_beta, sigma_beta, log = TRUE) +
    dnorm(psi, mu_psi, sigma_psi, log = TRUE)

  return(- out)

}

# gradient of minus l(beta, psi)
gradient_minl <- function(par, data, mu_beta = 0,
                          sigma_beta = 1, mu_psi = 0,
                          sigma_psi = 1) {

  # parameters
  beta <- par[1]
  psi <- par[2]

  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2

  partial_beta <- (y1 - (n1 - y1) * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi)) + (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi)) - (beta - mu_beta) / sigma_beta^2

  partial_psi <- .5 * (((n1 - y1) * exp(beta - .5 * psi) - y1) /
                         (1 + exp(beta - .5 * psi)) +
                         (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
                         (1 + exp(beta + .5 * psi))) -
    (psi - mu_psi) / sigma_psi^2

  return(- c(partial_beta, partial_psi))

}

# Hessian for minus l(beta, psi)
hessian_minl <- function(par, data, sigma_beta = 1, sigma_psi = 1) {

  # parameters
  beta <- par[1]
  psi <- par[2]

  # data
  n1 <- data$n1
  n2 <- data$n2

  partial_beta2 <-  (n1 * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi))^2 + (n2 * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi))^2 + 1 / sigma_beta^2

  partial_psi2 <- .25 * ((n1 * exp(beta - .5 * psi)) /
                           (1 + exp(beta - .5 * psi))^2 +
                           (n2 * exp(beta + .5 * psi)) /
                           (1 + exp(beta + .5 * psi))^2) +
    1 / sigma_psi^2

  partial_beta_psi <- - .5 * ((n1 * exp(beta - .5 * psi)) /
                                (1 + exp(beta - .5 * psi))^2 -
                                (n2 * exp(beta + .5 * psi)) /
                                (1 + exp(beta + .5 * psi))^2)

  hessian <- matrix(c(partial_beta2, rep(partial_beta_psi, 2),
                      partial_psi2), 2, 2, byrow = TRUE)
  return(hessian)

}

# determinant of Hessian for minus l(beta, psi)
det_hessian_minl <- function(par, data, sigma_beta = 1, sigma_psi = 1) {

  # parameters
  beta <- par[1]
  psi <- par[2]

  # data
  n1 <- data$n1
  n2 <- data$n2

  # out <- (n1 * n2 * exp(2 * beta)) /
  #   ((1 + exp(beta - .5 * psi))^2 * (1 + exp(beta + .5 * psi))^2) +
  #   (1 / sigma_psi^2 + 1 / (4 * sigma_beta^2)) *
  #   ((n1 * exp(beta - .5 * psi)) / (1 + exp(beta - .5 * psi))^2 +
  #      (n2 * exp(beta + .5 * psi)) / (1 + exp(beta + .5 * psi))^2) +
  #   1 / (sigma_beta^2 * sigma_psi^2)

  hessian <- hessian_minl(par = par, data = data,
                          sigma_beta = sigma_beta,
                          sigma_psi = sigma_psi)
  det_hessian <- hessian[1,1] * hessian[2,2] - hessian[1,2]^2

  return(det_hessian)

}

# inverse of Hessian for minus l(beta, psi)
inverse_hessian_minl <- function(par, data, sigma_beta = 1,
                                 sigma_psi = 1) {

  # parameters
  beta <- par[1]
  psi <- par[2]

  # data
  n1 <- data$n1
  n2 <- data$n2

  det_hessian <- det_hessian_minl(par = par, data = data,
                                  sigma_beta = sigma_beta,
                                  sigma_psi = sigma_beta)

  partial_beta2 <-  (n1 * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi))^2 + (n2 * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi))^2 + 1 / sigma_beta^2

  partial_psi2 <- .25 * ((n1 * exp(beta - .5 * psi)) /
                           (1 + exp(beta - .5 * psi))^2 +
                           (n2 * exp(beta + .5 * psi)) /
                           (1 + exp(beta + .5 * psi))^2) +
    1 / sigma_psi^2

  minus_partial_beta_psi <- .5 * ((n1 * exp(beta - .5 * psi)) /
                                    (1 + exp(beta - .5 * psi))^2 -
                                    (n2 * exp(beta + .5 * psi)) /
                                    (1 + exp(beta + .5 * psi))^2)

  inv_hessian <- 1 / det_hessian *
    matrix(c(partial_psi2, rep(minus_partial_beta_psi, 2),
             partial_beta2), 2, 2, byrow = TRUE)

  return(inv_hessian)

}

# minus l+(beta, eta)
minlplus <- function(par, data, mu_beta = 0, sigma_beta = 1,
                     mu_psi = 0, sigma_psi = 1) {


  # parameters
  beta <- par[1]
  xi <- par[2]
  psi <- exp(xi)

  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2

  # log probabilities
  logp1 <- beta - .5 * psi - log1pexp(beta - .5 * psi)
  log1minp1 <- - log1pexp(beta - .5 * psi)
  logp2 <- beta + .5 * psi - log1pexp(beta + .5 * psi)
  log1minp2 <- - log1pexp(beta + .5 * psi)

  out <- y1 * logp1 + (n1 - y1) * log1minp1 +
    y2 * logp2 + (n2 - y2) * log1minp2 +
    dnorm(beta, mu_beta, sigma_beta, log = TRUE) +
    dnorm(psi, mu_psi, sigma_psi, log = TRUE) -
    pnorm(0, mu_psi, sigma_psi, lower.tail = FALSE, log.p = TRUE) + xi

  return(- out)

}

# gradient of minus l+(beta, xi)
gradient_minlplus <- function(par, data, mu_beta = 0,
                              sigma_beta = 1, mu_psi = 0,
                              sigma_psi = 1) {

  # parameters
  beta <- par[1]
  xi <- par[2]
  psi <- exp(xi)

  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2

  partial_beta <- (y1 - (n1 - y1) * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi)) + (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi)) - (beta - mu_beta) / sigma_beta^2

  partial_xi <- .5 * psi * (((n1 - y1) * exp(beta - .5 * psi) - y1) /
                              (1 + exp(beta - .5 * psi)) +
                              (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
                              (1 + exp(beta + .5 * psi))) -
    psi * (psi - mu_psi) / sigma_psi^2 + 1

  return(- c(partial_beta, partial_xi))

}

# Hessian for minus l+(beta, xi)
hessian_minlplus <- function(par, data, sigma_beta = 1,
                             mu_psi = 0, sigma_psi = 1) {

  # parameters
  beta <- par[1]
  xi <- par[2]
  psi <- exp(xi)

  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2

  partial_beta2 <-  (n1 * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi))^2 + (n2 * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi))^2 + 1 / sigma_beta^2

  partial_xi2 <- - .5 * psi * (((n1 - y1) * exp(beta - .5 * psi) - y1) /
                                 (1 + exp(beta - .5 * psi)) +
                                 (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
                                 (1 + exp(beta + .5 * psi)) -
                                 .5 * psi * n1 * exp(beta - .5 * psi) /
                                 (1 + exp(beta - .5 * psi))^2 -
                                 .5 * psi * n2 * exp(beta + .5 * psi) /
                                 (1 + exp(beta + .5 * psi))^2) +
    psi * (2 * psi - mu_psi) / sigma_psi^2

  partial_beta_xi <- - .5 * psi * ((n1 * exp(beta - .5 * psi)) /
                                     (1 + exp(beta - .5 * psi))^2 -
                                     (n2 * exp(beta + .5 * psi)) /
                                     (1 + exp(beta + .5 * psi))^2)

  hessian <- matrix(c(partial_beta2, rep(partial_beta_xi, 2),
                      partial_xi2), 2, 2, byrow = TRUE)

  return(hessian)

}

# determinant of Hessian for minus l+(beta, xi)
det_hessian_minlplus <- function(par, data, sigma_beta = 1, mu_psi = 0,
                                 sigma_psi = 1) {

  hessian <- hessian_minlplus(par = par, data = data,
                              sigma_beta = sigma_beta,
                              mu_psi = mu_psi,
                              sigma_psi = sigma_psi)
  det_hessian <- hessian[1,1] * hessian[2,2] - hessian[1,2]^2

  return(det_hessian)

}

# inverse of Hessian for minus l+(beta, xi)
inverse_hessian_minlplus <- function(par, data, sigma_beta = 1, mu_psi = 0,
                                     sigma_psi = 1) {


  det_hessian <- det_hessian_minlplus(par = par, data = data,
                                      sigma_beta = sigma_beta,
                                      mu_psi = mu_psi,
                                      sigma_psi = sigma_beta)

  hessian <- hessian_minlplus(par = par, data = data,
                              sigma_beta = sigma_beta,
                              mu_psi = mu_psi,
                              sigma_psi = sigma_beta)

  inv_hessian <- 1 / det_hessian *
    matrix(c(hessian[2,2], rep(- hessian[1,2], 2),
             hessian[1,1]), 2, 2, byrow = TRUE)

  return(inv_hessian)

}

# minus l-(beta, xi)
minlminus <- function(par, data, mu_beta = 0, sigma_beta = 1,
                      mu_psi = 0, sigma_psi = 1) {

  # parameters
  beta <- par[1]
  xi <- par[2]
  psi <- - exp(xi)

  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2

  # log probabilities
  logp1 <- beta - .5 * psi - log1pexp(beta - .5 * psi)
  log1minp1 <- - log1pexp(beta - .5 * psi)
  logp2 <- beta + .5 * psi - log1pexp(beta + .5 * psi)
  log1minp2 <- - log1pexp(beta + .5 * psi)

  out <- y1 * logp1 + (n1 - y1) * log1minp1 +
    y2 * logp2 + (n2 - y2) * log1minp2 +
    dnorm(beta, mu_beta, sigma_beta, log = TRUE) +
    dnorm(psi, mu_psi, sigma_psi, log = TRUE) -
    pnorm(0, mu_psi, sigma_psi, log.p = TRUE) + xi

  return(- out)

}

# gradient of minus l+(beta, xi)
gradient_minlminus <- function(par, data, mu_beta = 0,
                               sigma_beta = 1, mu_psi = 0,
                               sigma_psi = 1) {

  # parameters
  beta <- par[1]
  xi <- par[2]
  psi <- - exp(xi)

  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2

  partial_beta <- (y1 - (n1 - y1) * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi)) + (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi)) - (beta - mu_beta) / sigma_beta^2

  partial_xi <- - .5 * psi * ((y1 - (n1 - y1) * exp(beta - .5 * psi)) /
                                (1 + exp(beta - .5 * psi)) +
                                ((n2 - y2) * exp(beta + .5 * psi) - y2) /
                                (1 + exp(beta + .5 * psi))) -
    psi * (psi - mu_psi) / sigma_psi^2 + 1

  return(- c(partial_beta, partial_xi))

}

# Hessian for minus l-(beta, xi)
hessian_minlminus <- function(par, data, sigma_beta = 1,
                              mu_psi = 0, sigma_psi = 1) {

  # parameters
  beta <- par[1]
  xi <- par[2]
  psi <- - exp(xi)

  # data
  y1 <- data$y1
  y2 <- data$y2
  n1 <- data$n1
  n2 <- data$n2

  partial_beta2 <-  (n1 * exp(beta - .5 * psi)) /
    (1 + exp(beta - .5 * psi))^2 + (n2 * exp(beta + .5 * psi)) /
    (1 + exp(beta + .5 * psi))^2 + 1 / sigma_beta^2

  partial_xi2 <- - .5 * psi * (((n1 - y1) * exp(beta - .5 * psi) - y1) /
                                 (1 + exp(beta - .5 * psi)) +
                                 (y2 - (n2 - y2) * exp(beta + .5 * psi)) /
                                 (1 + exp(beta + .5 * psi)) -
                                 .5 * psi * n1 * exp(beta - .5 * psi) /
                                 (1 + exp(beta - .5 * psi))^2 -
                                 .5 * psi * n2 * exp(beta + .5 * psi) /
                                 (1 + exp(beta + .5 * psi))^2) +
    psi * (2 * psi - mu_psi) / sigma_psi^2

  partial_beta_xi <- - .5 * psi * ((n1 * exp(beta - .5 * psi)) /
                                     (1 + exp(beta - .5 * psi))^2 -
                                     (n2 * exp(beta + .5 * psi)) /
                                     (1 + exp(beta + .5 * psi))^2)

  hessian <- matrix(c(partial_beta2, rep(partial_beta_xi, 2),
                      partial_xi2), 2, 2, byrow = TRUE)

  return(hessian)

}

# determinant of Hessian for minus l-(beta, xi)
det_hessian_minlminus <- function(par, data, sigma_beta = 1, mu_psi = 0,
                                  sigma_psi = 1) {

  hessian <- hessian_minlminus(par = par, data = data,
                               sigma_beta = sigma_beta,
                               mu_psi = mu_psi,
                               sigma_psi = sigma_psi)

  det_hessian <- hessian[1,1] * hessian[2,2] - hessian[1,2]^2

  return(det_hessian)

}

# inverse of Hessian for minus l-(beta, xi)
inverse_hessian_minlminus <- function(par, data, sigma_beta = 1,
                                      mu_psi = 0, sigma_psi = 1) {

  det_hessian <- det_hessian_minlminus(par = par, data = data,
                                       sigma_beta = sigma_beta,
                                       mu_psi = mu_psi,
                                       sigma_psi = sigma_beta)

  hessian <- hessian_minlminus(par = par, data = data,
                               sigma_beta = sigma_beta,
                               mu_psi = mu_psi,
                               sigma_psi = sigma_beta)

  inv_hessian <- 1 / det_hessian *
    matrix(c(hessian[2,2], rep(- hessian[1,2], 2),
             hessian[1,1]), 2, 2, byrow = TRUE)

  return(inv_hessian)

}
