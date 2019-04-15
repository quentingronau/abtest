
#--------------------------------------------------------------------------
# logit and inverse logit functions
#--------------------------------------------------------------------------

logit <- function(x) {
  log(x) - log1p(- x)
}

inv_logit <- function(x) {
  1/(1 + exp(-x))
}

#--------------------------------------------------------------------------
# log odds ratio
#--------------------------------------------------------------------------

# psi *is* the log odds ratio!

#' @importFrom truncnorm dtruncnorm ptruncnorm qtruncnorm
#' @importFrom stats dlnorm integrate plnorm
dlogor <- function(x, mu_psi, sigma_psi, hypothesis = "H1") {

  bounds <- switch(hypothesis,
                   "H1" =  c(-Inf, Inf),
                   "H+" =  c(0, Inf),
                   "H-" =  c(-Inf, 0))

  return(dtruncnorm(x, a = bounds[1], b = bounds[2],
                    mean = mu_psi, sd = sigma_psi))

}

plogor <- function(q, mu_psi, sigma_psi, hypothesis = "H1") {

  bounds <- switch(hypothesis,
                   "H1" =  c(-Inf, Inf),
                   "H+" =  c(0, Inf),
                   "H-" =  c(-Inf, 0))

  return(ptruncnorm(q, a = bounds[1], b = bounds[2],
                    mean = mu_psi, sd = sigma_psi))

}


#--------------------------------------------------------------------------
# odds ratio
#--------------------------------------------------------------------------

# distribution

dor <- function(x, mu_psi, sigma_psi, hypothesis = "H1") {

  bounds <- switch(hypothesis,
                   "H1" =  c(0, Inf),
                   "H+" =  c(exp(0), Inf),
                   "H-" =  c(0, exp(0)))

  return(ifelse(x < bounds[1] | x > bounds[2], 0,
                dlnorm(x, meanlog = mu_psi, sdlog = sigma_psi) /
                (plnorm(bounds[2], meanlog = mu_psi, sdlog = sigma_psi) -
                plnorm(bounds[1], meanlog = mu_psi, sdlog = sigma_psi))))



}

por <- function(q, mu_psi, sigma_psi, hypothesis = "H1") {

  bounds <- switch(hypothesis,
                   "H1" =  c(0, Inf),
                   "H+" =  c(exp(0), Inf),
                   "H-" =  c(0, exp(0)))

  num <- plnorm(pmax(pmin(q, bounds[2]), bounds[1]),
         meanlog = mu_psi, sdlog = sigma_psi) -
    plnorm(bounds[1], meanlog = mu_psi, sdlog = sigma_psi)
  denom <- plnorm(bounds[2], meanlog = mu_psi, sdlog = sigma_psi) -
    plnorm(bounds[1], meanlog = mu_psi, sdlog = sigma_psi)

  return(num / denom)

}

#--------------------------------------------------------------------------
# relative risk
#--------------------------------------------------------------------------

# pdf of relative risk
drrisk <- function(x, mu_psi, sigma_psi, mu_beta, sigma_beta,
                   hypothesis = "H1") {

  bounds <- switch(hypothesis,
                   "H1" =  c(-Inf, Inf),
                   "H+" =  c(0, Inf),
                   "H-" =  c(-Inf, 0))

  out <- vapply(x, FUN = function(y, mu_beta, sigma_beta,
                                       mu_psi, sigma_psi, bounds) {

    suppressWarnings({

      integrate(function(beta, mu_beta, sigma_beta,
                         mu_psi, sigma_psi, y, bounds) {

        out <- dnorm(beta, mu_beta, sigma_beta) *
          dtruncnorm(2 * log((-(1 - y) * exp(beta) +
                                sqrt((1 - y)^2 * exp(2 * beta) +
                                       4 * y)) / 2),
                     a = bounds[1], b = bounds[2],
                     mean = mu_psi, sd = sigma_psi) *
          2 * (exp(beta) + (2 - (1 - y) * exp(2 * beta)) /
                 sqrt((1 - y)^2 * exp(2 * beta) + 4 * y)) /
          (-(1 - y) * exp(beta) +
             sqrt((1 - y)^2 * exp(2 * beta) + 4 * y))

        out[is.na(out)] <- 0

        return(out)

      }, -Inf,  Inf,
      mu_beta = mu_beta, sigma_beta = sigma_beta,
      mu_psi = mu_psi, sigma_psi = sigma_psi,
      y = y, bounds = bounds)$value

    })

  }, FUN.VALUE = 0, mu_beta = mu_beta, sigma_beta = sigma_beta,
  mu_psi = mu_psi, sigma_psi = sigma_psi, bounds = bounds)

  return(out)

}

# cdf of relative risk
prrisk <- function(q, mu_psi, sigma_psi, mu_beta, sigma_beta,
                   hypothesis = "H1") {

  bounds <- switch(hypothesis,
                   "H1" =  c(-Inf, Inf),
                   "H+" =  c(0, Inf),
                   "H-" =  c(-Inf, 0))

  out <- vapply(q, FUN = function(x, mu_beta, sigma_beta,
                                  mu_psi, sigma_psi, bounds) {

    suppressWarnings({

      integrate(function(beta, mu_beta, sigma_beta,
                         mu_psi, sigma_psi, x, bounds) {

        out <- dnorm(beta, mu_beta, sigma_beta) *
          ptruncnorm(2 * log((-(1 - x) * exp(beta) +
                           sqrt((1 - x)^2 * exp(2 * beta) + 4 * x)) / 2),
                     a = bounds[1], b = bounds[2],
                     mean = mu_psi, sd = sigma_psi)

        out[is.na(out)] <- 0

        return(out)

      }, -Inf,  Inf,
      mu_beta = mu_beta, sigma_beta = sigma_beta,
      mu_psi = mu_psi, sigma_psi = sigma_psi,
      x = x, bounds = bounds)$value

    })

  }, FUN.VALUE = 0, mu_beta = mu_beta, sigma_beta = sigma_beta,
  mu_psi = mu_psi, sigma_psi = sigma_psi, bounds = bounds)

  return(out)

}

#--------------------------------------------------------------------------
# absolute risk
#--------------------------------------------------------------------------

# pdf of absolute risk
darisk <- function(x, mu_psi, sigma_psi,  mu_beta, sigma_beta,
                   hypothesis = "H1") {

  if (hypothesis == "H+") {
    x[x == 0] <- 1e-5
  } else if (hypothesis == "H-") {
    x[x == 0] <- -1e-5
  }

  bounds <- switch(hypothesis,
                   "H1" =  c(-Inf, Inf),
                   "H+" =  c(0, Inf),
                   "H-" =  c(-Inf, 0))

  out <- vapply(x, function(y, mu_beta, sigma_beta,
                           mu_psi, sigma_psi, bounds) {

    suppressWarnings({

      integrate(function(beta, mu_beta, sigma_beta, mu_psi,
                         sigma_psi, x, y, bounds) {

        out <- dnorm(beta, mu_beta, sigma_beta) *
          dtruncnorm(2 * log((y * (exp(2 * beta) + 1) +
                                sqrt(y^2 * (exp(2 * beta) - 1)^2 +
                                       4 * exp(2 * beta))) /
                               (2 * exp(beta) * (1 - y))),
                     a = bounds[1], b = bounds[2],
                     mean = mu_psi, sd = sigma_psi) *
          2 * ( (exp(2 * beta) + (y * (exp(2 * beta) - 1)^2) /
                   sqrt(y^2 * (exp(2 * beta) - 1)^2 +
                          4 * exp(2 * beta)) + 1 ) /
                  (y * (exp(2 * beta) + 1) +
                     sqrt(y^2 * (exp(2 * beta) - 1)^2 +
                            4 * exp(2 * beta))) + 1 / (1 - y) )

        out[is.na(out)] <- 0

        return(out)

      }, -Inf,  Inf,
      mu_beta = mu_beta, sigma_beta = sigma_beta,
      mu_psi = mu_psi, sigma_psi = sigma_psi,
      x = x, y = y, bounds = bounds)$value

    })

  }, FUN.VALUE = 0, mu_beta = mu_beta, sigma_beta = sigma_beta,
  mu_psi = mu_psi, sigma_psi = sigma_psi, bounds = bounds)

  return(out)

}

# cdf of absolute risk
parisk <- function(q, mu_psi, sigma_psi, mu_beta, sigma_beta,
                   hypothesis = "H1") {

  bounds <- switch(hypothesis,
                   "H1" =  c(-Inf, Inf),
                   "H+" =  c(0, Inf),
                   "H-" =  c(-Inf, 0))

  out <- vapply(q, function(upsilon, mu_beta, sigma_beta,
                            mu_psi, sigma_psi, bounds) {

    suppressWarnings({

      integrate(function(beta, mu_beta, sigma_beta, mu_psi,
                         sigma_psi, upsilon, bounds) {
        out <- dnorm(beta, mu_beta, sigma_beta) *
          ptruncnorm(2 * log((upsilon * (exp(2 * beta) + 1) +
                                sqrt(upsilon^2 * (exp(2 * beta) - 1)^2 +
                                       4 * exp(2 * beta))) /
                               (2 * exp(beta) * (1 - upsilon))),
                     a = bounds[1], b = bounds[2],
                     mean = mu_psi, sd = sigma_psi)

        out[is.na(out)] <- 0

        return(out)

      }, -Inf,  Inf,
      mu_beta = mu_beta, sigma_beta = sigma_beta,
      mu_psi = mu_psi, sigma_psi = sigma_psi,
      upsilon = upsilon, bounds = bounds)$value

    })

  }, FUN.VALUE = 0, mu_beta = mu_beta, sigma_beta = sigma_beta,
  mu_psi = mu_psi, sigma_psi = sigma_psi, bounds = bounds)

  return(out)

}

#--------------------------------------------------------------------------
# joint distribution p1 and p2
#--------------------------------------------------------------------------

# implied joint pdf for p1 and p2
dp1p2 <- function(p1, p2, mu_psi = 0, sigma_psi = 1,
                  mu_beta = 0, sigma_beta = 1,
                  hypothesis = "H1") {

  bounds <- switch(hypothesis,
                   "H1" =  c(-Inf, Inf),
                   "H+" =  c(0, Inf),
                   "H-" =  c(-Inf, 0))

  out <- mapply(FUN = function(x1, x2, mu_beta, sigma_beta,
                               mu_psi, sigma_psi, bounds) {

    if (x1 >= 1 || x1 <= 0 || x2 >= 1 || x2 <= 0) {
      out <- 0
    } else {
      out <-
        dnorm(.5 * (logit(x1) + logit(x2)), mean = mu_beta,
              sd = sigma_beta) *
        dtruncnorm(logit(x2) - logit(x1), a = bounds[1], b = bounds[2],
                   mean = mu_psi, sd = sigma_psi) *
        1/(x1 * x2 * (1 - x1) * (1 - x2)) # Jacobian
    }

  }, p1, p2, MoreArgs = list(mu_psi = mu_psi, sigma_psi = sigma_psi,
  mu_beta = mu_beta, sigma_beta = sigma_beta, bounds = bounds))

  return(out)

}

#--------------------------------------------------------------------------
# marginal distribution of p1
#--------------------------------------------------------------------------

dp1 <- function(x, mu_psi = 0, sigma_psi = 1,
                mu_beta = 0, sigma_beta = 1,
                hypothesis = "H1") {

  out <- vapply(x, function(y, mu_psi, sigma_psi,
                            mu_beta, sigma_beta,
                            hypothesis) {

    tryCatch({

      suppressWarnings({

        integrate(function(p2, p1, mu_psi, sigma_psi,
                           mu_beta, sigma_beta, hypothesis) {

          dp1p2(p1 = p1, p2 = p2, mu_psi = mu_psi, sigma_psi = sigma_psi,
                mu_beta = mu_beta, sigma_beta = sigma_beta,
                hypothesis = hypothesis)

        }, lower = 0, upper = 1, p1 = y, mu_psi = mu_psi,
        sigma_psi = sigma_psi, mu_beta = mu_beta,
        sigma_beta = sigma_beta,
        hypothesis = hypothesis)$value

      })
    }, error = function(e) 0)

  }, FUN.VALUE = 0, mu_psi = mu_psi, sigma_psi = sigma_psi,
  mu_beta = mu_beta, sigma_beta = sigma_beta,
  hypothesis = hypothesis)

  return(out)

}

#--------------------------------------------------------------------------
# marginal distribution of p2
#--------------------------------------------------------------------------

dp2 <- function(x, mu_psi = 0, sigma_psi = 1,
                mu_beta = 0, sigma_beta = 1,
                hypothesis = "H1") {

  out <- vapply(x, function(y, mu_psi, sigma_psi,
                            mu_beta, sigma_beta,
                            hypothesis) {

    tryCatch({

      suppressWarnings({

        integrate(dp1p2, lower = 0, upper = 1, p2 = y, mu_psi = mu_psi,
                  sigma_psi = sigma_psi, mu_beta = mu_beta,
                  sigma_beta = sigma_beta,
                  hypothesis = hypothesis)$value

      })
    }, error = function(e) 0)

  }, FUN.VALUE = 0, mu_psi = mu_psi, sigma_psi = sigma_psi,
  mu_beta = mu_beta, sigma_beta = sigma_beta,
  hypothesis = hypothesis)

  return(out)

}

#--------------------------------------------------------------------------
# conditional distribution of p2 given p1
#--------------------------------------------------------------------------

dp2givenp1 <- function(p2, p1 = 0.5, mu_psi = 0, sigma_psi = 1,
                       mu_beta = 0, sigma_beta = 1,
                       hypothesis = "H1") {

  bounds <- switch(hypothesis,
                   "H1" =  c(-Inf, Inf),
                   "H+" =  c(0, Inf),
                   "H-" =  c(-Inf, 0))

  num <- vapply(p2, function(x2, p1, mu_psi, sigma_psi,
                             mu_beta, sigma_beta, hypothesis) {

    dp1p2(p1, x2, mu_psi = mu_psi, sigma_psi = sigma_psi,
          mu_beta = mu_beta, sigma_beta = sigma_beta,
          hypothesis = hypothesis)

  }, FUN.VALUE = 0, p1 = p1, mu_psi = mu_psi, sigma_psi = sigma_psi,
  mu_beta = mu_beta, sigma_beta = sigma_beta,
  hypothesis = hypothesis)

  denom <- suppressWarnings({

    integrate(dp1p2, 0, 1, p1 = p1, mu_psi = mu_psi, sigma_psi = sigma_psi,
              mu_beta = mu_beta, sigma_beta = sigma_beta,
              hypothesis = hypothesis)$value

  })

  out <- num / denom

  return(out)

}
