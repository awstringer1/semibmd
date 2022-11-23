# Helper functions for semibmd

#' Reflection for bounded newton's method
#'
#' Implement the reflection of Coleman and Li (1994)
#'
#' @inheritParams bounded_newton
#' @param yt Point to reflect
#' @return Numeric value representing the reflected point
#'
#' @details Implements the formula \code{min(wt,2*(ub - lb) - wt) + lb} where
#' \code{wt = abs(yt - lb) mod (2*(ub-lb))}
#' and \code{lb = min(bounds),ub = max(bounds)}.
reflect <- function(yt,bounds) {
  lb <- min(bounds)
  ub <- max(bounds)
  wt <- abs(yt - lb) %% (2*(ub-lb))
  min(wt,2*(ub - lb) - wt) + lb
}

#' Get the benchmark dose and lower limit
#'
#' Helper function to extract the bmd and bmdl from an object of class \code{semibmd}
#' returned by \code{benchmark_dose}.
#'
#' @param object Object of class \code{semibmd}
#' returned by \code{benchmark_dose}
#'
#' @return Numeric vector of length 2 containing \code{c(bmd,bmdl)}.
#'
#' @export
get_bmd <- function(object) with(object,c(bmd,bmdl))

#' Get the dose-response model
#'
#' Helper function to extract the dose-response model object, of class \code{scam}
#' or \code{gam}, from an object of class \code{semibmd}
#' returned by \code{benchmark_dose}.
#'
#' @param object Object of class \code{semibmd}
#' returned by \code{benchmark_dose}
#'
#' @return Object of class \code{scam} or \code{gam} containing the fitted dose-response
#' model used to compute the bmd(l).
#'
#' @export
get_model <- function(object) object$model


