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
#' @param ... Not used
#'
#' @return Numeric vector of length 2 containing \code{c(bmd,bmdl)}.
#'
#' @export
get_bmd <- function(object,...) with(object,c(bmd,bmdl))

#' Get the dose-response model
#'
#' Helper function to extract the dose-response model object, of class \code{scam}
#' or \code{gam}, from an object of class \code{semibmd}
#' returned by \code{benchmark_dose}.
#'
#' @param object Object of class \code{semibmd}
#' returned by \code{benchmark_dose}
#' @param ... Not used
#'
#' @return Object of class \code{scam} or \code{gam} containing the fitted dose-response
#' model used to compute the bmd(l).
#'
#' @export
get_model <- function(object,...) object$model

#' Check the estimated BMD
#'
#' Check the value of the nonlinear equation defining the estimated benchmark dose
#' at the estimated benchmark dose. Used for debugging and diagnosing BMD estimation
#' failure.
#'
#' @param object Object of class \code{semibmd}
#' returned by \code{benchmark_dose}
#' @param ... Not used
#'
#' @return Value of the nonlinear equation \code{U(x)} evaluated at \code{x = } estimated benchmark dose.
#' If the BMD is estimated correctly, this should equal \code{0}.
#'
#' @export
get_uxb <- function(object,...) object$info$Uxb

#' Check for errors in BMD estimation
#'
#' Check whether errors were returned by the estimation procedure for the model fitting,
#' BMD or BMDL calculation.
#'
#' @param object Object of class \code{semibmd}
#' returned by \code{benchmark_dose}
#' @param verbose Logical, print any errors obtained? Default \code{FALSE}.
#' @param ... Not used
#'
#' @return Logical, were any errors obtained during the fitting process? If \code{TRUE} and \code{verbose=TRUE},
#' will print the errors
#'
#' @export
get_errors <- function(object,verbose=FALSE,...) {
  out <- FALSE
  errors <- object$info$errors
  errornames <- names(errors)
  if (length(errors) > 0) {
    out <- TRUE
    if (verbose) {
      for (i in 1:length(errors)) {
        cat("Error in",errornames[i],"estimation:",errors[i]$message,"\n")
      }
    }
  }
  out
}
