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
get_bmd <- function(object,...) {
  if (!exists("bmd",object)) stop("No BMD estimated for this object.")
  out <- object$bmd
  if (exists("bmdl",object)) {
    out <- c(out,object$bmdl)
  }
  out
}

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

#' Check the estimated BMD/L
#'
#' Check the value of the nonlinear equation defining the estimated benchmark dose (lower)
#' at the estimated benchmark dose (lower). Used for debugging and diagnosing BMD/L estimation
#' failure.
#'
#' @param object Object of class \code{semibmd}
#' returned by \code{benchmark_dose}
#' @param ... Not used
#'
#' @return Value of the nonlinear equation \code{U(x)} evaluated at \code{x = } estimated benchmark dose.
#' If the BMD is estimated correctly, this should equal \code{0}.
#' Similar for \code{Psi}.
#'
#' @rdname getuxb
#' @export
get_uxb <- function(object,...) object$info$Uxb
#' @rdname getuxb
#' @export
get_psixl <- function(object,...) object$info$Psixl

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
  if (inherits(object,'condition')) return(TRUE)
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

#' Get all benchmark dose lower limits
#'
#' Helper function to extract all BMDL's estimated in an object of class \code{semibmd}
#' returned by \code{benchmark_dose}.
#'
#' @param object Object of class \code{semibmd}
#' returned by \code{benchmark_dose}
#' @param ... Not used
#'
#' @return Named numeric vector containing the BMDLs. Currently this includes \code{score}
#' and \code{delta}.
#'
#' @export
get_all_bmdl <- function(object,...) {
  out <- numeric()
  if (exists('bmdl',object)) {
    out[1] <- get_bmd(object)[2]
    names(out) <- 'score'
  }
  possible_bmdl <- c('delta','bootstrap','bmdl_bayes')
  for (j in 1:length(possible_bmdl)) {
    nm <- possible_bmdl[j]
    if (nm %in% names(object$info$bmdl_alternatives)) {
      out <- c(out,unname(object$info$bmdl_alternatives[[nm]]))
      names(out) <- c(names(out)[names(out)!= ""],nm)
    }
  }
  out
}

#' Get quantities related to the asymptotic approximations
#'
#' Helper function to extract variance and estimated derivative used in the construction
#' of BMDLs.
#'
#' @param object Object of class \code{semibmd}
#' returned by \code{benchmark_dose}
#' @param ... Not used
#'
#' @return Named numeric vector containing the approximate quantities. Currently this includes \code{Vn} (variance of \code{U(xb)})
#' and \code{Upn} (derivative of \code{U} at estimated BMD).
#'
#' @export
get_approximations <- function(object,...) {
  out <- numeric()
  if ('delta' %in% names(object$info$bmdl_alternatives)) {
    out <- Reduce(c,object$info$approximations)
    names(out) <- names(object$info$approximations)
  }
  out
}

#' Get computation times
#'
#' Helper function to extract the computation time in seconds for each step of what's
#' estimated in an object of class \code{semibmd}
#' returned by \code{benchmark_dose}.
#'
#' @param object Object of class \code{semibmd}
#' returned by \code{benchmark_dose}
#' @param ... Not used
#'
#' @return Named numeric vector containing the computation times in seconds.
#'
#' @rdname comptimes
#'
#' @export
get_computation_times <- function(object,...) {
  times <- object$info$computation_time
  timenames <- names(times)
  ll <- length(times)
  if (ll > 0) {
    out <- numeric(ll)
    for (i in 1:ll) {
      out[i] <- times[[i]]
    }
    names(out) <- timenames
    out['total'] <- sum(out)
    return(out)
  }
  numeric()
}
#' @rdname comptimes
#'
#' @export
get_relative_computation_times <- function(object,...) {
  tms <- get_computation_times(object,...)
  m <- length(tms)
  scales::percent(tms[1:(m-1)]/tms[m],accuracy=.01)
}

# Sampling-based intervals
# Not exported.
get_samples <- function(beta,alpha,H,tmbdata,M=1000) {
  # beta: coefficients
  # alpha: intercept
  # H: full joint hessian
  # M: number of samples to return
  # Returns a d x M matrix of samples, where d = dim(H)
  Hchol <- Matrix::Cholesky(H,LDL=FALSE,perm=FALSE)
  p <- ncol(H)
  d <- length(beta)
  Z <- matrix(stats::rnorm(p*M),p,M)
  Z <- Matrix::solve(Hchol,Z,system="Lt")
  betasamps <- tmbdata$U %*% Z[1:d, ]
  betasamps <- sweep(betasamps,1,beta,'+')
  # gammasamps <- apply(betasamps,2,get_gamma)
  alphasamps <- Z[d+1, ] + alpha
  # fitted_samps <- tmbdata$X %*% gammasamps
  # colmeans <- colMeans(fitted_samps)
  # fitted_samps <- sweep(fitted_samps,2,colmeans,"-")
  # fitted_samps <- sweep(fitted_samps,2,alphasamps,"+")
  list(beta = betasamps,alpha=alphasamps)
}



