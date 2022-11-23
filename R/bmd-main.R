### Benchmark dosing ###
# Core functions for benchmark dosing
# Alex Stringer
# 11/2022


#' Reflected Newton's method
#'
#' R function implementing Newton's method for root finding with reflection
#'
#' @param g Function to find the root of. Passed to \code{rlang::as_function}.
#' @param bounds Numeric vector of length 2 containing lower and upper bounds for where the root is
#' @param eps Numerical tolerance
#'
#' @details This function solves \code{g(x)=0} for a solution satisfying \code{x>=min(bounds)},
#' \code{x<=max(bounds)}, and \code{|g(x)|<eps}.
#'
#' @return A numeric value representing the root.
#'
#' @examples
#' bounded_newton(function(x) x,c(-1,1))
#' bounded_newton(~.x,c(-1,1))
#' bounded_newton("force",c(-1,1))
#'
#' @export
bounded_newton <- function(g,bounds,eps = 1e-06) {
  stopifnot(length(bounds)==2)
  g <- rlang::as_function(g)
  gprime <- function(u) numDeriv::grad(g,u,method = "simple")
  # t <- 0
  xt <- mean(bounds) # Starting value
  while(abs(g(xt)) > eps) {
    # t <- t+1
    xt <- xt - g(xt)/gprime(xt)
    xt <- reflect(xt,bounds)
  }
  xt
}

#' Semi-parametric benchmark dosing
#'
#' Benchmark dosing based on semi-parametric (monotone) generalized additive models
#'
#' @param formula formula passed to \code{scam::scam} or \code{mgcv::gam}
#' @param data \code{data.frame} passed to \code{scam::scam} or \code{mgcv::gam}
#' @param exposure Character vector naming the exposure variable in \code{formula}
#' @param x0 Baseline exposure
#' @param p0 Baseline response
#' @param BMR Benchmark response
#' @param monotone Logical: should a monotone generalized additive model be used as the dose-response model?
#' Default \code{TRUE} uses \code{scam::scam}; setting to \code{FALSE} uses \code{mgcv::gam}.
#'
#' @export
estimate_bmd <- function(formula,data,exposure,x0=0,p0=.05,BMR=.05,monotone = TRUE) {
  ## Setup BMD related quantities ##
  A <- stats::qnorm(p0+BMR) - stats::qnorm(p0)
  ## Fit model ##
  if (monotone) {
    mod <- tryCatch(scam::scam(formula=formula,data=data,family = 'gaussian'),error = function(e) e)
  } else {
    mod <- tryCatch(mgcv::gam(formula=formula,data=data,family='gaussian'),error = function(e) e)
  }
  if (inherits(mod,'condition')) return(mod)

  # scam and gam largely share the same API for prediction etc. For the remainder of the function
  # differences are commented.

  # Set up prediction
  predframe <- data.frame(as.list(colMeans(data)))
  predframe[exposure] <- x0
  smooth_terms <- mgcv::interpret.gam(formula)$smooth.spec
  smooth_term_names <- Reduce(c,sapply(smooth_terms,function(x) x['term'],simplify = TRUE))
  which_is_bmd <- which(smooth_term_names==exposure)
  if (length(which_is_bmd)==0)
    stop(paste0("Exposure variable not found in smooth terms. exposure = ",exposure,", but smooth terms = ",smooth_term_names))
  bmd_term_label <- smooth_terms[[which_is_bmd]]$label
  # TODO: check same for scam and gam


  f0 <- stats::predict(mod,newdata = predframe,type = 'terms')[ ,bmd_term_label]
  ss <- sqrt(summary(mod)$dispersion)

  # This gives the design matrix for the whole predictor with all terms
  # This is ok, because the non-"x" terms will cancel when subtracted,
  # and not contribute to the quadratic form defining the variance
  Bmat_x0 <- stats::predict(mod,newdata = predframe,type='lpmatrix')
  # Root finding
  U <- function(x) {
    tmppredframe <- rbind(predframe,predframe)
    tmppredframe[2,exposure] <- x
    preds <- stats::predict(mod,newdata = tmppredframe,type='response')
    (1/ss)*(preds[1] - preds[2]) - A
  }
  # Confidence interval
  Psi <- function(x) {
    UU <- U(x)
    tmppredframe <- predframe
    tmppredframe[ ,exposure] <- x
    Bmat_xb <- stats::predict(mod,newdata = tmppredframe,type='lpmatrix')
    if (monotone) {
      V <- mod$Vp.t
    } else {
      V <- stats::vcov(mod,unconditional = TRUE)
    }
    varest <- (Bmat_x0 - Bmat_xb) %*% V %*% t(Bmat_x0 - Bmat_xb)/ss^2
    UU^2 - varest * stats::qchisq(.95,1)
  }

  ## Newton ##
  # Get bounds
  xmax <- max(data[ ,exposure])
  if (U(xmax) < 0) return(simpleError("U(x_max) < 0; bmd estimate does not exist"))
  bounds_u <- c(x0,xmax)

  bmd_est <- bounded_newton(U,bounds_u)

  bounds_l <- c(x0,bmd_est)
  bmd_l_est <- bounded_newton(Psi,bounds_l)

  ## Create output object ##

  out <- list(bmd = bmd_est,
              bmdl = bmd_l_est)
  class(out) <- "semibmd"
  out
}




