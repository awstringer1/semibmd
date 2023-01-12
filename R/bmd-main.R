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
#' @param maxitr Maximum number of iterations. Default is fairly small, since when the problem is bounded it usually converges quickly.
#'
#' @details This function solves \code{g(x)=0} for a solution satisfying \code{x>=min(bounds)},
#' \code{x<=max(bounds)}, and \code{|g(x)|<eps}. It uses numerical derivatives.
#'
#' @return A numeric value representing the root.
#'
#' @examples
#' bounded_newton(function(x) x,c(-1,1))
#' bounded_newton(~.x,c(-1,1))
#' bounded_newton("force",c(-1,1))
#'
#' @export
bounded_newton <- function(g,bounds,eps = 1e-06,maxitr=10) {
  stopifnot(length(bounds)==2)
  g <- rlang::as_function(g)
  gprime <- function(u) numDeriv::grad(g,u,method = "simple")
  # t <- 0
  xt <- mean(bounds) # Starting value
  itr <- 1
  while(abs(g(xt)) > eps & itr < maxitr) {
    # t <- t+1
    xt <- xt - g(xt)/gprime(xt)
    xt <- reflect(xt,bounds)
    itr <- itr+1
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
#' @param x0 Baseline exposure, default \code{0}
#' @param p0 Baseline response
#' @param BMR Benchmark response
#' @param BMDL character: what kind of BMDL to return? \code{score} gives the score-based inversion lower bound;
#' \code{delta} gives the delta method-based lower bound. \code{all} gives all available options.
#' @param monotone Logical: should a monotone generalized additive model be used as the dose-response model?
#' Default \code{TRUE} uses \code{scam::scam}; setting to \code{FALSE} uses \code{mgcv::gam}.
#' @param verbose Logical, print progress and diagnostic information for debugging? Default \code{FALSE}.
#' @inheritParams bounded_newton
#' @param boot Nonegative integer, how many bootstrap iterations to do for the parametric bootstrap BMDL? Default of \code{0}
#' means no bootstrapping.
#'
#' @note The BMDL assigned to the \code{bmdl} output object slot, and accessed via \code{get_bmd()} etc,
#' is always the score based one; the others exist only for comparison and diagnostics, and will be included
#' in the object's \code{info} slot.
#'
#' @export
benchmark_dose <- function(formula,data,exposure,x0=0,p0=.05,BMR=.05,BMDL=c("all","score","delta"),monotone = TRUE,verbose = FALSE,eps=1e-06,maxitr=100,boot=0) {
  BMDL <- BMDL[1]

  ## Create output object ##
  out <- list(
    info = list(errors = list(),bmdl_alternatives = list(),computation_time = list())
  )
  class(out) <- "semibmd"

  ## Setup BMD related quantities ##
  A <- stats::qnorm(p0+BMR) - stats::qnorm(p0)
  ## Fit model ##
  if (verbose) cat("Fitting model...\n")
  if (monotone) {
    tm <- Sys.time()
    mod <- tryCatch(scam::scam(formula=formula,data=data,family = 'gaussian'),error = function(e) e)
    dt <- as.numeric(difftime(Sys.time(),tm,units='secs'))
  } else {
    tm <- Sys.time()
    mod <- tryCatch(mgcv::gam(formula=formula,data=data,family='gaussian'),error = function(e) e)
    dt <- as.numeric(difftime(Sys.time(),tm,units='secs'))
  }
  out$info$computation_time$model <- dt
  if (inherits(mod,'condition')) {
    if (verbose) cat("Received the following error when fitting model:",mod$message,".\n")
    out$info$errors$model <- mod
    return(out)
  }
  out$model <- mod
  if (verbose) cat("Finished fitting model.\n")

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
  check_Umax <- tryCatch(U(xmax),error=function(e)e)
  if (inherits(check_Umax,'condition')) {
    if (verbose) cat(check_Umax$message,".\n")
    out$info$errors$Umax <- check_Umax
    return(out)
  }
  if (check_Umax < 0) {
    e <- simpleError("Error: U(x_max) < 0; bmd estimate does not exist. Please make a note of it.")
    if (verbose) cat(e$message,".\n")
    out$info$errors$Umax <- e
    return(out)
  }
  bounds_u <- c(x0,xmax)

  if (verbose) cat("Running Newton for BMD...\n")
  tm <- Sys.time()
  bmd_est <- tryCatch(bounded_newton(U,bounds_u,eps=eps,maxitr=maxitr),error = function(e) e)
  dt <- as.numeric(difftime(Sys.time(),tm,units='secs'))
  out$info$computation_time$bmd <- dt
  if (inherits(bmd_est,'condition')) {
    if (verbose) cat("Received the following error when estimating BMD:",bmd_est$message,".\n")
    out$info$errors$bmd <- bmd_est
    return(out)
  }
  out$bmd <- bmd_est
  out$info$Uxb <- U(bmd_est)
  if (verbose) cat("Finished running Newton for BMD.\n")

  if (BMDL=="none") return(out)

  if (BMDL %in% c("all","score")) {
    bounds_l <- c(x0,bmd_est)
    if (verbose) cat("Running Newton for BMDL...\n")
    tm <- Sys.time()
    bmd_l_est <- tryCatch(bounded_newton(Psi,bounds_l,eps=eps,maxitr=maxitr),error = function(e) e)
    dt <- as.numeric(difftime(Sys.time(),tm,units='secs'))
    out$info$computation_time$bmdl_score <- dt
    if (inherits(bmd_l_est,'condition')) {
      if (verbose) cat("Received the following error when estimating BMD:",bmd_l_est$message,".\n")
      out$info$errors$bmdl <- bmd_l_est
      return(out)
    }
    out$bmdl <- bmd_l_est
    out$info$Psixl <- Psi(bmd_l_est)
    if (verbose) cat("Finished running Newton for BMD.\n")
  }

  if (BMDL %in% c("all","delta")) {
    tm <- Sys.time()
    tmppredframe <- predframe
    tmppredframe[ ,exposure] <- bmd_est
    Bmat_xb <- stats::predict(mod,newdata = tmppredframe,type='lpmatrix')
    if (monotone) {
      V <- mod$Vp.t
    } else {
      V <- stats::vcov(mod,unconditional = TRUE)
    }
    Vn <- (Bmat_x0 - Bmat_xb) %*% V %*% t(Bmat_x0 - Bmat_xb)/ss^2
    Upn <- abs(numDeriv::grad(U,bmd_est))
    bmd_l_delta_est <- bmd_est - stats::qnorm(.975)*sqrt(Vn)/Upn
    dt <- as.numeric(difftime(Sys.time(),tm,units='secs'))
    out$info$computation_time$bmdl_delta <- dt
    out$info$bmdl_alternatives$delta <- bmd_l_delta_est
    out$info$approximations$Vn <- Vn
    out$info$approximations$Upn <- Upn
  }

  ## Bootstrapping ##
  if (boot > 0) {
    # Generate a new dataset from the fitted model
    # Calculate the BMD (no BMDL)
    # Do this "boot" times
    # Return lower 2.5%ile of these as the bootstrapped BMDL

    bootbmd <- numeric(boot)
    n <- nrow(data)
    predmean <- stats::predict(mod,type='response')
    response_var <- mgcv::interpret.gam(formula)$response
    if (verbose) cat ("Bootstrapping with B = ",boot," samples...\n",sep="")
    tm <- Sys.time()
    for (b in 1:boot) {
      # Generate a new dataset
      newdat <- data
      newdat[ ,response_var] <- stats::rnorm(n,predmean,ss)
      # Fit the model
      bootmod <- benchmark_dose(
        formula = formula,
        data = newdat,
        exposure = exposure,
        x0 = x0,p0 = p0,BMR = BMR,
        BMDL = "none",
        monotone = monotone,verbose = FALSE,eps = eps,maxitr = maxitr,boot = 0
      )
      # Get the bmd
      bootbmd[b] <- get_bmd(bootmod)[1]
    }
    dt <- as.numeric(difftime(Sys.time(),tm,units='secs'))
    out$info$computation_time$bmdl_boot <- dt
    if (verbose) cat("Finished bootstrapping.\n")
    out$info$bmdl_alternatives$bootstrap <- unname(stats::quantile(bootbmd,.025))
  }

  out
}

#' Summary method for \code{semibmd} objects
#'
#' Summary method for \code{semibmd} objects returned by \code{semibmd::benchmark_dose}.
#' Combines the summary output of \code{scam/gam} with the benchmark dosing information
#' computed by \code{benchmark_dose}.
#'
#' @param object Object of class \code{semibmd} returned by \code{semibmd::benchmark_dose}.
#' @param x Object of class \code{summary.semibmd} to print.
#' @param digits Number of digits to round bmd(l) for printing.
#' @param ... Not used.
#'
#' @details This function first calls the appropriate \code{summary} method for the
#' dose-response model, and then adds information about the benchmark dose and BMDL.
#'
#' @return Object of class \code{summary.semibmd}: a list with summary information, for printing
#'
#' @rdname summary.semibmd
#'
#' @export
summary.semibmd <- function(object,...) {
  modsummary <- list()
  modsummary$modelsummary <- summary(object$model)
  thebmd <- get_bmd(object)
  modsummary$bmd <- thebmd[1]
  modsummary$bmdl <- thebmd[2]
  class(modsummary) <- 'summary.semibmd'
  modsummary
}
#' @rdname summary.semibmd
#' @export
print.summary.semibmd <- function(x,digits = 4,...) {
  cat("---\n")
  cat("Dose-response model summary:\n")
  cat("---\n")
  print(x$modelsummary)
  cat("---\n")

  cat("Benchmark dose summary:\n")
  cat("---\n")
  print(data.frame(bmd = round(x$bmd,digits),bmdl = round(x$bmdl,digits)))
  cat("---\n")
}

#' Plot method for \code{semibmd} objects
#'
#' Plot method for \code{semibmd} objects returned by \code{semibmd::benchmark_dose}.
#' Combines the standard plot from \code{scam/gam} with the benchmark dosing information
#' computed by \code{benchmark_dose} and \code{summary.semibmd}.
#'
#' @param x Object of class \code{semibmd} returned by \code{semibmd::benchmark_dose}.
#' @param ... Not used.
#'
#' @details This function first calls \code{summary.semibmd}, then calls \code{plot.scam/gam}
#' as appropriate and modifies the resulting plot with the benchmark dosing information.
#'
#' @return Prints a plot with the dose-response and benchmark dose information on it.
#'
#' @rdname plot.semibmd
#'
#' @export
plot.semibmd <- function(x,...) {
  mod <- get_model(x)
  bmd <- get_bmd(x)
  plotinfo <- plot(mod)
  graphics::abline(v = bmd[1],lty='dotdash')
  graphics::abline(v = bmd[2],lty='dotdash')
}
