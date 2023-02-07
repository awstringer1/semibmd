### BMD estimation with TMB and Rcpp backend ###

#' Semi-parametric benchmark dosing
#'
#' Benchmark dosing based on a custom implementation of monotone B-spline-based additive models
#' using TMB and Rcpp
#'
#' @param monosmooths A list of \code{mgcv}-compatible formulas of the type s(x,bs='bs') to be treated as monotone smooths. Must
#' use B-Spline bases. Must contain the variable specified in \code{exposure}.
#' @param smooths A list of \code{mgcv}-compatible formulas of the type s(x,bs='bs') to be treated as non-monotone smooths. Must
#' use B-Spline bases.
#' @param data \code{data.frame} containing the variables in \code{smooths}, \code{monosmooths}, and \code{response}.
#' @param exposure Character vector naming the exposure variable in the \code{monosmooths} formulas. Will always be treated as monotone decreasing.
#' @param response Character vector naming the response variable in \code{data}.
#' @param x0 Baseline exposure, default \code{0}
#' @param p0 Baseline response
#' @param BMR Benchmark response
#' @param verbose Logical, print progress and diagnostic information for debugging? Default \code{FALSE}.
#' @param eps Numerical tolerance for Newton's method
#' @param maxitr Maximum number of iterations for Newton's method.
#' @param bayes_boot Nonegative integer, how many Bayesian "bootstrap" iterations to do; means how many posterior draws to base
#' sample-based inferences for the BMD(L) off of.
#' Must be a positive integer.
#'
#' @note The BMDL assigned to the \code{bmdl} output object slot, and accessed via \code{get_bmd()} etc,
#' is based on Bayesian bootstrapping
#'
#' @importFrom mgcv s
#'
#' @rawNamespace useDynLib(semibmd, .registration=TRUE); useDynLib(semibmd_TMBExports)
#'
#' @export
benchmark_dose_tmb <- function(monosmooths,smooths,data,exposure,response,x0,p0,BMR,verbose=FALSE,eps=1e-06,maxitr=10,bayes_boot=1e03) {
  ## Create output object ##
  out <- list(
    info = list(errors = list(),bmdl_alternatives = list(),computation_time = list())
  )
  class(out) <- "semibmd"

  ## Setup BMD related quantities ##
  A <- stats::qnorm(p0+BMR) - stats::qnorm(p0)

  ## Setup Model Fitting Quantities ##
  # NON-monotone smooths
  # TODO
  smoothobj <- list()
  if (!is.null(smooths)) {
    # TODO: list processing for multiple variables
    smoothobj <- mgcv::smoothCon(smooths[[1]],data=data,absorb.cons=FALSE,scale.penalty=TRUE)[[1]]
    S <- smoothobj$S[[1]]
    X <- smoothobj$X
    cc <- colSums(X)
    EE <- eigen(S)
    r <- sum(EE$values>1e-08)
    Dp <- EE$values[1:r]
  }

  # Monotone smooths
  if (is.null(monosmooths)) stop("At least one monotone smoothing variable must be specified.")
  monosmoothobj <- mgcv::smoothCon(monosmooths[[1]],data=data,absorb.cons=FALSE,scale.penalty=TRUE)[[1]]
  Smono <- monosmoothobj$S[[1]]
  Xmono <- monosmoothobj$X
  ccmono <- colSums(Xmono)
  EEmono <- eigen(Smono)
  rmono <- sum(EEmono$values>1e-08)
  dmono <- length(EEmono$values)
  Dpmono <- EEmono$values[1:rmono]
  Umono <- EEmono$vectors

  # TODO: properly add the non-monotone smooths
  tmbdata <- list(
    y = data[[response]],
    xcov = data[[exposure]],
    X = Xmono,
    S = Smono, # Only required for starting values
    Dp = methods::as(methods::as(Matrix::Diagonal(x=Dpmono),'generalMatrix'),'TsparseMatrix'),
    U = Umono,
    c = ccmono,
    smoothobj = monosmoothobj
  )

  # Parameters for TMB
  alpha <- mean(tmbdata$y)
  vr <- stats::var(tmbdata$y)
  sm <- 1/mean(Dpmono)
  beta <- as.numeric(with(tmbdata,Matrix::solve((Matrix::crossprod(X) + sm*S),Matrix::crossprod(X,y))))

  UR <- tmbdata$U[ ,1:rmono]
  UF <- tmbdata$U[ ,(rmono+1):dmono]
  betaR <- as.numeric(Matrix::crossprod(UR,beta))
  betaF <- as.numeric(Matrix::crossprod(UF,beta))

  tmbparams <- 	list(
    alpha = alpha,
    betaR = betaR,
    betaF = betaF,
    logprec = -log(vr),
    logsmoothing = log(sm)
  )

  template_inner <- tryCatch(TMB::MakeADFun(
    data = c(model = "monotonesmoothing",tmbdata),
    parameters = tmbparams,
    silent = TRUE,
    DLL = "semibmd_TMBExports",
    random = c("alpha","betaR","betaF")
  ),error = function(e) e)
  if (inherits(template_inner,'condition')) {
    if (verbose) cat("Received the following error when creating TMB template:",template_inner$message,".\n")
    out$info$errors$tmb_template <- template_inner
    return(out)
  }

  opt1 <- tryCatch(with(template_inner,stats::optim(par,fn,gr,method='BFGS',hessian = FALSE)),error = function(e) e)

  if (inherits(opt1,'condition')) {
    if (verbose) cat("Received the following error optimizing TMB template:",opt1$message,".\n")
    out$info$errors$optimization <- opt1
    return(out)
  }

  sigmaest <- exp(-.5*opt1$par[1])
  lambdaest <- exp(opt1$par[2])

  val <- tryCatch(template_inner$fn(opt1$par),error = function(e) e)
  if (inherits(val,'condition')) {
    if (verbose) cat("Received the following error evaluating optimized TMB template:",val$message,".\n")
    out$info$errors$optimization <- val
    return(out)
  }

  if (is.na(val) | is.nan(val)) {
    if (verbose) cat("Evaluating optimized TMB template, got value",val,".\n")
    return(out)
  }

  randest <- with(template_inner$env,last.par[random])
  betaRFest <- randest[names(randest) %in% c("betaR","betaF")]
  betaest <- as.numeric(tmbdata$U %*% betaRFest)
  gammaest <- get_gamma(betaest)
  alphaest <- randest['alpha']

  ## Standard Errors ##
  fullprec <- tryCatch(TMB::sdreport(template_inner,getJointPrecision=TRUE)$jointPrecision,error = function(e) e)
  if (inherits(fullprec,'condition')) {
    if (verbose) cat("Received the following error when obtaining the joint Hessian:",fullprec$message,".\n")
    out$info$errors$jointhessian <- fullprec
    return(out)
  }
  randomidx <- template_inner$env$random
  paramdimfull <- ncol(fullprec)
  paramdimrandom <- length(randomidx)
  hyperidx <- (paramdimrandom+1):(paramdimfull)
  randprec <- fullprec[randomidx,randomidx]

  samps <- tryCatch(get_samples(betaest,alphaest,randprec,tmbdata,M=bayes_boot),error = function(e) e)
  if (inherits(samps,'condition')) {
    if (verbose) cat("Received the following error when drawing posterior samples:",samps$message,".\n")
    out$info$errors$posteriorsamples <- samps
    return(out)
  }
  # Get empirical coverage of samples
  # Compute upper and lower quantiles
  samp_lower <- tryCatch(apply(samps$fitted,1,stats::quantile,probs=.025),error = function(e) e)
  samp_upper <- tryCatch(apply(samps$fitted,1,stats::quantile,probs=.975),error = function(e) e)
  if (inherits(samp_lower,'condition')) {
    if (verbose) cat("Received the following error when computing lower quantiles:",samp_lower$message,".\n")
    out$info$errors$lower_quantiles <- samp_lower
    out$info$samps <- samps # Return stuff if there was an error in them
    out$info$betaest <- betaest
    out$info$alphaest <- alphaest
    out$info$randprec <- randprec
    out$info$tmbdata <- tmbdata
    return(out)
  }
  if (inherits(samp_upper,'condition')) {
    if (verbose) cat("Received the following error when computing upper quantiles:",samp_upper$message,".\n")
    out$info$errors$lower_quantiles <- samp_upper
    out$info$samps <- samps # Return stuff if there was an error in them
    out$info$betaest <- betaest
    out$info$alphaest <- alphaest
    out$info$randprec <- randprec
    out$info$tmbdata <- tmbdata
    return(out)
  }

  xmax <- max(data[[exposure]])
  bmd_samps <- numeric(ncol(samps$beta))
  for (b in 1:ncol(samps$beta)) {
    tmp <- tryCatch(get_bmd_cpp(samps$beta[ ,b],tmbdata$smoothobj$knots,c(x0,xmax),x0,sigmaest,A,1e-06,100),error = function(e) e)
    if (inherits(tmp,'condition')) {
      bmd_samps[b] <- -1
    } else {
      bmd_samps[b] <- tmp
    }
  }
  bmd_samp_errors <- sum(bmd_samps == -1)
  bmd_samps_clean <- bmd_samps[bmd_samps > -1]

  bmd_est <- tryCatch(stats::median(bmd_samps_clean),error = function(e) e)
  if (inherits(bmd_est,'condition')) {
    if (verbose) cat("Received the following error when estimating BMD:",samp_upper$bmd_est,".\n")
    out$info$errors$bmd_est <- bmd_est
    return(out)
  }
  out$bmd <- bmd_est

  bmdl_est <- tryCatch(stats::quantile(bmd_samps_clean,probs=.025),error = function(e) e)
  if (inherits(bmdl_est,'condition')) {
    if (verbose) cat("Received the following error when computing BMDL:",bmdl_est$message,".\n")
    out$info$errors$lower_quantiles <- bmdl_est
    out$info$bmd_samps <- bmd_samps
    out$info$bmd_samps_clean <- bmd_samps_clean
    out$info$samps <- samps # Return stuff if there was an error in them
    out$info$betaest <- betaest
    out$info$alphaest <- alphaest
    out$info$randprec <- randprec
    out$info$tmbdata <- tmbdata
    return(out)
  }
  out$bmdl <- bmdl_est

  # TODO: add in plot and summary information about the dose-response model

  out
}

