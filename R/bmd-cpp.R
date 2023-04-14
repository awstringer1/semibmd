### BMD estimation with TMB and Rcpp backend ###

#' Semi-parametric benchmark dosing
#'
#' Benchmark dosing based on a custom implementation of monotone B-spline-based additive models
#' using TMB and Rcpp
#'
#' @param monosmooths A list of \code{mgcv}-compatible smooth formula terms of the type \code{s(x,bs='bs')} to be treated as monotone smooths. Must
#' use B-Spline bases. Must contain the variable specified in \code{exposure}.
#' @param smooths A list of \code{mgcv}-compatible smooth formula terms of the type \code{s(x,bs='bs')} to be treated as non-monotone smooths. Must
#' use B-Spline bases.
#' @param linearterms A one-sided formula of the form \code{~x1+x2} compatible with \code{model.matrix}, defining any linear terms in the model.
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
#' @param scale_data Logical, default \code{FALSE}: should the response variable be scaled to have mean 0 and standard deviation 1? See Details.
#'
#' @note The BMDL assigned to the \code{bmdl} output object slot, and accessed via \code{get_bmd()} etc,
#' is based on Bayesian bootstrapping
#'
#' @details Scaling the response variable is a
#' typical operation in regression with continuous response, but it does change the estimation of the BMD(L).
#' It is recommended to input your data on whatever scale is scientifically meaningful, and then choose
#' \code{scale_data=TRUE} to scale the data internally and account for this operation in the estimation of BMD(L).
#'
#' @importFrom mgcv s
#'
#' @rawNamespace useDynLib(semibmd, .registration=TRUE); useDynLib(semibmd_TMBExports)
#'
#' @export
benchmark_dose_tmb <- function(monosmooths,
                               smooths,
                               linearterms,
                               data,
                               exposure,
                               response,
                               x0,
                               p0,
                               BMR,
                               verbose=FALSE,
                               eps=1e-06,
                               maxitr=10,
                               bayes_boot=1e03,
                               scale_data = FALSE
) {
  ## Create output object ##
  out <- list(
    info = list(errors = list(),bmdl_alternatives = list(),computation_time = list(),data = data)
  )
  class(out) <- "semibmd"

  tm <- Sys.time()
  ## Setup BMD related quantities ##
  A <- stats::qnorm(p0+BMR) - stats::qnorm(p0)


  p <- 4 # B-Spline order

  if (scale_data) {
    mean_response <- mean(data[[response]])
    sd_response <- stats::sd(data[[response]])
    data[[response]] <- (data[[response]] - mean_response)/sd_response
  }

  ## Setup Model Fitting Quantities ##
  # NON-monotone smooths
  # TODO
  smoothobj <- S <- X <- cc <- EE <- r <- d <- Dp <- U <- list()
  nonmono <- 0
  if (!is.null(smooths)) {
    nonmono <- 1L
    for (j in 1:length(smooths)) {
      # TODO: list processing for multiple variables
      smoothobj[[j]] <- mgcv::smoothCon(smooths[[j]],data=data,absorb.cons=TRUE,scale.penalty=TRUE)
      S[[j]] <- Matrix::bdiag(Map('[[',Map('[[',smoothobj[[j]],'S'),1))
      X[[j]] <- Reduce(cbind,Map('[[',smoothobj[[j]],'X'))
      cc[[j]] <- colSums(X[[j]])
      EE[[j]] <- eigen(S[[j]])
      r[[j]] <- sum(EE[[j]]$values>1e-08)
      d[[j]] <- length(EE[[j]]$values)
      Dp[[j]] <- EE[[j]]$values[1:r[[j]]]
      U[[j]] <- EE[[j]]$vectors
    }
    ## Create the combined smooth matrices ##
    # X: column binded. Also use to create indices
    Xsmooth <- Reduce(cbind,X)
    dsmooth <- Reduce(c,Map(ncol,X))
    rsmooth <- Reduce(c,r)
    # U: block diagonal, sparse
    Usmooth <- Matrix::bdiag(U)
    Dpsmooth <- Reduce(c,Dp)
  } else {
    rsmooth <- dsmooth <- 1 # Placeholder
  }


  # Monotone smooths
  if (is.null(monosmooths)) stop("At least one monotone smoothing variable must be specified.")
  if (length(monosmooths)>1) stop("Only one monotone smooth supported at this time.")
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
    Xmono = Xmono,
    c = colSums(Xmono),
    S = Smono, # Only required for starting values
    Dpmono = methods::as(methods::as(Matrix::Diagonal(x=Dpmono),'generalMatrix'),'TsparseMatrix'),
    Umono = Umono,
    smoothobj = monosmoothobj,
    nonmono = nonmono
  )
  if (nonmono)
      tmbdata <- c(tmbdata,list(
        Xsmooth = Xsmooth,
        Dpsmooth = methods::as(methods::as(Matrix::Diagonal(x=Dpsmooth),'generalMatrix'),'TsparseMatrix'),
        Usmooth = Usmooth,
        rs = rsmooth,
        ds = dsmooth
      ))

  # Starting values for TMB
  alpha <- mean(tmbdata$y)
  vr <- stats::var(tmbdata$y)
  sm <- 1/mean(EEmono$values[1:rmono])
  betamono <- as.numeric(with(tmbdata,Matrix::solve((Matrix::crossprod(Xmono) + sm*Smono),Matrix::crossprod(Xmono,y))))

  URmono <- tmbdata$Umono[ ,1:rmono]
  UFmono <- tmbdata$Umono[ ,(rmono+1):dmono]
  betaRmono <- as.numeric(Matrix::crossprod(URmono,betamono))
  betaFmono <- as.numeric(Matrix::crossprod(UFmono,betamono))
  # Unconstrained smooths
  betasmooth <- numeric(sum(dsmooth))
  betaRsmooth <- numeric(sum(rsmooth))
  betaFsmooth <- numeric(sum(dsmooth-rsmooth))
  logsmoothingsmooth <- numeric(length(rsmooth))
  if (nonmono) {
    betaRsmooth <- betaFsmooth <- c()
    for (j in 1:length(smooths)) {
      logsmoothingsmooth[j] <- -log(mean(EE[[j]]$values[1:r[[j]]]))
      tmpbeta <- as.numeric(Matrix::solve((Matrix::crossprod(X[[j]]) + exp(logsmoothingsmooth[j])*S[[j]]),Matrix::crossprod(X[[j]],tmbdata$y)))
      tmpUR <- EE[[j]]$vectors[ ,1:r[[j]]]
      tmpUF <- EE[[j]]$vectors[ ,(r[[j]]+1):d[[j]]]
      betaRsmooth <- c(betaRsmooth,as.numeric(Matrix::crossprod(tmpUR,tmpbeta)))
      betaFsmooth <- c(betaFsmooth,as.numeric(Matrix::crossprod(tmpUF,tmpbeta)))
    }
  }

  tmbparams <- 	list(
    alpha = alpha,
    betaRmono = betaRmono,
    betaFmono = betaFmono,
    logprec = -log(vr),
    logsmoothingmono = log(sm),
    betaRsmooth = betaRsmooth,
    betaFsmooth = betaFsmooth,
    logsmoothingsmooth = logsmoothingsmooth
  )

  if (nonmono) {
    whichmodel <- "monotonesmoothingmulti"
    whichrandom <- c("alpha","betaRmono","betaFmono","betaRsmooth","betaFsmooth")
  } else {
    whichmodel <- "monotonesmoothing"
    whichrandom <- c("alpha","betaRmono","betaFmono")
  }

  template_inner <- tryCatch(TMB::MakeADFun(
    data = c(model = whichmodel,tmbdata),
    parameters = tmbparams,
    silent = TRUE,
    DLL = "semibmd_TMBExports",
    random = whichrandom
    # random = NULL,
    # map = list(logsmoothingsmooth=factor(rep(NA,2)),logprec=factor(NA),logsmoothingmono=factor(NA))
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
  # if (scale_data) sigmaest <- sigmaest*sd_response
  lambdaest <- exp(opt1$par[2])
  if (nonmono) {
    lambdasmoothest <- exp(opt1$par[3:length(opt1$par)])
  }

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

  ## Estimates ##
  ## Everything here still pertains only to the monotone smooth ##
  randest <- with(template_inner$env,last.par[random])
  alphaest <- randest['alpha']
  # if (scale_data) {
  #   alphaest <- sd_response*alphaest+mean_response # Add the mean back onto the intercept, only
  #   randest <- sd_response*randest
  # }
  betaRFest <- randest[names(randest) %in% c("betaRmono","betaFmono")]
  betaest <- as.numeric(tmbdata$Umono %*% betaRFest)
  gammaest <- get_gamma(betaest) # Do this after scaling

  ## Standard Errors ##
  fullprec <- tryCatch(TMB::sdreport(template_inner,getJointPrecision=TRUE)$jointPrecision,error = function(e) e)
  if (inherits(fullprec,'condition')) {
    if (verbose) cat("Received the following error when obtaining the joint Hessian:",fullprec$message,".\n")
    out$info$errors$jointhessian <- fullprec
    return(out)
  }
  if (scale_data) {
    # divide PRECISION matrix by square of sd_response, i.e. multiply variance matrix by variance of response.
    # fullprec <- (1/sd_response^2)*fullprec
  }

  randomidx <- template_inner$env$random
  paramdimfull <- ncol(fullprec)
  paramdimrandom <- length(randomidx)
  hyperidx <- (paramdimrandom+1):(paramdimfull)
  randprec <- fullprec[randomidx,randomidx]

  # UPDATE: index out the monotone smooths and the intercept, not all random effects
  monoidx <- which(names(randest) %in% c("betaRmono","betaFmono","alpha"))
  # This works better if the monotone smooths are the first ones in the vector, followed by the intercept
  # This is determined only by the implementation, not the user, but check for it to prevent a mistake by future me
  # TODO: use proper indexing, will need to do for plotting all the curves
  stopifnot(all(monoidx == 1:length(monoidx)))
  dt <- as.numeric(difftime(Sys.time(),tm,units='secs'))
  out$info$computation_time$model <- dt

  tm <- Sys.time()
  # samps <- tryCatch(get_samples(betaest,alphaest,randprec,tmbdata,M=bayes_boot),error = function(e) e)
  if (nonmono) {
    nonmonocoefest <- randest[names(randest) %in% c("betaRsmooth","betaFsmooth")]
    samps <- tryCatch(get_samples(betaest,alphaest,randprec,tmbdata,M=bayes_boot,nonmonocoef=nonmonocoefest),error = function(e) e)
  } else {
    samps <- tryCatch(get_samples(betaest,alphaest,randprec,tmbdata,M=bayes_boot),error = function(e) e)
  }
  if (inherits(samps,'condition')) {
    if (verbose) cat("Received the following error when drawing posterior samples:",samps$message,".\n")
    out$info$errors$posteriorsamples <- samps
    return(out)
  }
  dt <- as.numeric(difftime(Sys.time(),tm,units='secs'))
  out$info$computation_time$posterior_samples <- dt

  # Posterior of BMD
  xmax <- max(data[[exposure]])
  bmd_samps <- numeric(ncol(samps$beta))
  tm <- Sys.time()
  for (b in 1:ncol(samps$beta)) {
    tmp <- tryCatch(get_bmd_cpp(samps$beta[ ,b],tmbdata$smoothobj$knots,c(x0,xmax),x0,sigmaest,A,1e-06,100),error = function(e) e)
    if (inherits(tmp,'condition')) {
      bmd_samps[b] <- -1
    } else if (tmp < x0) {
      # if get_bmd_cpp returns less than lower bound, it's an error
      bmd_samps[b] <- -1 
    } else {
      bmd_samps[b] <- tmp
    }
  }
  bmd_samp_errors <- sum(bmd_samps == -1)
  bmd_samps_clean <- bmd_samps[bmd_samps > -1 & !is.na(bmd_samps) & !is.nan(bmd_samps)]
  out$info$bmd_samps <- bmd_samps
  dt <- as.numeric(difftime(Sys.time(),tm,units='secs'))
  out$info$computation_time$bmd_samples <- dt

  tm <- Sys.time()
  bmd_est <- tryCatch(get_bmd_cpp(betaest,tmbdata$smoothobj$knots,c(x0,xmax),x0,sigmaest,A,1e-06,100),error = function(e) e)
  if (inherits(bmd_est,'condition')) {
    if (verbose) cat("Received the following error when estimating BMD:",bmd_est$bmd_est,".\n")
    out$info$errors$bmd_est <- bmd_est
    out$info$bmd_samps <- bmd_samps
    out$info$bmd_samps_clean <- bmd_samps_clean
    return(out)
  } else if (bmd_est < x0) {
    if (verbose) cat("Received error value",bmd_est,"from get_bmd_cpp when estimating BMD.\n")
    out$info$errors$bmd_est <- bmd_est
    out$info$bmd_samps <- bmd_samps
    out$info$bmd_samps_clean <- bmd_samps_clean
    return(out)

  }
  out$bmd <- bmd_est
  dt <- as.numeric(difftime(Sys.time(),tm,units='secs'))
  out$info$computation_time$bmd_estimate <- dt

  ## BMDL several ways ##

  bmdl_est_bayesboot <- tryCatch(stats::quantile(bmd_samps_clean,probs=.025),error = function(e) e)
  if (inherits(bmdl_est_bayesboot,'condition')) {
    if (verbose) cat("Received the following error when computing BMDL:",bmdl_est_bayesboot$message,".\n")
    out$info$errors$lower_quantiles <- bmdl_est_bayesboot
    out$info$bmd_samps <- bmd_samps
    out$info$bmd_samps_clean <- bmd_samps_clean
    out$info$samps <- samps # Return stuff if there was an error in them
    out$info$betaest <- betaest
    out$info$alphaest <- alphaest
    out$info$fullprec <- fullprec
    out$info$tmbdata <- tmbdata
    return(out)
  }
  out$info$bmdl_alternatives$bmdl_bayes <- bmdl_est_bayesboot

  ## Delta ##
  tm <- Sys.time()
  gammasamps <- apply(samps$beta,2,get_gamma)
  # gammameans <- colMeans(gammasamps)
  # gammasampsC <- sweep(gammasamps,2,gammameans,'-')
  V <- stats::cov(t(gammasamps)) # This is surprisingly fast
  # V <- stats::cov(t(gammasampsC)) # This is surprisingly fast

  kx <- knotindex(bmd_est,tmbdata$smoothobj$knots)
  bx0 <- Bsplinevec(x0,tmbdata$smoothobj$knots,4)
  Vn <- Vx_cpp(bmd_est,V,tmbdata$smoothobj$knots,bx0,sigmaest)
  # if (scale_data) Vn <- Vn * (sd_response^2)
  dt_tmp <- as.numeric(difftime(Sys.time(),tm,units='secs'))

  tm <- Sys.time()
  gammadiff <- (p-1)*c(0,diff(gammaest)[1:(length(gammaest)-1)]) / (tmbdata$smoothobj$knots[(p+1):(length(gammaest)+p)] - tmbdata$smoothobj$knots[2:(length(gammaest)+1)])
  #Upn <- abs(Uxd_cpp(bmd_est,gammadiff,tmbdata$smoothobj$knots,kx,sigmaest))
  Upn <- abs(Uxd_cpp(bmd_est,gammaest,tmbdata$smoothobj$knots,kx,sigmaest))
  bmd_l_delta_est <- bmd_est - stats::qnorm(.975)*sqrt(Vn)/Upn
  dt <- as.numeric(difftime(Sys.time(),tm,units='secs'))
  out$info$computation_time$bmdl_delta <- dt + dt_tmp
  out$info$bmdl_alternatives$delta <- bmd_l_delta_est
  out$info$approximations$Vn <- Vn
  out$info$approximations$Upn <- Upn
  ## Score ##
  tm <- Sys.time()
  bmd_l_score_est <- tryCatch(get_score_cpp(betaest,V,tmbdata$smoothobj$knots,c(x0,bmd_est),x0,sigmaest,A,1e-06,10),error = function(e) e)
  dt <- as.numeric(difftime(Sys.time(),tm,units='secs'))
  out$info$computation_time$bmdl_score <- dt + dt_tmp
  if (inherits(bmd_l_score_est,'condition')) {
    if (verbose) cat("Received the following error when estimating score BMDL:",bmd_l_score_est$bmd_l_score_est,".\n")
    out$info$errors$bmd_l_score_est <- bmd_l_score_est
    return(out)
  } else if (bmd_l_score_est < x0) {
    if (verbose) cat("Received error value",bmd_l_score_est,"when estimating BMDL via get_score_cpp.\n")
    out$info$errors$bmd_l_score_est <- bmd_l_score_est
    return(out)
  }
  out$bmdl <- bmd_l_score_est

  # Plot information
  tm <- Sys.time()
  out$model <- list(
    plotinfo = list(
      minx = min(data[[exposure]]),
      maxx = max(data[[exposure]]),
      monosmoothobj = monosmoothobj,
      samps = samps
    )
  )
  if (nonmono) {
    out$model$plotinfosmooth <- list()
    for(j in 1:length(smooths)) {
      tmpsmooth <- smooths[[j]]
      tmpvar <- tmpsmooth$term
      out$model$plotinfosmooth[[j]] <- list(
        minx = min(data[[tmpvar]]),
        maxx = max(data[[tmpvar]]),
        smoothobj = smoothobj[[j]],
        U = U[[j]],
        r = r[[j]],
        d = d[[j]]
      )
    }
  }
  dt <- as.numeric(difftime(Sys.time(),tm,units='secs'))
  out$info$computation_time$plot_information <- dt

  # Add the U(x) and Psi(x) functions to the output
  Ux <- function(x) {
    knots <- tmbdata$smoothobj$knots
    kx <- knotindex(x,knots)
    k0 <- knotindex(x0,knots)
    fx0 <- deBoor(x0,k0,knots,gammaest,4)
    Ux_cpp(x,gammaest,knots,kx,fx0,sigmaest,A)
  }
  Uxv <- function(x) {
    out <- numeric(length(x))
    for (i in 1:length(out)) out[i] <- Ux(x[i])
    out
  }

  Psix <- function(x) {
    knots <- tmbdata$smoothobj$knots
    kx <- knotindex(x,knots)
    k0 <- knotindex(x0,knots)
    fx0 <- deBoor(x0,k0,knots,gammaest,4)
    bx0 <- Bsplinevec(x0,knots,4)

    Psix_cpp(x,gammaest,V,knots,kx,fx0,bx0,sigmaest,A)
  }
  Psixv <- function(x) {
    out <- numeric(length(x))
    for (i in 1:length(out)) out[i] <- Psix(x[i])
    out
  }
  out$info$functions <- list(Ux = Uxv,Psix = Psixv)

  out
}

