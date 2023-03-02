# devtools::load_all()
# Setup file for semibmd unit tests
set_parameters <- function(xmin,xmax,x0=.05,sigma=1,p0=.025,BMR=.1) {
  # True function f
  f <- function(x) 1/sqrt(x) # Nonlinear, monotonic, invertible
  finv <- function(x) 1/(x^2) # For calculating true benchmark dose

  # x0 <- xmin + (xmax-xmin)*.1 # Reference dose
  tau0 <- f(x0) + sigma * qnorm(p0) # Response level
  stopifnot(abs(pnorm(tau0,f(x0),sigma) - p0) < .Machine$double.eps) # Definition of tau

  A <- qnorm(p0+BMR) - qnorm(p0)
  # Benchmark dose
  xb <- finv(f(x0) - sigma*A)
  stopifnot(abs(pnorm(tau0,f(xb),sigma) - (p0+BMR)) < .Machine$double.eps) # Definition of benchmark dose

  # Check the derivative at the true BMD
  fdxb <- numDeriv::grad(f,xb)

  list(
    xmin = xmin,xmax = xmax,
    sigma = sigma,
    p0 = p0,BMR = BMR,
    f = f,
    x0 = x0,
    tau0 = tau0,
    A = A,
    xb = xb,fdxb = fdxb
  )
}

## Simulate Data ##
simulate_data <- function(n,params) {
  xcov <- with(params,seq(xmin,xmax,length.out=n))
  x1 <- with(params,runif(n,xmin,xmax))
  x2 <- with(params,runif(n,xmin,xmax))

  with(params,data.frame(y = rnorm(n,f(xcov)+2*x1-x2,sigma),x = xcov,x1=x1,x2=x2))
}
params <- set_parameters(.05,.2)
set.seed(3798)
dat <- simulate_data(100,params)

## Benchmark dosing ##
# formula <- y~s(x,bs='mpd')
# data <- dat
# exposure <- 'x'
# x0 <- .05
# p0 <- .05
# BMR <- .05
# BMDL = c("all", "score", "delta")
# monotone = TRUE
# verbose = TRUE
# eps = 1e-06
# maxitr = 100
# boot <- 0
# bayes_boot <- 1000


mod1 <- benchmark_dose(y~s(x,bs='mpd'),data=dat,exposure = 'x',x0=.05)
mod2 <- benchmark_dose(y~s(x,bs='mpd')+x1+x2,data=dat,exposure = 'x',x0=.05)
mod3 <- tryCatch(benchmark_dose(y~s(x,bs='mpd')+x1+x2,data=dat,exposure = 'x1',x0=.05),error=function(e)e)

mod1v <- benchmark_dose(y~s(x,bs='mpd'),data=dat,exposure = 'x',x0=.05,verbose=TRUE)
mod2v <- benchmark_dose(y~s(x,bs='mpd')+x1+x2,data=dat,exposure = 'x',x0=.05,verbose=TRUE)
mod3v <- tryCatch(benchmark_dose(y~s(x,bs='mpd')+x1+x2,data=dat,exposure = 'x1',x0=.05),error=function(e)e,verbose=TRUE)

mod1g <- benchmark_dose(y~s(x,bs='bs'),data=dat,exposure = 'x',x0=.05,monotone = FALSE)
mod2g <- benchmark_dose(y~s(x,bs='bs')+x1+x2,data=dat,exposure = 'x',x0=.05,monotone = FALSE)
mod3g <- tryCatch(benchmark_dose(y~s(x,bs='bs')+x1+x2,data=dat,exposure = 'x1',x0=.05),error=function(e)e,monotone = FALSE)


# Only single BMDL calculation
mod1_score <- benchmark_dose(y~s(x,bs='mpd'),data=dat,exposure = 'x',x0=.05,BMDL='score')
mod1_delta <- benchmark_dose(y~s(x,bs='mpd'),data=dat,exposure = 'x',x0=.05,BMDL='delta')
mod1_none <- benchmark_dose(y~s(x,bs='mpd'),data=dat,exposure = 'x',x0=.05,BMDL='none')


# Bootstrapping
mod1_boot <- benchmark_dose(y~s(x,bs='mpd'),data=dat,exposure = 'x',x0=.05,boot=10)


# Bayesian bootstrapping
mod1_bayes <- benchmark_dose(y~s(x,bs='mpd'),data=dat,exposure = 'x',x0=.05,bayes_boot=100)
mod1_bothboot <- benchmark_dose(y~s(x,bs='mpd'),data=dat,exposure = 'x',x0=.05,boot=10,bayes_boot=100)




## TMB ##
set_parameters_mono <- function(scale=1,sigma=1,p0=GLOBAL_P0,BMR=.025,plot_it=FALSE,knots=10) {
  xmin <- 0
  xmax <- 1
  # True function f
  stopifnot(scale>0)
  f <- function(x) exp(-x*scale)
  finv <- function(x) -log(x)/scale

  x0 <- 0
  stopifnot(x0>=xmin)
  tau0 <- f(x0) + sigma * qnorm(p0) # Response level
  stopifnot(abs(pnorm(tau0,f(x0),sigma) - p0) < sqrt(.Machine$double.eps)) # Definition of tau

  A <- qnorm(p0+BMR) - qnorm(p0)
  # Benchmark dose
  xb <- finv(f(x0) - sigma*A)
  stopifnot(abs(pnorm(tau0,f(xb),sigma) - (p0+BMR)) < sqrt(.Machine$double.eps)) # Definition of benchmark dose

  if (plot_it) {
    curve(f,xmin,xmax,ylim=c(0,1))
    abline(v=xb,lty='dashed')
    cat("xb = ",xb,"; f'(xb)/sigma = ",numDeriv::grad(f,xb)/sigma,".\n",sep="")
  } else {
    list(
      xmin = xmin,xmax = xmax,
      sigma = sigma,
      p0 = p0,BMR = BMR,
      f = f,
      x0 = x0,
      tau0 = tau0,
      A = A,
      xb = xb,
      knots=knots
    )
  }
}

simulate_data_mono <- function(n,params) {
  xcov <- with(params,seq(xmin,xmax,length.out=n))
  with(params,data.frame(y = rnorm(n,f(xcov),sigma),x = xcov))
}

set.seed(472389)
params2 <- set_parameters_mono(scale=1,sigma=.5,p0=.01,BMR=.01,knots=10)
dat2 <- simulate_data_mono(500,params2)

mod1_tmb <- benchmark_dose_tmb(
  monosmooths = list(s(x,bs='bs')),
  smooths = NULL,
  linearterms = NULL,
  data = dat2,
  exposure = 'x',
  response = 'y',
  x0 = 0,
  p0 = .01,
  BMR = .01,
  verbose = TRUE,
  eps = 1e-06,
  maxitr = 10,
  bayes_boot = 1e03
)

# monosmooths <- list(s(x,bs='bs'))
# smooths <- NULL
# params2 <- set_parameters_mono(scale=1,sigma=.5,p0=.01,BMR=.01,knots=10)
# data <- dat2
# exposure <- 'x'
# response <- 'y'
# x0 <- 0
# p0 <- .01
# BMR <- .01
# verbose=FALSE
# eps=1e-06
# maxitr=10
# bayes_boot=1e03


## Multiple Smooth components ##

simulate_data_multi <- function(n,params) {
  xcov <- with(params,seq(xmin,xmax,length.out=n))
  z1 <- with(params,runif(n,xmin,xmax))
  z2 <- with(params,runif(n,xmin,xmax))


  fz1 <- function(z) sin(2*pi*z)
  fz2 <- function(z) cos(2*pi*z)
  # fz2 <- function(z) dnorm(2*z-.5)/dnorm(0)



  with(params,data.frame(y = rnorm(n,f(xcov)+fz1(z1)+fz2(z2),sigma),x = xcov,z1=z1,z2=z2))
}


set.seed(472398)
params3 <- set_parameters_mono(scale=5,sigma=.5,p0=.01,BMR=.01,knots=10)
dat3 <- simulate_data_multi(500,params3)

mod2_tmb <- benchmark_dose_tmb(
  monosmooths = list(s(x,bs='bs')),
  smooths = list(s(z1,bs='bs'),s(z2,bs='bs')),
  linearterms = NULL,
  data = dat3,
  exposure = 'x',
  response = 'y',
  x0 = 0,
  p0 = .01,
  BMR = .01,
  verbose = TRUE,
  eps = 1e-06,
  maxitr = 10,
  bayes_boot = 1e03
)


monosmooths <- list(s(x,bs='bs'))
smooths <- list(s(z1,bs='bs'),s(z2,bs='bs'))
linearterms <- NULL
data <- dat3
exposure <- 'x'
response <- 'y'
x0 <- 0
p0 <- .01
BMR <- .01
verbose=FALSE
eps=1e-06
maxitr=10
bayes_boot=1e03

library(scam)
testmod <- scam(y ~ s(x,bs='mpd')+s(z1,bs='bs')+s(z2,bs='bs'),data=dat3,family='gaussian')
summary(testmod)



## By variables ##

simulate_data_by <- function(n,params) {
  # Will return data of size 2n
  xcov <- rep(with(params,seq(xmin,xmax,length.out=n)),2)
  z <- with(params,runif(2*n,xmin,xmax))
  id <- rep(c(1,2),each=n)


  fz1 <- function(z) sin(2*pi*z)
  fz2 <- function(z) cos(2*pi*z)

  dat <- data.frame(
    x = xcov,
    z = z,
    id = factor(id)
  )
  mu <- params$f(xcov)
  mu[id==1] <- mu[id==1] + fz1(z[id==1])
  mu[id==2] <- mu[id==2] + fz2(z[id==2])
  dat$y = rnorm(2*n,mu,params$sigma)

  dat
}

set.seed(8732)
dat4 <- simulate_data_by(500,params3)


monosmooths <- list(s(x,bs='bs'))
smooths <- list(s(z,bs='bs',by=id))
linearterms <- NULL
data <- dat4
exposure <- 'x'
response <- 'y'
x0 <- 0
p0 <- .01
BMR <- .01
verbose=FALSE
eps=1e-06
maxitr=10
bayes_boot=1e03


mod3_tmb <- benchmark_dose_tmb(
  monosmooths = list(s(x,bs='bs')),
  smooths = list(s(z,bs='bs',by=id)),
  linearterms = NULL,
  data = dat4,
  exposure = 'x',
  response = 'y',
  x0 = 0,
  p0 = .01,
  BMR = .01,
  verbose = TRUE,
  eps = 1e-06,
  maxitr = 10,
  bayes_boot = 1e03
)
summary(mod3_tmb)













