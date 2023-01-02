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
mod1 <- benchmark_dose(y~s(x,bs='mpd'),data=dat,exposure = 'x',x0=.05)
mod2 <- benchmark_dose(y~s(x,bs='mpd')+x1+x2,data=dat,exposure = 'x',x0=.05)
mod3 <- tryCatch(benchmark_dose(y~s(x,bs='mpd')+x1+x2,data=dat,exposure = 'x1',x0=.05),error=function(e)e)

mod1v <- benchmark_dose(y~s(x,bs='mpd'),data=dat,exposure = 'x',x0=.05,verbose=TRUE)
mod2v <- benchmark_dose(y~s(x,bs='mpd')+x1+x2,data=dat,exposure = 'x',x0=.05,verbose=TRUE)
mod3v <- tryCatch(benchmark_dose(y~s(x,bs='mpd')+x1+x2,data=dat,exposure = 'x1',x0=.05),error=function(e)e,verbose=TRUE)

mod1g <- benchmark_dose(y~s(x,bs='bs'),data=dat,exposure = 'x',x0=.05,monotone = FALSE)
mod2g <- benchmark_dose(y~s(x,bs='bs')+x1+x2,data=dat,exposure = 'x',x0=.05,monotone = FALSE)
mod3g <- tryCatch(benchmark_dose(y~s(x,bs='bs')+x1+x2,data=dat,exposure = 'x1',x0=.05),error=function(e)e,monotone = FALSE)

