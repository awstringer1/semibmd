Semi-parametric Benchmark Dosing with semibmd
================
Alex Stringer
2023-03-01

# Install the `semibmd` package

The package is hosted at <https://github.com/awstringer1/semibmd>. It’s
currently a private repo so you can’t install using
e.g. `remotes::install_github`. Instead, clone the repository to a local
directory, say `~/work/projects/benchmark-dose/code/semibmd` for
example, and then in `R` do:

``` r
install.packages(pkgs='~/work/projects/benchmark-dose/code/semibmd',repos=NULL,type='source')
library(semibmd)
```

# Simulate some data

Simulate some data to which the dose-response model is to be fit:

``` r
n <- 1000
f <- function(x) exp(-5*x)
xmin <- 0
xmax <- 1
xcov <- seq(xmin,xmax,length.out=n)
x1 <- runif(n,xmin,xmax)
x2 <- runif(n,xmin,xmax)
set.seed(43798)
dat <- data.frame(y = rnorm(n,f(xcov)+2*x1-x2,1),x = xcov,x1=x1,x2=x2)
head(dat)
```

    ##            y           x        x1        x2
    ## 1  1.9994748 0.000000000 0.2244994 0.4272809
    ## 2  3.5659719 0.001001001 0.6792357 0.5067146
    ## 3  2.0950913 0.002002002 0.5415578 0.4952741
    ## 4  2.0352977 0.003003003 0.9846344 0.6402956
    ## 5 -0.2218419 0.004004004 0.6608238 0.4663065
    ## 6  0.5493814 0.005005005 0.6450404 0.8538231

The exposure variable is `x` and there are two additional variables `x1`
and `x2` that can be included in the model.

# Benchmark dosing

## Fit the model

Fit the dose-response model and calculate the benchmark dose
information:

``` r
mod <- benchmark_dose(y~s(x,bs='mpd')+x1+x2,
                      data=dat,
                      exposure = 'x',
                      x0=0,
                      p0=.05,
                      BMR=.05,
                      monotone = TRUE)
```

The first argument to `benchmark_dose` is a formula compatible with
`scam::scam` (if `monotone=TRUE`) or with `mgcv::gam` (if
`monotone=FALSE`). The `scam` package accepts formulas just like `gam`,
with the addition that terms you want to be monotone should have their
bases specified manually to be one of the special monotone spline basis.
See `?scam::scam` for a list. Here I chose `bs='mpd'` for “monotone
p-spline, decreasing”, to enforce a monotone-decreasing estimated curve.

If you choose `monotone = TRUE` then `scam` is used; I added the option
to choose `monotone=FALSE` and fit an ordinary `gam`, based on our
discussions. However, I have not yet tested this feature.

The rest of the arguments are fairly self-explanatory. You have to tell
it which variable is the exposure variable (that you want the BMD(L)
for) using `exposure=...`. You can look at `?benchmark_dose` for the
full documentation.

## Get the summaries

Summaries and plots are obtained the usual way:

``` r
summary(mod)
```

    ## ---
    ## Dose-response model summary:
    ## ---
    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## y ~ s(x, bs = "mpd") + x1 + x2
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   0.8613     0.3026   2.846  0.00451 ** 
    ## x1            1.9638     0.1081  18.174  < 2e-16 ***
    ## x2           -0.9853     0.1087  -9.063  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##        edf Ref.df     F  p-value    
    ## s(x) 2.006  2.455 14.83 7.62e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.3108   Deviance explained = 31.4%
    ## GCV score = 0.99827  Scale est. = 0.99328   n = 1000
    ## 
    ## ---
    ## Benchmark dose summary:
    ## ---
    ##      bmd   bmdl
    ## 1 0.2588 0.1149
    ## ---

``` r
plot(mod)
```

![](README_files/figure-gfm/summary_plot-1.png)<!-- -->

The `summary` method just appends the estimated BMD(L) onto the summary
from the dose-response model. The `plot` method just adds vertical lines
to the plotted `scam/gam` at the BMD and BMDL.

For more detailed access, you can get the actual fitted model using
`get_model(mod)` and the estimated BMD(L) using `get_bmd(mod)`.

## Regular non-monotone GAM

We can also fit a regular `gam` by setting `montone = FALSE`. In this
case, you change the basis in the formula to be something compatible
with `mgcv`:

``` r
mod <- benchmark_dose(y~s(x,bs='bs')+x1+x2, # B-spline
                      data=dat,
                      exposure = 'x',
                      x0=0,
                      p0=.05,
                      BMR=.05,
                      monotone = FALSE)
```

Summaries and plots are obtained the same way:

``` r
summary(mod)
```

    ## ---
    ## Dose-response model summary:
    ## ---
    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## y ~ s(x, bs = "bs") + x1 + x2
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.16284    0.08297   1.963     0.05 *  
    ## x1           1.96751    0.10818  18.188   <2e-16 ***
    ## x2          -0.98329    0.10869  -9.046   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##        edf Ref.df     F  p-value    
    ## s(x) 2.196  2.731 12.28 1.29e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.311   Deviance explained = 31.4%
    ## GCV = 0.99835  Scale est. = 0.99317   n = 1000
    ## ---
    ## Benchmark dose summary:
    ## ---
    ##      bmd   bmdl
    ## 1 0.3158 0.1758
    ## ---

``` r
plot(mod)
```

![](README_files/figure-gfm/summary_plot_gam-1.png)<!-- -->

# Alternative BMDL Calculation

I have added functionality for computing BMDLs using the delta method
and bootstrapping. The argument `BMDL` in the call to `benchmark_dose`
controls which BMDL(s) is/are calculated; current options are `all`
(default), `score` and `delta`. Note that the contents of the `bmdl`
slot in the output object is never anything other than the `score` BMDL;
the `delta` BMDL is stored elsewhere and accessed using a getter
function, see below.

The argument `boot` controls the number of (parametric) bootstrap
iterations, with default `0` meaning don’t bootstrap. Like the `delta`
BMDL, a bootstrapped BMDL will be stored elsewhere in the output object
and accessed using a getter.

Here is how to fit it:

``` r
mod <- benchmark_dose(y~s(x,bs='bs')+x1+x2,
                      data=dat,
                      exposure = 'x',
                      x0=0,
                      p0=.05,
                      BMR=.05,
                      monotone = TRUE,
                      BMDL = 'all',
                      boot = 200
)
```

The `summary` method still returns just the `score` BMDL. Here are the
various getter functions to get everything:

``` r
# BMD and score BMDL:
get_bmd(mod)
```

    ## [1] 0.3158169 0.1757731

``` r
# All BMDLs:
get_all_bmdl(mod)
```

    ##      score      delta  bootstrap 
    ## 0.17577313 0.09277874 0.18861372

# Errors and diagnostics

The `benchmark_dose` function does error checking/handling at each step
of the procedure and attempts to collate and report the errors in a safe
manner. If you get an actual exception thrown when using it, let me know
what it is and I’ll try and add it. You can see what errors were thrown
with the following:

``` r
get_errors(mod)
```

    ## [1] FALSE

This should return `FALSE` if no errors were thrown.

You can see diagnostics and approximation quantities too:

``` r
# Diagnostics: U at estimated BMD and Psi at estimated score BMDL.
# Both should be zero:
get_uxb(mod)
```

    ##            1 
    ## 1.799314e-09

``` r
get_psixl(mod)
```

    ##               1
    ## 1 -4.356204e-10

``` r
# Approximation quantities: variance and derivative of U at BMD (used for delta method)
get_approximations(mod)
```

    ##         Vn        Upn 
    ## 0.01040667 0.89644726

``` r
# Computation times of each step
get_computation_times(mod)
```

    ##        model          bmd   bmdl_score   bmdl_delta    bmdl_boot        total 
    ##  0.032787085  0.009737968  0.014474869  0.005775928 10.452320099 10.515095949

# New Smoothing Method using TMB

The method also is implemented from scratch using the Template Model
Builder (TMB; <https://kaskr.github.io/adcomp/>). The `mgcv` package is
still used to set up the bases and penalties, but then all the
computations are done using custom code. The following methodological
differences are present over `scam`:

- The use of Laplace-approximate marginal likelihood to estimate
  smoothing parameters instead of generalized cross-validation,
- The use of full integrated second derivative penalties for the
  monotone smooths (and other smooths) instead of difference-based
  P-spline penalties, and
- The direct application of linear sum-to-zero constraints on the
  monontone smooth, as opposed to constraint-absorbing
  reparameterizations. The latter are still used for the non-monotone
  smooth terms.

Simulations in paper show substantially improved performance of doing
this way compared to `scam`.

Here are two examples. Currently the method supports (a) a single
monotone smooth, plus (b) any number of non-monotone smooths, and (c) no
linear terms (although these could probably be provided as part of the
`smooths` argument, I have not tested this).

## A single monotone dose response curve only

Here is an example of a model with only a single monotone smooth and
nothing else. The interface is a bit different and is documented.
**Note**: important, you use the regular `bs='bs'` basis now for the
monotone dose-response smooth, **NOT** the `scam`-type `bs='mpd'`. The
monotonicity is handled internally and is controlled by you passing this
formula to the `monosmooths` argument of `benchmark_dose_tmb`.

``` r
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
library(mgcv) # Need this for s()
params2 <- set_parameters_mono(scale=1,sigma=.5,p0=.01,BMR=.01,knots=10)
dat2 <- simulate_data_mono(1000,params2)

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
  verbose = FALSE,
  eps = 1e-06,
  maxitr = 10,
  bayes_boot = 1e03
)
class(mod1_tmb)
```

    ## [1] "semibmd"

The `benchmark_dose_tmb` function is documented and does the work here.
Because the `benchmark_dose_tmb` function returns an object of class
`semibmd`, the same `plot`, `summary`, and other methods above work for
it:

``` r
summary(mod1_tmb)
```

    ## ---
    ## Dose-response model summary:
    ## ---
    ## [1] "Monotone Additive Model fit using Laplace-Approximate Marginal Likelihood in TMB"
    ## ---
    ## Benchmark dose summary:
    ## ---
    ##      bmd   bmdl
    ## 1 0.0541 0.0233
    ## ---

``` r
plot(mod1_tmb)
```

![](README_files/figure-gfm/monosummary-1.png)<!-- -->

The output of the plot method is expanded a bit: I include a spaghetti
plot of posterior samples of the fitted curve to enhance it a bit.

We can compare to the same model fit to the same data the old way:

``` r
mod1_scam <- benchmark_dose(
  y ~ s(x,bs='mpd'),
  data = dat2,
  exposure = 'x',
  x0 = 0,
  p0 = .01,
  BMR = .01,
  verbose = FALSE,
  eps = 1e-06,
  maxitr = 10
)
summary(mod1_scam)
```

    ## ---
    ## Dose-response model summary:
    ## ---
    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## y ~ s(x, bs = "mpd")
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   1.4326     0.2865       5 6.78e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##       edf Ref.df     F p-value    
    ## s(x) 3.11  3.681 45.17  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.1383   Deviance explained = 14.1%
    ## GCV score = 0.24044  Scale est. = 0.23945   n = 1000
    ## 
    ## ---
    ## Benchmark dose summary:
    ## ---
    ##      bmd   bmdl
    ## 1 0.0683 0.0353
    ## ---

``` r
plot(mod1_scam)
```

![](README_files/figure-gfm/monosimscam-1.png)<!-- -->

We can the the BMD/L are a bit different:

``` r
get_bmd(mod1_tmb)
```

    ## [1] 0.05413111 0.02328782

``` r
get_bmd(mod1_scam)
```

    ## [1] 0.06834004 0.03527283

``` r
get_all_bmdl(mod1_tmb)
```

    ##       score       delta  bmdl_bayes 
    ##  0.02328782 -0.01016350  0.02155320

``` r
get_all_bmdl(mod1_scam)
```

    ##      score      delta 
    ## 0.03527283 0.01106405

Again, simulations show the `TMB` approach is apparently better. Note
that in the `plot` method for `scam` the curve is shown **without** the
intercept, where for `TMB` it is shown **with** the intercept. This can
be changed if necessary and is pretty arbitrary. The sum of the curve at
the exposure values equals $0$ in what `scam` is plotting and equals the
estimated intercept in what `TMB` is plotting; again, pretty arbitrary
choice and also the BMD doesn’t depend on this.

## A monotone dose response curve and two arbitrary smooths

Here is an example with a monotone dose-response smooth and two other
smooths. **Note** the use of the `bs='bs'` for the monotone
dose-response curve, **not** the `scam`-type `bs='mpd'`.

``` r
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

summary(mod2_tmb)
```

    ## ---
    ## Dose-response model summary:
    ## ---
    ## [1] "Monotone Additive Model fit using Laplace-Approximate Marginal Likelihood in TMB"
    ## ---
    ## Benchmark dose summary:
    ## ---
    ##     bmd   bmdl
    ## 1 0.026 0.0142
    ## ---

``` r
params3$xb # True BMD
```

    ## [1] 0.02930584

The `plot` method works a little bit different: since there is more than
one curve to plot, it will use the same `Hit <Return> to see next plot`
thing that other packages do.

``` r
plot(mod2_tmb)
```

![](README_files/figure-gfm/multismoothplot-1.png)<!-- -->![](README_files/figure-gfm/multismoothplot-2.png)<!-- -->![](README_files/figure-gfm/multismoothplot-3.png)<!-- -->

You will probably want finer control over the plot data, to make them
look nice. You can call `plot` with `plot=FALSE` (sorry for the stupid
naming) to get the data for the plot, and then make the plot(s) yourself
so you can add titles and so on:

``` r
plotdata <- plot(mod2_tmb,plot=FALSE)
length(plotdata)
```

    ## [1] 3

``` r
Map(names,plotdata)
```

    ## $mono
    ## [1] "x"        "fitted"   "estimate" "lower"    "upper"   
    ## 
    ## $smooth_1
    ## [1] "x"        "fitted"   "estimate" "lower"    "upper"   
    ## 
    ## $smooth_2
    ## [1] "x"        "fitted"   "estimate" "lower"    "upper"

``` r
# Plot the monontone smooth
with(plotdata[[1]],plot(x,fitted[ ,1],type='l',col=scales::alpha('lightgrey',0.2),xlim = range(x),ylim = c(min(lower),max(upper)),main = "My custom plot title",xlab = "Exposure",ylab = "Response"))
M <- 100 # Or whatever
for (i in 2:M)
  with(plotdata[[1]],graphics::lines(x,fitted[ ,i],col=scales::alpha('lightgrey',0.2)))

with(plotdata[[1]],graphics::lines(x,estimate))
with(plotdata[[1]],graphics::lines(x,lower,lty='dashed'))
with(plotdata[[1]],graphics::lines(x,upper,lty='dashed'))
with(plotdata[[1]],graphics::abline(v = get_bmd(mod2_tmb)[1],lty='longdash'))
with(plotdata[[1]],graphics::abline(v = get_bmd(mod2_tmb)[2],lty='dotdash'))
```

![](README_files/figure-gfm/multismoothplotmanual-1.png)<!-- -->

The plot data for the other smooths is identical.

## By Variables

`by` variables in smooths are handled by `mgcv` internally, and the
necessary modifications have been made to `semibmd::benchmark_dose_tmb`
to handle them. Here is an example.

``` r
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
dat4 <- simulate_data_by(500,params3) # 2 groups, 500 per group
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
```

    ## ---
    ## Dose-response model summary:
    ## ---
    ## [1] "Monotone Additive Model fit using Laplace-Approximate Marginal Likelihood in TMB"
    ## ---
    ## Benchmark dose summary:
    ## ---
    ##      bmd   bmdl
    ## 1 0.0229 0.0146
    ## ---

``` r
params3$xb
```

    ## [1] 0.02930584

The `plot` method will plot the smooths for each factor level completely
separately. I can change this, I just did it this way because I didn’t
have a better idea.

``` r
plot(mod3_tmb)
```

![](README_files/figure-gfm/plot-by-1.png)<!-- -->![](README_files/figure-gfm/plot-by-2.png)<!-- -->![](README_files/figure-gfm/plot-by-3.png)<!-- -->
