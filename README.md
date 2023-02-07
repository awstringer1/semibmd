Semi-parametric Benchmark Dosing with semibmd
================
Alex Stringer
2023-01-31

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
n <- 100
f <- function(x) 1/sqrt(x+.05)
xmin <- 0
xmax <- .2
xcov <- seq(xmin,xmax,length.out=n)
x1 <- runif(n,xmin,xmax)
x2 <- runif(n,xmin,xmax)
set.seed(43798)
dat <- data.frame(y = rnorm(n,f(xcov)+2*x1-x2,1),x = xcov,x1=x1,x2=x2)
head(dat)
```

    ##          y           x          x1         x2
    ## 1 5.353954 0.000000000 0.036759392 0.16945822
    ## 2 5.988536 0.002020202 0.002334835 0.11977962
    ## 3 5.048516 0.004040404 0.132281661 0.03496355
    ## 4 3.918865 0.006060606 0.036003663 0.09785578
    ## 5 2.107215 0.008080808 0.025234051 0.03527543
    ## 6 3.312403 0.010101010 0.116874519 0.13823451

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
    ## (Intercept)   4.3253     0.4442   9.738 5.56e-16 ***
    ## x1            0.2075     1.8659   0.111    0.912    
    ## x2           -1.4292     1.9617  -0.729    0.468    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##       edf Ref.df     F  p-value    
    ## s(x) 1.23  1.413 32.88 4.66e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.2855   Deviance explained = 30.9%
    ## GCV score = 1.2327  Scale est. = 1.1805    n = 100
    ## 
    ## ---
    ## Benchmark dose summary:
    ## ---
    ##      bmd   bmdl
    ## 1 0.0293 0.0184
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
    ## (Intercept)   2.6722     0.2915   9.167  9.3e-15 ***
    ## x1            0.1978     1.8658   0.106    0.916    
    ## x2           -1.4285     1.9623  -0.728    0.468    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##        edf Ref.df     F p-value    
    ## s(x) 1.186  1.348 32.47  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.285   Deviance explained = 30.8%
    ## GCV = 1.2328  Scale est. = 1.1812    n = 100
    ## ---
    ## Benchmark dose summary:
    ## ---
    ##      bmd   bmdl
    ## 1 0.0302 0.0202
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

    ## [1] 0.03018814 0.02023496

``` r
# All BMDLs:
get_all_bmdl(mod)
```

    ##      score      delta  bootstrap 
    ## 0.02023496 0.01583446 0.01197126

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
    ## 1.470831e-09

``` r
get_psixl(mod)
```

    ##               1
    ## 1 -8.403449e-08

``` r
# Approximation quantities: variance and derivative of U at BMD (used for delta method)
get_approximations(mod)
```

    ##          Vn         Upn 
    ##  0.00762825 11.92607537

``` r
# Computation times of each step
get_computation_times(mod)
```

    ##       model         bmd  bmdl_score  bmdl_delta   bmdl_boot       total 
    ## 0.012759924 0.009022951 0.017812014 0.007297993 8.015168905 8.062061787
