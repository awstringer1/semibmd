Semi-parametric Benchmark Dosing with semibmd
================
Alex Stringer
2023-01-12

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

    ##          y           x         x1         x2
    ## 1 5.567913 0.000000000 0.14743287 0.17684554
    ## 2 6.073063 0.002020202 0.01475311 0.06008934
    ## 3 5.103770 0.004040404 0.17848509 0.07211622
    ## 4 3.949342 0.006060606 0.03116702 0.05770544
    ## 5 2.025384 0.008080808 0.02020776 0.10705371
    ## 6 3.121305 0.010101010 0.04657504 0.18873361

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
    ## (Intercept)  4.10929    0.54164   7.587 2.15e-11 ***
    ## x1           2.22784    1.89863   1.173    0.244    
    ## x2          -0.04648    2.04971  -0.023    0.982    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##        edf Ref.df     F  p-value    
    ## s(x) 1.459  1.758 12.71 4.35e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.302   Deviance explained = 32.6%
    ## GCV score = 1.2406  Scale est. = 1.1853    n = 100
    ## 
    ## ---
    ## Benchmark dose summary:
    ## ---
    ##      bmd   bmdl
    ## 1 0.0258 0.0143
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
    ## (Intercept)   2.3323     0.2614   8.923 3.15e-14 ***
    ## x1            2.2274     1.9011   1.172    0.244    
    ## x2           -0.0341     2.0511  -0.017    0.987    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##        edf Ref.df    F p-value    
    ## s(x) 1.351  1.621 26.4 2.5e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.301   Deviance explained = 32.5%
    ## GCV =  1.241  Scale est. = 1.187     n = 100
    ## ---
    ## Benchmark dose summary:
    ## ---
    ##      bmd   bmdl
    ## 1 0.0281 0.0176
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

    ## [1] 0.02814321 0.01761193

``` r
# All BMDLs:
get_all_bmdl(mod)
```

    ##      score      delta  bootstrap 
    ## 0.01761193 0.01199496 0.01134902

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

    ##           1 
    ## 1.93638e-08

``` r
get_psixl(mod)
```

    ##               1
    ## 1 -1.010736e-08

``` r
# Approximation quantities: variance and derivative of U at BMD (used for delta method)
get_approximations(mod)
```

    ##          Vn         Upn 
    ##  0.01095053 12.70106418
