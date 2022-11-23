Semi-parametric Benchmark Dosing with semibmd
================
Alex Stringer
2022-11-23

# Install the `semibmd` package

The package is hosted at
(<https://github.com/awstringer1/semibmd>)\[<https://github.com/awstringer1/semibmd>\].
It’s currently a private repo so you can’t install using
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
f <- function(x) 1/sqrt(x)
xmin <- .05
xmax <- .2
xcov <- seq(xmin,xmax,length.out=n)
x1 <- runif(n,xmin,xmax)
x2 <- runif(n,xmin,xmax)
set.seed(43798)
dat <- data.frame(y = rnorm(n,f(xcov)+2*x1-x2,1),x = xcov,x1=x1,x2=x2)
head(dat)
```

    ##          y          x         x1         x2
    ## 1 5.500323 0.05000000 0.10184548 0.15326083
    ## 2 6.236828 0.05151515 0.08778823 0.06383406
    ## 3 4.947664 0.05303030 0.09106139 0.09415035
    ## 4 4.279017 0.05454545 0.17273698 0.06942820
    ## 5 2.273126 0.05606061 0.12235470 0.13770802
    ## 6 3.616912 0.05757576 0.18324410 0.05495769

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
                      x0=.05,
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
discussions.

The rest of the arguments are fairly self-explanatory. You have to tell
it which variable is the exposure variable (that you want the BMD(L)
for) using `exposure=`. You can look at `?benchmark_dose` for the full
documentation.

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
    ## (Intercept)   4.6199     0.5382   8.584 1.62e-13 ***
    ## x1           -1.5349     2.4776  -0.620    0.537    
    ## x2           -0.6121     2.4001  -0.255    0.799    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##      edf Ref.df     F  p-value    
    ## s(x)   1      1 40.27 5.08e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.2806   Deviance explained = 30.2%
    ## GCV score = 1.2043  Scale est. = 1.1561    n = 100
    ## 
    ## ---
    ## Benchmark dose summary:
    ## ---
    ##      bmd   bmdl
    ## 1 0.0744 0.0686
    ## ---

``` r
plot(mod)
```

![](README_files/figure-gfm/summary_plot-1.png)<!-- -->

For more detailed access, you can get the actual fitted model using
`get_model(mod)` and the estimated BMD(L) using `get_bmd(mod)`.
