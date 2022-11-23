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
