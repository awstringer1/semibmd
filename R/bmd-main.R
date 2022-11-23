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
#'
#' @details This function solves \code{g(x)=0} for a solution satisfying \code{x>=min(bounds)},
#' \code{x<=max(bounds)}, and \code{|g(x)|<eps}.
#'
#' @return A numeric value representing the root.
#'
#' @examples
#' bounded_newton(function(x) x,c(-1,1))
#' bounded_newton(~.x,c(-1,1))
#' bounded_newton("force",c(-1,1))
#'
#' @export
bounded_newton <- function(g,bounds,eps = 1e-06) {
  stopifnot(length(bounds)==2)
  g <- rlang::as_function(g)
  gprime <- function(u) numDeriv::grad(g,u,method = "simple")
  # t <- 0
  xt <- mean(bounds) # Starting value
  while(abs(g(xt)) > eps) {
    # t <- t+1
    xt <- xt - g(xt)/gprime(xt)
    xt <- reflect(xt,bounds)
  }
  xt
}
