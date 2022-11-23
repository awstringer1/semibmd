# Helper functions for semibmd

#' Reflection for bounded newton's method
#'
#' Implement the reflection of Coleman and Li (1994)
#'
#' @inheritParams bounded_newton
#' @param yt Point to reflect
#' @return Numeric value representing the reflected point
#'
#' @details Implements the formula \code{min(wt,2*(ub - lb) - wt) + lb} where
#' \code{wt = abs(yt - lb) mod (2*(ub-lb))}
#' and \code{lb = min(bounds),ub = max(bounds)}.
reflect <- function(yt,bounds) {
  lb <- min(bounds)
  ub <- max(bounds)
  wt <- abs(yt - lb) %% (2*(ub-lb))
  min(wt,2*(ub - lb) - wt) + lb
}
