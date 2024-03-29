% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{get_data}
\alias{get_data}
\title{Get data used to fit the dose-response model}
\usage{
get_data(object, ...)
}
\arguments{
\item{object}{Object of class \code{semibmd}
returned by \code{benchmark_dose_tmb}}

\item{...}{Not used}
}
\value{
A list containing elements \code{data}, the \code{data.frame} passed to \code{benchmark_dose_tmb};
\code{response}, the name of the response variable; \code{exposure}, the name of the exposure variable; and
\code{smooths}, a list with smoothing formulas used for the monotone and nonmonotone parts of the model.
}
\description{
Helper function to extract the data used to fit the dose-response model in an object of class \code{semibmd}
returned by \code{benchmark_dose_tmb}.
}
