% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/controls.R
\name{trafo_control}
\alias{trafo_control}
\title{Options for transformation models}
\usage{
trafo_control(
  order_bsp = 10L,
  support = function(y) range(y),
  y_basis_fun = NULL,
  y_basis_fun_lower = NULL,
  y_basis_fun_prime = NULL,
  penalize_bsp = 0,
  order_bsp_penalty = 2,
  tf_bsps = FALSE,
  response_type = c("continuous", "ordered", "survival", "count"),
  atm_toplayer = function(x) layer_dense(x, units = 1L, name = "atm_toplayer", use_bias =
    FALSE),
  basis = c("bernstein", "ordered", "shiftscale")
)
}
\arguments{
\item{order_bsp}{The order of Bernstein polynomials in case \code{y_basis_fun}
is a Bernstein polynomial defined by \code{eval_bsp} or (one less than)
the number of classes of an ordinal outcome.}

\item{support}{A function returning a vector with two elements, namely
the support for the basis of y.}

\item{y_basis_fun}{Function; basis function for Y}

\item{y_basis_fun_lower}{Function; basis function for lower bound of interval
censored response}

\item{y_basis_fun_prime}{Function; basis function derivative}

\item{penalize_bsp}{Scalar value > 0; controls amount of penalization of
Bernstein polynomials.}

\item{order_bsp_penalty}{Integer; order of Bernstein polynomial penalty. 0
results in a penalty based on integrated squared second order derivatives,
values >= 1 in difference penalties.}

\item{tf_bsps}{Logical; whether to use a TensorFlow implementation of the
Bernstein polynomial functions.}

\item{response_type}{Character; type of response can be continuous, ordered,
survival, or count.}

\item{atm_toplayer}{Function; a function specifying the layer on top of ATM
lags.}

\item{basis}{Character or function; implemented options are
\code{"bernstein"} (a Bernstein polynomial basis), \code{"ordered"}
(for ordinal responses), or \code{"shiftscale"} for (log-) linear bases}
}
\value{
Returns a named \code{list} with all options, basis functions,
    support, and penalties.
}
\description{
Options for transformation models
}
