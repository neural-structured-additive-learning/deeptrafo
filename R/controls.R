#' Options for transformation models
#'
#' @param order_bsp The order of Bernstein polynomials in case \code{y_basis_fun}
#'     is a Bernstein polynomial defined by \code{eval_bsp} or (one less than)
#'     the number of classes of an ordinal outcome.
#' @param support A function returning a vector with two elements, namely
#'     the support for the basis of y.
#' @param penalize_bsp Scalar value > 0; controls amount of penalization of
#'     Bernstein polynomials.
#' @param order_bsp_penalty Integer; order of Bernstein polynomial penalty. 0
#'     results in a penalty based on integrated squared second order derivatives,
#'     values >= 1 in difference penalties.
#' @param tf_bsps Logical; whether to use a TensorFlow implementation of the
#'     Bernstein polynomial functions.
#' @param response_type Character; type of response can be continuous, ordered,
#'     survival, or count.
#' @param atm_toplayer Function; a function specifying the layer on top of ATM
#'     lags.
#' @param y_basis_fun Function; basis function for Y
#' @param y_basis_fun_lower Function; basis function for lower bound of interval
#'     censored response
#' @param y_basis_fun_prime Function; basis function derivative
#' @param basis Character or function; implemented options are
#'     \code{"bernstein"} (a Bernstein polynomial basis), \code{"ordered"}
#'     (for ordinal responses), or \code{"shiftscale"} for (log-) linear bases
#'
#' @return Returns a named \code{list} with all options, basis functions,
#'     support, and penalties.
#'
#' @export
#'
trafo_control <- function(order_bsp = 10L,
                          support = function(y) range(y),
                          y_basis_fun = NULL,
                          y_basis_fun_lower = NULL,
                          y_basis_fun_prime = NULL,
                          penalize_bsp = 0,
                          order_bsp_penalty = 2,
                          tf_bsps = FALSE,
                          response_type = c("continuous", "ordered", "survival", "count"),
                          atm_toplayer = function(x) {
                            layer_dense(x,
                              units = 1L,
                              name = "atm_toplayer",
                              use_bias = FALSE
                            )
                          },
                          basis = c("bernstein", "ordered", "shiftscale")) {
  if (!reticulate::py_available()) {
    message(
      "No Python Environemt available. Use `check_and_install()` ",
      "to install recommended environment."
    )
    invisible(return(NULL))
  }
  if (!reticulate::py_module_available("tensorflow")) {
    message("Tensorflow not available. Use `install_tensorflow()`.")
    invisible(return(NULL))
  }

  response_type <- match.arg(response_type)

  trafo <- if (is.function(basis)) {
    basis
  } else {
    basis <- match.arg(basis)

    switch(basis,
      "bernstein" = mono_trafo_multi,
      "ordered" = mono_trafo_multi,
      "shiftscale" = shift_scale_trafo_multi
    )
  }


  # define support (either function or dummy function outputting the supplied range)
  if (!is.function(support)) {
    supp <- function(x) support
  } else {
    supp <- support
  }

  if (response_type == "survival") {
    supp <- function(x) support(x[, 1])
  }

  # define bsp functions
  if (tf_bsps & is.null(y_basis_fun) & is.null(y_basis_fun_prime)) {
    if (is.function(supp)) {
      stop("Support must be given if TensorFlow Bernstein Basis Polynomials are used.")
    }

    eval_bsp <- eval_bsp_tf(order_bsp, supp)
    eval_bsp_prime <- eval_bsp_prime_tf(order_bsp, supp)
  }

  if (is.null(y_basis_fun)) {
    y_basis_fun <- switch(response_type,
      "continuous" = function(y, orderbsp = order_bsp, suppy = supp(y)) {
        eval_bsp(y, order = orderbsp, supp = suppy)
      },
      "survival" = function(y, orderbsp = order_bsp, suppy = supp(y)) {
        eval_bsp(y[, 1], order = orderbsp, supp = suppy)
      },
      "count" = function(y, orderbsp = order_bsp, suppy = supp(y)) {
        eval_bsp(y, order = orderbsp, supp = suppy)
      },
      "ordered" = function(y, orderbsp = order_bsp, suppy = suppy) {
        eval_ord_upr(y)
      }
    )
  }

  if (is.null(y_basis_fun_lower)) {
    y_basis_fun_lower <- switch(response_type,
      "continuous" = function(y, orderbsp = order_bsp, suppy = supp(y)) {
        ret <- eval_bsp(y, order = orderbsp, supp = suppy)
        ret[] <- 1e-8
        ret
      },
      "survival" = function(y, orderbsp = order_bsp, suppy = supp(y)) {
        ret <- eval_bsp(y[, 1], order = orderbsp, supp = suppy)
      },
      "count" = function(y, orderbsp = order_bsp, suppy = supp(y)) {
        eval_bsp(y - 1L, order = orderbsp, supp = suppy)
      },
      "ordered" = function(y, orderbsp = order_bsp, suppy = suppy) {
        eval_ord_lwr(y)
      }
    )
  }

  if (is.null(y_basis_fun_prime)) {
    y_basis_fun_prime <- switch(response_type,
      "continuous" = function(y, orderbsp = order_bsp,
                              suppy = supp(y)) { # / diff(supp(y))) {
        eval_bsp_prime(y, order = orderbsp, supp = suppy)
      },
      "survival" = function(y, orderbsp = order_bsp,
                            suppy = supp(y)) { # / diff(supp(y))) {
        eval_bsp_prime(y[, 1], order = orderbsp, supp = suppy)
      },
      "count" = function(y, orderbsp = order_bsp, suppy = supp(y)) {
        ret <- eval_bsp(y - 1L, order = orderbsp, supp = suppy)
        ret[] <- 1e-8
        ret
      },
      "ordered" = function(y, orderbsp = order_bsp, suppy = suppy) {
        ret <- eval_ord_lwr(y)
        ret[] <- 1e-8
        ret
      }
    )
  }


  if (response_type == "ordered") {
    penalize_bsp <- 0
    order_bsp_penalty <- 1
  }

  return(
    list(
      y_basis_fun = y_basis_fun,
      y_basis_fun_lower = y_basis_fun_lower,
      y_basis_fun_prime = y_basis_fun_prime,
      penalize_bsp = penalize_bsp,
      order_bsp_penalty = order_bsp_penalty,
      response_type = response_type,
      order_bsp = order_bsp,
      atm_toplayer = atm_toplayer,
      supp = supp,
      trafo = trafo
    )
  )
}
