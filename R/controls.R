#' Options for transformation models
#'
#' @param order_bsp the order of Bernstein polynomials in case \code{y_basis_fun}
#'     is a Bernstein polynomial defined by \code{eval_bsp}
#' @param support a function yielding a vector with two elements, namely
#'     the support for the basis of y;
#' @param y_basis_fun,y_basis_fun_prime basis functions for y transformation
#' @param penalize_bsp scalar value > 0; amount of penalization of Bernstein polynomials
#' @param order_bsp_penalty integer; order of Bernstein polynomial penalty. 0 results in a
#'     penalty based on integrated squared second order derivatives, values >= 1 in difference
#'     penalties
#' @param tf_bsps logical; whether to use a TensorFlow implementation of the Bernstein polynomial
#'     functions
#' @export
#'
trafo_control <- function(order_bsp = 10L,
                          support = function(y) range(y),
                          y_basis_fun = NULL,
                          y_basis_fun_prime = NULL,
                          penalize_bsp = 0,
                          order_bsp_penalty = 2,
                          tf_bsps = FALSE,
                          ordered = FALSE) {
  # define support (either function or dummy function outputting the supplied range)
  if (!is.function(support)) {

    supp <- function(x) support

  } else {

    supp <- support

  }

  # define bsp functions
  if (tf_bsps & is.null(y_basis_fun) & is.null(y_basis_fun_prime)) {

    if (is.function(supp))
      stop("Support must be given if TensorFlow Bernstein Basis Polynomials are used.")

    eval_bsp <- eval_bsp_tf(order_bsp, supp)
    eval_bsp_prime <- eval_bsp_prime_tf(order_bsp, supp)

  }

  if (is.null(y_basis_fun)) {

    y_basis_fun <- function(y, orderbsp = order_bsp, suppy = supp(y)) {
      eval_bsp(y, order = orderbsp, supp = suppy)
    }

  }

  if (is.null(y_basis_fun_prime)) {

    y_basis_fun_prime <- function(y, orderbsp = order_bsp,
                                  suppy = supp(y) / diff(supp(y))) {
      eval_bsp_prime(y, order = orderbsp, supp = suppy)
    }

  }

  if (ordered) {
    y_basis_fun <- function(y, orderbsp = order_bsp, suppy = suppy) {
      eval_ord_upr(y)
    }
    y_basis_fun_prime <- function(y, orderbsp = order_bsp, suppy = suppy) {
      eval_ord_lwr(y)
    }
    penalize_bsp <- 0
    order_bsp_penalty <- 1
  }

  return(
    list(y_basis_fun = y_basis_fun,
         y_basis_fun_prime = y_basis_fun_prime,
         penalize_bsp = penalize_bsp,
         order_bsp_penalty = order_bsp_penalty,
         ordered = ordered)
  )

}
