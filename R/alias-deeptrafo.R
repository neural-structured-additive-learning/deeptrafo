
#' Deep conditional transformation models with alternative formula interface
#'
#' @param response formula for the response; e.g., \code{~ y}
#' @param intercept formula for the intercept function; e.g., \code{~ x},
#'     for which interacting bases with the response will be set up
#' @param shift formula for the shift part of the model; e.g., \code{~ s(x)}
#' @param shared formula for sharing weights between predictors in the intercept
#'     and shift part of the model
#' @inheritParams deeptrafo
#'
#' @return See return statement of \code{\link[deeptrafo]{deeptrafo}}
#'
#' @export
#'
dctm <- function(
  response, intercept = NULL, shift = NULL,
  shared = NULL, data, lag_formula = NULL,
  response_type = get_response_type(data[[all.vars(response)[1]]]),
  order = get_order(response_type, data[[all.vars(fml)[1]]]),
  addconst_interaction = 0, family = "logistic", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  fml <- forms2form(response, intercept, shift, shared)

  deeptrafo(formula = fml, data = data, lag_formula = lag_formula,
            response_type = response_type, order = order,
            addconst_interaction = addconst_interaction, family = family,
            monitor_metrics = monitor_metrics, trafo_options = trafo_options,
            ... = ...)

}

# Helpers -----------------------------------------------------------------

formula_parts <- function(x) {
  if (!is.null(x))
    gsub("~", "", Reduce(paste, deparse(x)))
}

forms2form <- function(response, intercept = NULL, shift = NULL, shared = NULL) {

  stopifnot(!is.null(response))

  rsp <- formula_parts(response)
  int <- formula_parts(intercept)
  shi <- formula_parts(shift)
  sha <- formula_parts(shared)

  if (is.null(int))
    lhs <- rsp
  else
    lhs <- paste0(rsp, "|", int)

  if (is.null(shi))
    shi <- 1

  if (is.null(sha))
    rhs <- shi
  else
    rhs <- paste0(shi, "|", sha)

  as.formula(paste0(lhs, "~", rhs))

}
