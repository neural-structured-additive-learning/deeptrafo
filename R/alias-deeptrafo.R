
#' Deep conditional transformation models with alternative formula interface
#'
#' @param response Formula for the response; e.g. \code{~ y}
#' @param intercept Formula for the intercept function; e.g., \code{~ x},
#'     for which interacting bases with the response will be set up
#' @param shift Formula for the shift part of the model; e.g., \code{~ s(x)}
#' @param shared Formula for sharing weights between predictors in the intercept
#'     and shift part of the model
#' @inheritParams deeptrafo
#'
#' @return See return statement of \code{\link[deeptrafo]{deeptrafo}}
#'
#' @examples
#' df <- data.frame(y = rnorm(50), x = rnorm(50))
#' m <- dctm(response = ~ y, shift = ~ x, data = df)
#' fit(m, epochs = 2L)
#'
#' @export
#'
dctm <- function(
  response, intercept = NULL, shift = NULL,
  shared = NULL, data, lag_formula = NULL,
  response_type = get_response_type(data[[all.vars(response)[1]]]),
  order = get_order(response_type, data[[all.vars(response)[1]]]),
  addconst_interaction = 0, family = "logistic", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  fml <- forms2form(response, intercept, shift, shared)

  ret <- deeptrafo(formula = fml, data = data, lag_formula = lag_formula,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, family = family,
                   monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                   ... = ...)

  class(ret) <- c("dctm", class(ret))

  ret

}

#' Ordinal neural network transformation models
#'
#' @param response Formula for the response; e.g., \code{~ y}
#' @param intercept Formula for the intercept function; e.g., \code{~ x},
#'     for which interacting bases with the response will be set up
#' @param shift Formula for the shift part of the model; e.g., \code{~ s(x)}
#' @param shared Formula for sharing weights between predictors in the intercept
#'     and shift part of the model
#' @inheritParams deeptrafo
#'
#' @return See return statement of \code{\link[deeptrafo]{deeptrafo}}
#'
#' @references Kook, L. & Herzog, L., Hothorn, T., Dürr, O., & Sick, B. (2022).
#'     Deep and interpretable regression models for ordinal outcomes.
#'     Pattern Recognition, 122, 108263. DOI 10.1016/j.patcog.2021.108263
#'
#' @examples
#' df <- data.frame(y = ordered(sample.int(6, 50, TRUE)), x = rnorm(50))
#' m <- ontram(response = ~ y, shift = ~ x, data = df)
#' fit(m, epochs = 2L)
#'
#' @export
#'
ontram <- function(
  response, intercept = NULL, shift = NULL,
  shared = NULL, data, lag_formula = NULL,
  response_type = "ordered",
  order = get_order(response_type, data[[all.vars(response)[1]]]),
  addconst_interaction = 0, family = "logistic", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  stopifnot(is.ordered(data[[all.vars(response)[1]]][1]))

  ret <- dctm(response = response, intercept = intercept, shift = shift,
              shared = shared, data = data, lag_formula = lag_formula,
              response_type = response_type, order = order,
              addconst_interaction = addconst_interaction, family = family,
              monitor_metrics = monitor_metrics, trafo_options = trafo_options,
              ... = ...)

  class(ret) <- c("ontram", class(ret))

  ret

}

#' Deep continuous outcome logistic regression
#'
#' @inheritParams deeptrafo
#'
#' @return See return statement of \code{\link[deeptrafo]{deeptrafo}}
#'
#' @examples
#' df <- data.frame(y = rnorm(50), x = rnorm(50))
#' m <- ColrNN(y ~ x, data = df)
#' fit(m, epochs = 2L)
#'
#' @export
#'
ColrNN <- function(
  formula, data, lag_formula = NULL,
  response_type = get_response_type(data[[all.vars(formula)[1]]]),
  order = get_order(response_type, data[[all.vars(formula)[1]]]),
  addconst_interaction = 0, family = "logistic", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- deeptrafo(formula = formula, data = data, lag_formula = lag_formula,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, family = family,
                   monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                   ... = ...)

  class(ret) <- c("ColrNN", class(ret))

  ret

}

#' Lehmann-type neural network transformation models
#'
#' @inheritParams deeptrafo
#'
#' @return See return statement of \code{\link[deeptrafo]{deeptrafo}}
#'
#' @examples
#' df <- data.frame(y = rnorm(50), x = rnorm(50))
#' m <- LehmannNN(y ~ x, data = df)
#' fit(m, epochs = 2L)
#'
#' @export
#'
LehmanNN <- function(
  formula, data, lag_formula = NULL,
  response_type = get_response_type(data[[all.vars(formula)[1]]]),
  order = get_order(response_type, data[[all.vars(formula)[1]]]),
  addconst_interaction = 0, family = "gumbel", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- deeptrafo(formula = formula, data = data, lag_formula = lag_formula,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, family = family,
                   monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                   ... = ...)

  class(ret) <- c("LehmanNN", class(ret))

  ret

}

#' BoxCox-type neural network transformation models
#'
#' @inheritParams deeptrafo
#'
#' @return See return statement of \code{\link[deeptrafo]{deeptrafo}}
#'
#' @examples
#' df <- data.frame(y = rnorm(50), x = rnorm(50))
#' m <- BoxCoxNN(y ~ x, data = df)
#' fit(m, epochs = 2L)
#'
#' @export
#'
BoxCoxNN <- function(
  formula, data, lag_formula = NULL,
  response_type = get_response_type(data[[all.vars(formula)[1]]]),
  order = get_order(response_type, data[[all.vars(formula)[1]]]),
  addconst_interaction = 0, family = "normal", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- deeptrafo(formula = formula, data = data, lag_formula = lag_formula,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, family = family,
                   monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                   ... = ...)

  class(ret) <- c("BoxCoxNN", class(ret))

  ret

}

#' Deep (proportional odds) logistic regression
#'
#' @inheritParams deeptrafo
#'
#' @return See return statement of \code{\link[deeptrafo]{deeptrafo}}
#'
#' @examples
#' df <- data.frame(y = ordered(sample.int(5, 50, replace = TRUE)),
#'      x = rnorm(50))
#' m <- PolrNN(y ~ x, data = df)
#' fit(m, epochs = 2L)
#'
#' @export
#'
PolrNN <- function(
  formula, data, lag_formula = NULL,
  response_type = get_response_type(data[[all.vars(formula)[1]]]),
  order = get_order(response_type, data[[all.vars(formula)[1]]]),
  addconst_interaction = 0, family = "logistic", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  stopifnot(response_type == "ordered")

  ret <- deeptrafo(formula = formula, data = data, lag_formula = lag_formula,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, family = family,
                   monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                   ... = ...)

  class(ret) <- c("PolrNN", class(ret))

  ret

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
