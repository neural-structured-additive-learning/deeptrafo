
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
#' if (reticulate::py_module_available("tensorflow") &
#'     reticulate::py_module_available("keras") &
#'     reticulate::py_module_available("tensorflow_probability")) {
#'    m <- dctm(response = ~ y, shift = ~ 0 + x, data = df)
#'    coef(m)
#' }
#'
#' @export
#'
dctm <- function(
  response, intercept = NULL, shift = NULL,
  shared = NULL, data,
  response_type = get_response_type(data[[all.vars(response)[1]]]),
  order = get_order(response_type, data[[all.vars(response)[1]]]),
  addconst_interaction = 0, latent_distr = "logistic", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  fml <- forms2form(response, intercept, shift, shared)

  ret <- deeptrafo(formula = fml, data = data,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, latent_distr = latent_distr,
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
#' if (reticulate::py_module_available("tensorflow") &
#'     reticulate::py_module_available("keras") &
#'     reticulate::py_module_available("tensorflow_probability")) {
#'     m <- ontram(response = ~ y, shift = ~ x, data = df)
#'    coef(m)
#' }
#'
#' @export
#'
ontram <- function(
  response, intercept = NULL, shift = NULL,
  shared = NULL, data,
  response_type = "ordered",
  order = get_order(response_type, data[[all.vars(response)[1]]]),
  addconst_interaction = 0, latent_distr = "logistic", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  stopifnot(is.ordered(data[[all.vars(response)[1]]][1]))

  ret <- dctm(response = response, intercept = intercept, shift = shift,
              shared = shared, data = data,
              response_type = response_type, order = order,
              addconst_interaction = addconst_interaction, latent_distr = latent_distr,
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
#' if (reticulate::py_module_available("tensorflow") &
#'     reticulate::py_module_available("keras") &
#'     reticulate::py_module_available("tensorflow_probability")) {
#'     m <- ColrNN(y ~ x, data = df)
#'     coef(m)
#' }
#'
#' @export
#'
ColrNN <- function(
  formula, data,
  response_type = get_response_type(data[[all.vars(formula)[1]]]),
  order = get_order(response_type, data[[all.vars(formula)[1]]]),
  addconst_interaction = 0, latent_distr = "logistic", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- deeptrafo(formula = formula, data = data,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                   monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                   ... = ...)

  class(ret) <- c("ColrNN", class(ret))

  ret

}

#' Cox proportional hazards type neural network transformation models
#'
#' @inheritParams deeptrafo
#'
#' @return See return statement of \code{\link[deeptrafo]{deeptrafo}}
#'
#' @examples
#' df <- data.frame(y = rnorm(50), x = rnorm(50))
#' if (reticulate::py_module_available("tensorflow") &
#'     reticulate::py_module_available("keras") &
#'     reticulate::py_module_available("tensorflow_probability")) {
#'     m <- CoxphNN(y ~ x, data = df)
#'     coef(m)
#' }
#'
#' @export
#'
CoxphNN <- function(
  formula, data,
  response_type = get_response_type(data[[all.vars(formula)[1]]]),
  order = get_order(response_type, data[[all.vars(formula)[1]]]),
  addconst_interaction = 0, latent_distr = "gompertz", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- deeptrafo(formula = formula, data = data,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                   monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                   ... = ...)

  class(ret) <- c("CoxphNN", class(ret))

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
#' if (reticulate::py_module_available("tensorflow") &
#'     reticulate::py_module_available("keras") &
#'     reticulate::py_module_available("tensorflow_probability")) {
#'     m <- LehmanNN(y ~ 0 + x, data = df)
#'     coef(m)
#' }
#'
#' @export
#'
LehmanNN <- function(
  formula, data,
  response_type = get_response_type(data[[all.vars(formula)[1]]]),
  order = get_order(response_type, data[[all.vars(formula)[1]]]),
  addconst_interaction = 0, latent_distr = "gumbel", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- deeptrafo(formula = formula, data = data,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, latent_distr = latent_distr,
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
#' if (reticulate::py_module_available("tensorflow") &
#'     reticulate::py_module_available("keras") &
#'     reticulate::py_module_available("tensorflow_probability")) {
#'     m <- BoxCoxNN(y ~ x, data = df)
#'     coef(m)
#' }
#'
#' @export
#'
BoxCoxNN <- function(
  formula, data,
  response_type = get_response_type(data[[all.vars(formula)[1]]]),
  order = get_order(response_type, data[[all.vars(formula)[1]]]),
  addconst_interaction = 0, latent_distr = "normal", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- deeptrafo(formula = formula, data = data,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, latent_distr = latent_distr,
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
#' if (reticulate::py_module_available("tensorflow") &
#'     reticulate::py_module_available("keras") &
#'     reticulate::py_module_available("tensorflow_probability")) {
#'     m <- PolrNN(y ~ x, data = df)
#'     coef(m)
#' }
#'
#' @export
#'
PolrNN <- function(
  formula, data,
  response_type = get_response_type(data[[all.vars(formula)[1]]]),
  order = get_order(response_type, data[[all.vars(formula)[1]]]),
  addconst_interaction = 0, latent_distr = "logistic", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order, response_type = response_type),
  ...
) {

  stopifnot(response_type == "ordered")

  ret <- deeptrafo(formula = formula, data = data,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                   monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                   ... = ...)

  class(ret) <- c("PolrNN", class(ret))

  ret

}

#' Deep normal linear regression
#'
#' @inheritParams deeptrafo
#'
#' @return See return statement of \code{\link[deeptrafo]{deeptrafo}}
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = 10 + rnorm(50), x = rnorm(50))
#' if (reticulate::py_module_available("tensorflow") &
#'     reticulate::py_module_available("keras") &
#'     reticulate::py_module_available("tensorflow_probability")) {
#'     m <- LmNN(y ~ 0 + x, data = df)
#' \donttest{
#'     optimizer <- optimizer_adam(learning_rate = 0.01, decay = 4e-4)
#'     m <- LmNN(y ~ 0 + x, data = df, optimizer = optimizer)
#'     library(tram)
#'     fit(m, epochs = 900L, validation_split = 0)
#'     logLik(mm <- Lm(y ~ x, data = df)); logLik(m)
#'     coef(mm, with_baseline = TRUE); unlist(c(coef(m, which = "interacting"),
#'                                              coef(m, which = "shifting")))
#' }
#' }
#'
#' @export
#'
LmNN <- function(
  formula, data,
  response_type = get_response_type(data[[all.vars(formula)[1]]]),
  order = get_order(response_type, data[[all.vars(formula)[1]]]),
  addconst_interaction = 0, latent_distr = "normal", monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = 1L,
                                response_type = response_type,
                                y_basis_fun = eval_lin,
                                y_basis_fun_lower = .empty_fun(eval_lin),
                                y_basis_fun_prime = eval_lin_prime,
                                basis = "shiftscale"),
  ...
) {

  stopifnot(response_type == "continuous")

  ret <- deeptrafo(formula = formula, data = data,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                   monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                   ... = ...)

  class(ret) <- c("LmNN", class(ret))

  ret

}

#' Deep parametric survival regression
#'
#' @inheritParams deeptrafo
#'
#' @return See return statement of \code{\link[deeptrafo]{deeptrafo}}
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = abs(1 + rnorm(50)), x = rnorm(50))
#' if (reticulate::py_module_available("tensorflow") &
#'     reticulate::py_module_available("keras") &
#'     reticulate::py_module_available("tensorflow_probability")) {
#'     m <- SurvregNN(y ~ 0 + x, data = df)
#' \donttest{
#'     optimizer <- optimizer_adam(learning_rate = 0.01, decay = 4e-4)
#'     m <- SurvregNN(y ~ 0 + x, data = df, optimizer = optimizer)
#'     library(tram)
#'     fit(m, epochs = 500L, validation_split = 0)
#'     logLik(mm <- Survreg(y ~ x, data = df, dist = "loglogistic")); logLik(m)
#'     coef(mm, with_baseline = TRUE); unlist(c(coef(m, which = "interacting"),
#'                                              coef(m, which = "shifting")))
#' }
#' }
#'
#' @export
#'
SurvregNN <- function(
  formula, data,
  response_type = get_response_type(data[[all.vars(formula)[1]]]),
  order = get_order(response_type, data[[all.vars(formula)[1]]]),
  addconst_interaction = 0, latent_distr = "gompertz", monitor_metrics = NULL,
  trafo_options = NULL,
  ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  if (response_type == "survival") {
    ybf <- function(y) eval_loglin(y[, 1])
    ybfl <- function(y) eval_loglin(y[, 1])
    ybfp <- function(y) eval_loglin_prime(y[, 1])
  } else {
    ybf <- eval_loglin
    ybfl <- .empty_fun(ybf)
    ybfp <- eval_loglin_prime
  }

  trafo_options <- trafo_control(
    order_bsp = 1L,
    response_type = response_type,
    y_basis_fun = ybf,
    y_basis_fun_lower = ybfl,
    y_basis_fun_prime = ybfp,
    basis = "shiftscale"
  )

  ret <- deeptrafo(formula = formula, data = data,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                   monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                   ... = ...)

  class(ret) <- c("SurvregNN", class(ret))

  ret

}

#' Deep distribution-free count regression
#'
#' @inheritParams deeptrafo
#'
#' @return See return statement of \code{\link[deeptrafo]{deeptrafo}}
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(y = as.integer(abs(1 + rnorm(50, sd = 10))), x = rnorm(50))
#' if (reticulate::py_module_available("tensorflow") &
#'     reticulate::py_module_available("keras") &
#'     reticulate::py_module_available("tensorflow_probability")) {
#'     m <- cotramNN(y ~ 0 + x, data = df, order = 6)
#' \donttest{
#'     optimizer <- optimizer_adam(learning_rate = 0.1, decay = 4e-4)
#'     m <- cotramNN(y ~ 0 + x, data = df, optimizer = optimizer, order = 6)
#'     library(cotram)
#'     fit(m, epochs = 800L, validation_split = 0)
#'     logLik(mm <- cotram(y ~ x, data = df, method = "logit")); logLik(m)
#'     coef(mm, with_baseline = TRUE); unlist(c(coef(m, which = "interacting"),
#'                                              coef(m, which = "shifting")))
#' }
#' }
#'
#' @export
#'
cotramNN <- function(
  formula, data,
  response_type = get_response_type(data[[all.vars(formula)[1]]]),
  order = get_order(response_type, data[[all.vars(formula)[1]]]),
  addconst_interaction = 0, latent_distr = "logistic", monitor_metrics = NULL,
  ...
) {

  stopifnot(response_type == "count")

  tsupp <- range(data[[all.vars(formula)[1]]])
  trafo_options <- trafo_control(
    order_bsp = order,
    response_type = response_type,
    y_basis_fun = .get_eval_cotram(order, tsupp),
    y_basis_fun_lower = .get_eval_cotram_lower(order, tsupp),
    y_basis_fun_prime = .empty_fun(.get_eval_cotram(order, tsupp))
  )

  ret <- deeptrafo(formula = formula, data = data,
                   response_type = response_type, order = order,
                   addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                   monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                   ... = ...)

  class(ret) <- c("cotramNN", class(ret))

  ret

}

# Helpers -----------------------------------------------------------------

.empty_fun <- function(FUN) {
  function(...) {
    ret <- FUN(...)
    ret[] <- 1e-8
    ret
  }
}

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
