
context("Test deep ensembles")

# FUNs --------------------------------------------------------------------

check_ensemble <- function(formula = y ~ 1,
                           type = c("continuous", "ordinal",
                                    "count", "survival"),
                           ...) {

  type <- match.arg(type)

  dgp <- switch(
    type,
    "continuous" = dgp_cont,
    "ordinal" = dgp_ordinal,
    "count" = dgp_count,
    "survival" = dgp_surv
  )

  df <- dgp()

  m <- deeptrafo(formula, df, ... = ...)

  ens <- ensemble(m, n_ensemble = 2L)

  check_ensemble_methods(ens)

}

check_ensemble_methods <- function(object) {

  expect_is(object$ensemble_results[[1]], "keras_training_history")
  expect_is(coef(object, which = "interacting")[[1]], "matrix")
  expect_is(coef(object, which = "shifting")[[1]], "matrix")
  expect_is(fitted(object), "list")

  invisible(object)

}

dgp_cont <- function(n = 100) {
  data.frame(y = rnorm(n),
             x = abs(rnorm(n)),
             z = rnorm(n))
}

dgp_ordinal <- function(ncl = 6L, n = 100) {
  data.frame(y = ordered(sample.int(ncl, n, replace = TRUE)),
             x = abs(rnorm(n)), z = rnorm(n))
}

dgp_count <- function(n = 100) {
  data.frame(
    y = sample.int(50, size = n, replace = TRUE),
    x = abs(rnorm(n)),
    z = rnorm(n),
    f = factor(sample.int(2, size = n, replace = TRUE))
  )
}

dgp_surv <- function(n = 100) {
  data.frame(
    y = survival::Surv(abs(rnorm(n, sd = 10)), sample(0:1, n, TRUE)),
    x = abs(rnorm(n)),
    z = rnorm(n),
    f = factor(sample.int(2, size = n, replace = TRUE))
  )
}

# Tests -------------------------------------------------------------------

test_that("deep ensemble with continuous outcome", {

  dnn <- function(x) x %>%
    layer_dense(units = 2L, activation = "relu") %>%
    layer_dense(units = 1L, use_bias = FALSE)

  ens <- check_ensemble(y ~ 1, type = "continuous")
  ens <- check_ensemble(y ~ dnn(x), type = "continuous",
                        list_of_deep_models = list(dnn = dnn))

})

test_that("deep ensemble with ordinal outcome", {

  dnn <- function(x) x %>%
    layer_dense(units = 2L, activation = "relu") %>%
    layer_dense(units = 1L, use_bias = FALSE)

  ens <- check_ensemble(y ~ 1, type = "ordinal")
  ens <- check_ensemble(y ~ dnn(x), type = "ordinal",
                        list_of_deep_models = list(dnn = dnn))

})

test_that("deep ensemble with survival outcome", {

  dnn <- function(x) x %>%
    layer_dense(units = 2L, activation = "relu") %>%
    layer_dense(units = 1L, use_bias = FALSE)

  ens <- check_ensemble(y ~ 1, type = "survival")
  ens <- check_ensemble(y ~ dnn(x), type = "survival",
                        list_of_deep_models = list(dnn = dnn))

})

test_that("deep ensemble with count outcome", {

  dnn <- function(x) x %>%
    layer_dense(units = 2L, activation = "relu") %>%
    layer_dense(units = 1L, use_bias = FALSE)

  ens <- check_ensemble(y ~ 1, type = "count")
  ens <- check_ensemble(y ~ dnn(x), type = "count",
                        list_of_deep_models = list(dnn = dnn))

})
