
context("Test deep ensembles")

# FUNs --------------------------------------------------------------------

check_ensemble <- function(formula = y ~ 1,
                           type = c("continuous", "ordinal",
                                    "count", "survival"),
                           which_ensemble = c("dr", "trf"),
                           ...) {

  type <- match.arg(type)
  which_ensemble <- match.arg(which_ensemble)

  dgp <- switch(
    type,
    "continuous" = dgp_cont,
    "ordinal" = dgp_ordinal,
    "count" = dgp_count,
    "survival" = dgp_surv
  )

  df <- dgp()

  if (which_ensemble == "dr") {
    m <- deeptrafo(formula, df, ... = ...)
    ens <- ensemble(m, n_ensemble = 2L, verbose = FALSE)
  } else if (which_ensemble == "trf") {
    ens <- trafoensemble(formula, df, ...)
  }

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

  ens <- check_ensemble(y ~ 1, type = "continuous", which_ensemble = "trf")
  ens <- check_ensemble(y ~ dnn(x), type = "continuous",
                        list_of_deep_models = list(dnn = dnn),
                        which_ensemble = "trf")

})

test_that("deep ensemble with ordinal outcome", {

  dnn <- function(x) x %>%
    layer_dense(units = 2L, activation = "relu") %>%
    layer_dense(units = 1L, use_bias = FALSE)

  ens <- check_ensemble(y ~ 1, type = "ordinal")
  ens <- check_ensemble(y ~ dnn(x), type = "ordinal",
                        list_of_deep_models = list(dnn = dnn))
  ens <- check_ensemble(y ~ 1, type = "ordinal", which_ensemble = "trf")
  ens <- check_ensemble(y ~ dnn(x), type = "ordinal",
                        list_of_deep_models = list(dnn = dnn),
                        which_ensemble = "trf")

})

test_that("deep ensemble with survival outcome", {

  dnn <- function(x) x %>%
    layer_dense(units = 2L, activation = "relu") %>%
    layer_dense(units = 1L, use_bias = FALSE)

  ens <- check_ensemble(y ~ 1, type = "survival")
  ens <- check_ensemble(y ~ dnn(x), type = "survival",
                        list_of_deep_models = list(dnn = dnn))
  ens <- check_ensemble(y ~ 1, type = "survival", which_ensemble = "trf")
  ens <- check_ensemble(y ~ dnn(x), type = "survival",
                        list_of_deep_models = list(dnn = dnn),
                        which_ensemble = "trf")

})

test_that("deep ensemble with count outcome", {

  dnn <- function(x) x %>%
    layer_dense(units = 2L, activation = "relu") %>%
    layer_dense(units = 1L, use_bias = FALSE)

  ens <- check_ensemble(y ~ 1, type = "count")
  ens <- check_ensemble(y ~ dnn(x), type = "count",
                        list_of_deep_models = list(dnn = dnn))
  ens <- check_ensemble(y ~ 1, type = "count", which_ensemble = "trf")
  ens <- check_ensemble(y ~ dnn(x), type = "count",
                        list_of_deep_models = list(dnn = dnn),
                        which_ensemble = "trf")

})

test_that("ensembles with callbacks, custom optimizers work", {
  df <- dgp_cont()
  ens <- trafoensemble(y ~ 1, data = df, optimizer = optimizer_adam(1),
                       n_ensemble = 3, seed = rep(1, 3), tf_seed = 1)
  cfb <- unname(c(unlist(coef(ens))))
  expect_equal(cfb, rep(cfb[1], 3))
  expect_equal(ll <- unname(unlist(logLik(ens))), unname(rep(ll[1], 5)))
  ens <- trafoensemble(y ~ 1, data = df, optimizer = optimizer_adam(1),
                       n_ensemble = 3, seed = rep(1, 3), tf_seed = 1,
                       callbacks = list(callback_early_stopping()),
                       validation_split = 0.05)
  expect_length(ens$ensemble_results[[2]]$metrics$loss, 2)
})
