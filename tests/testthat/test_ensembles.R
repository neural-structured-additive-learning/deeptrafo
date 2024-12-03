context("Test deep ensembles")

if (.Platform$OS.type != "windows" &&
  reticulate::py_available() &&
  reticulate::py_module_available("tensorflow") &&
  reticulate::py_module_available("keras") &&
  reticulate::py_module_available("tensorflow_probability")) {
  # FUNs --------------------------------------------------------------------

  check_ensemble <- function(formula = y ~ 1,
                             type = c(
                               "continuous", "ordinal",
                               "count", "survival"
                             ),
                             which_ensemble = c("dr", "trf"),
                             ...) {
    type <- match.arg(type)
    which_ensemble <- match.arg(which_ensemble)

    dgp <- switch(type,
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
    data.frame(
      y = rnorm(n),
      x = abs(rnorm(n)),
      z = rnorm(n)
    )
  }

  dgp_ordinal <- function(ncl = 6L, n = 100) {
    data.frame(
      y = ordered(sample.int(ncl, n, replace = TRUE)),
      x = abs(rnorm(n)), z = rnorm(n)
    )
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
    dnn <- function(x) {
      x %>%
        layer_dense(units = 2L, activation = "relu") %>%
        layer_dense(units = 1L, use_bias = FALSE)
    }

    ens <- check_ensemble(y ~ 1, type = "continuous")
    ens <- check_ensemble(y ~ dnn(x),
      type = "continuous",
      list_of_deep_models = list(dnn = dnn)
    )

    ens <- check_ensemble(y ~ 1, type = "continuous", which_ensemble = "trf")
    ens <- check_ensemble(y ~ dnn(x),
      type = "continuous",
      list_of_deep_models = list(dnn = dnn),
      which_ensemble = "trf"
    )
  })

  test_that("deep ensemble with ordinal outcome", {
    dnn <- function(x) {
      x %>%
        layer_dense(units = 2L, activation = "relu") %>%
        layer_dense(units = 1L, use_bias = FALSE)
    }

    ens <- check_ensemble(y ~ 1, type = "ordinal")
    ens <- check_ensemble(y ~ dnn(x),
      type = "ordinal",
      list_of_deep_models = list(dnn = dnn)
    )
    ens <- check_ensemble(y ~ 1, type = "ordinal", which_ensemble = "trf")
    ens <- check_ensemble(y ~ dnn(x),
      type = "ordinal",
      list_of_deep_models = list(dnn = dnn),
      which_ensemble = "trf"
    )
  })

  test_that("deep ensemble with survival outcome", {
    dnn <- function(x) {
      x %>%
        layer_dense(units = 2L, activation = "relu") %>%
        layer_dense(units = 1L, use_bias = FALSE)
    }

    ens <- check_ensemble(y ~ 1, type = "survival")
    ens <- check_ensemble(y ~ dnn(x),
      type = "survival",
      list_of_deep_models = list(dnn = dnn)
    )
    ens <- check_ensemble(y ~ 1, type = "survival", which_ensemble = "trf")
    ens <- check_ensemble(y ~ dnn(x),
      type = "survival",
      list_of_deep_models = list(dnn = dnn),
      which_ensemble = "trf"
    )
  })

  test_that("deep ensemble with count outcome", {
    dnn <- function(x) {
      x %>%
        layer_dense(units = 2L, activation = "relu") %>%
        layer_dense(units = 1L, use_bias = FALSE)
    }

    ens <- check_ensemble(y ~ 1, type = "count")
    ens <- check_ensemble(y ~ dnn(x),
      type = "count",
      list_of_deep_models = list(dnn = dnn)
    )
    ens <- check_ensemble(y ~ 1, type = "count", which_ensemble = "trf")
    ens <- check_ensemble(y ~ dnn(x),
      type = "count",
      list_of_deep_models = list(dnn = dnn),
      which_ensemble = "trf"
    )
  })

  test_that("ensembles with callbacks, custom optimizers work", {
    df <- dgp_cont()
    ens <- trafoensemble(y ~ 1,
      data = df, optimizer = optimizer_adam(1),
      n_ensemble = 3, seed = rep(1, 3), tf_seeds = rep(1, 3)
    )
    cfb <- unname(c(unlist(coef(ens))))
    expect_equal(cfb, rep(cfb[1], 3))
    expect_equal(ll <- unname(unlist(logLik(ens))), unname(rep(ll[1], 5)))
    ens <- trafoensemble(y ~ 1,
      data = df, optimizer = optimizer_adam(1),
      n_ensemble = 3, seed = rep(1, 3), tf_seeds = rep(1, 3),
      callbacks = list(callback_early_stopping()),
      validation_split = 0.05
    )
    expect_length(ens$ensemble_results[[2]]$metrics$loss, 2)
  })

  test_that("reinit works as expected", {
    df <- dgp_cont()

    ### No special args
    ens <- trafoensemble(y ~ 1,
      data = df, n_ensemble = 3, seed = c(2, 10, 20),
      tf_seeds = c(2, 10, 20),
      epochs = 3, validation_split = 0.1
    )
    expect_false(any(ens$ensemble_results[[1]]$metrics$loss ==
      ens$ensemble_results[[2]]$metrics$loss))

    ### Custom optimizer
    ens <- trafoensemble(y ~ 1,
      data = df, optimizer = optimizer_adam(1e-4),
      n_ensemble = 3, seed = c(2, 10, 20), tf_seeds = c(2, 10, 20),
      epochs = 3,
      validation_split = 0.1
    )
    expect_false(any(ens$ensemble_results[[1]]$metrics$loss ==
      ens$ensemble_results[[2]]$metrics$loss))

    ### Custom callback
    ens <- trafoensemble(y ~ 1,
      data = df, optimizer = optimizer_adam(1e-4),
      n_ensemble = 3, seed = c(2, 10, 20), tf_seeds = c(2, 10, 20),
      epochs = 3, callbacks = callback_early_stopping(patience = 0),
      validation_split = 0.1
    )
    expect_false(any(ens$ensemble_results[[1]]$metrics$loss ==
      ens$ensemble_results[[2]]$metrics$loss))

    ### Custom neural network
    ens <- trafoensemble(y ~ nn(x),
      data = df, optimizer = optimizer_adam(1e-4),
      n_ensemble = 3, seed = c(2, 10, 20), tf_seeds = c(2, 10, 20),
      epochs = 3,
      callbacks = callback_early_stopping(patience = 0),
      validation_split = 0.1, list_of_deep_models = list(
        nn = \(x) x %>% layer_dense(1, activation = "relu")
      )
    )
    expect_false(any(ens$ensemble_results[[1]]$metrics$loss ==
      ens$ensemble_results[[2]]$metrics$loss))

    ### No special args + seed
    ens <- trafoensemble(y ~ 1,
      data = df, n_ensemble = 3, seed = rep(1, 3),
      tf_seeds = rep(1, 3), epochs = 3, validation_split = 0.1
    )
    expect_true(all(ens$ensemble_results[[1]]$metrics$loss ==
      ens$ensemble_results[[2]]$metrics$loss))

    ### Custom optimizer + seed
    ens <- trafoensemble(y ~ 1,
      data = df, optimizer = optimizer_adam(1e-4),
      n_ensemble = 3, seed = rep(1, 3), tf_seeds = rep(1, 3),
      epochs = 3, validation_split = 0.1
    )
    expect_true(all(ens$ensemble_results[[1]]$metrics$loss ==
      ens$ensemble_results[[2]]$metrics$loss))

    ### Custom callback + seed
    ens <- trafoensemble(y ~ 1,
      data = df, optimizer = optimizer_adam(1e-4),
      n_ensemble = 3, seed = rep(1, 3), tf_seeds = rep(1, 3),
      epochs = 3, callbacks = callback_early_stopping(patience = 0),
      validation_split = 0.1
    )
    expect_true(all(ens$ensemble_results[[1]]$metrics$loss ==
      ens$ensemble_results[[2]]$metrics$loss))

    ### Custom neural network + seed
    ens <- trafoensemble(y ~ nn(x),
      data = df, optimizer = optimizer_adam(1e-4),
      n_ensemble = 3, seed = rep(1, 3), tf_seeds = rep(1, 3),
      epochs = 3, callbacks = callback_early_stopping(patience = 0),
      validation_split = 0.1, list_of_deep_models = list(
        nn = \(x) x %>% layer_dense(1, activation = "relu")
      )
    )
    expect_true(all(ens$ensemble_results[[1]]$metrics$loss ==
      ens$ensemble_results[[2]]$metrics$loss))
  })

  test_that("weighted transformation ensemble", {
    set.seed(123)

    dgp <- function(n = 1e2) {
      x <- runif(n)
      y <- sin(x * pi * 3) + rnorm(n, sd = 0.1)
      data.frame(y = y, x = x)
    }

    ### Generate train, validation, test data
    train_data <- dgp()
    validation_data <- dgp()
    test_data <- dgp()

    ### Train ensemble with early stopping (separate validation split)
    ens <- trafoensemble(
      y ~ x,
      data = train_data, epochs = 10, n_ensemble = 3, verbose = FALSE,
      validation_split = 0.1, optimizer = optimizer_adam(learning_rate = 0.1),
      callbacks = list(callback_early_stopping(patience = 50))
    )

    ### Compute the optimal weights on the validation data (weights = NULL; default)
    tuned <- weighted_logLik(ens, newdata = validation_data)

    expect_true(all(tuned$weights <= 1))
    expect_true(all(tuned$weights >= 0))
    expect_length(tuned$weights, 3L)
    expect_lt(tuned$ensemble, tuned$mean)

    ### Use optimal weights and make test predictions (weights = tuned$weights)
    test_nll <- weighted_logLik(ens, weights = tuned$weights, newdata = test_data)

    expect_equal(tuned$weights, test_nll$weights)
    expect_lt(test_nll$ensemble, test_nll$mean)
  })
}
