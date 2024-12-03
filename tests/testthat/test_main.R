context("Test deeptrafo")

if (.Platform$OS.type != "windows" &&
  reticulate::py_available() &&
  reticulate::py_module_available("tensorflow") &&
  reticulate::py_module_available("keras") &&
  reticulate::py_module_available("tensorflow_probability")) {
  # source("tests/testthat/test-funs.R")
  source("test-funs.R")

  # Additive models ---------------------------------------------------------

  test_that("simple additive model", {
    dat <- data.frame(
      y = rnorm(100), x = rnorm(100), z = rnorm(100),
      f = factor(sample(0:1, 100, TRUE))
    )
    fml <- y | f ~ z + s(z)
    m <- deeptrafo(fml, dat)

    check_methods(m, newdata = dat, test_plots = FALSE)
  })

  test_that("unconditional additive model", {
    dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100))
    valdat <- data.frame(y = rcauchy(100), x = rcauchy(100), z = rcauchy(100))
    fml <- y ~ 1
    m <- deeptrafo(fml, dat)
    hist <- fit(m,
      epochs = 10, validation_data = list(x = valdat, y = valdat$y),
      verbose = FALSE
    )
    expect_false(any(is.nan(hist$metrics$loss)))

    check_methods(m, newdata = dat, test_plots = FALSE)
  })

  # Ordinal -----------------------------------------------------------------

  test_that("unconditional ordinal model", {
    test_models(y ~ 1)
  })

  test_that("ordinal model", {
    test_models(y ~ x)
  })

  test_that("ordinal model with smooth effects", {
    test_models(y ~ s(z))
  })

  test_that("ordinal model with response-varying effects", {
    test_models(y | x ~ s(z))
  })

  test_that("monotonicity problem (ordinal case)", {
    test_models(y | s(x) ~ z)
  })

  test_that("ordinal model with NN component", {
    nn <- keras_model_sequential() %>%
      layer_dense(input_shape = 1L, units = 6L, activation = "relu") %>%
      layer_dense(units = 1L)

    test_models(y ~ nn(x), list_of_deep_models = list(nn = nn))
  })

  test_that("ordinal NLL works", {
    df <- data.frame(y = ordered(rep(1:5, each = 5)))
    m <- deeptrafo(y ~ 1, data = df)
    fit(m,
      validation_split = NULL, epochs = 10, batch_size = nrow(df),
      verbose = FALSE
    )
    # coef(m); coef(m, "interacting")

    cf0 <- qlogis((1:4) / 5)
    ll0 <- -nrow(df) * log(1 / 5)

    sp_inv <- function(x) c(x[1], log(exp(diff(x)) - 1), -Inf)

    tmp <- get_weights(m$model)
    tmp[[2]][] <- 0.0
    tmp[[1]][] <- sp_inv(cf0)
    set_weights(m$model, tmp)

    cf <- coef(m, which = "interacting")

    tloss <- nll("logistic")
    ll <- tloss(m$init_params$y, fitted(m))$numpy()

    expect_equal(ll0, sum(ll), tolerance = 1e-5)
    expect_equal(cf0, unname(unlist(cf))[1:4], tol = 1e-4)
  })

  # Count models ------------------------------------------------------------

  test_that("unconditional count model", {
    test_models(y ~ 1, which = "count")
  })

  test_that("count model", {
    test_models(y ~ x, which = "count")
  })

  test_that("count model with smooth effects", {
    test_models(y ~ s(z), which = "count")
  })

  test_that("count model with response-varying effects", {
    test_models(y | f ~ s(z), which = "count")
  })

  test_that("monotonicity problem (count case)", {
    test_models(y | s(x) ~ z, which = "count")
  })

  test_that("count model with NN component", {
    nn <- keras_model_sequential() %>%
      layer_dense(input_shape = 1L, units = 6L, activation = "relu") %>%
      layer_dense(units = 1L)

    test_models(y ~ nn(x), list_of_deep_models = list(nn = nn), which = "count")
  })

  # Survival models --------------------------------------------------------

  test_that("unconditional survival model", {
    test_models(y ~ 1, which = "survival")
  })

  test_that("survival model", {
    test_models(y ~ x, which = "survival")
  })

  test_that("survival model with smooth effects", {
    test_models(y ~ s(z), which = "survival")
  })

  test_that("survival model with response-varying effects", {
    test_models(y | f ~ s(z), which = "survival")
  })

  test_that("monotonicity problem (survival case)", {
    test_models(y | s(x) ~ z, which = "survival")
  })

  test_that("survival model with NN component", {
    nn <- keras_model_sequential() %>%
      layer_dense(input_shape = 1L, units = 6L, activation = "relu") %>%
      layer_dense(units = 1L)

    test_models(y ~ nn(x), list_of_deep_models = list(nn = nn), which = "survival")
    test_models(y | nn(x) ~ 1, list_of_deep_models = list(nn = nn), which = "survival")
  })

  # Autoregressive models ---------------------------------------------------

  test_that("autoregressive transformation model", {
    dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100))
    fml <- y | s(x) ~ 0 + s(z) + atplag(1, 2) # max lag (i.e. 2) reduces data set (also at predict)
    m <- deeptrafo(fml, dat)

    expect_is(predict(m, newdata = dat[1:10, -1], K = 2, type = "pdf"), "list")
    expect_is(predict(m, newdata = dat[1:10, -1], q = c(-1, 1), type = "pdf"), "list")

    check_methods(m, newdata = dat)

    cf <- coef(m, which_param = "autoregressive")
    expect_equal(length(cf), 2)
  })

  test_that("autoregressive count transformation model", {
    dat <- data.frame(y = round(rnorm(100, mean = 1e3, sd = 20)), x = rnorm(100), z = rnorm(100))
    fml <- y | s(x) ~ 0 + s(z) + atplag(1, 2) # max lag (i.e. 2) reduces data set (also at predict)
    m <- deeptrafo(fml, dat)

    expect_is(predict(m, newdata = dat[1:10, -1], K = 2, type = "pdf"), "list")
    expect_is(predict(m, newdata = dat[1:10, -1], q = range(dat$y), type = "pdf"), "list")

    check_methods(m, newdata = dat)

    cf <- coef(m, which_param = "autoregressive")
    expect_equal(length(cf), 2)
  })

  test_that("autoregressive transformation model specification", {
    dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100))
    dat <- na.omit(dat)
    expect_length(coef(deeptrafo(y ~ atplag(1), data = dat), which = "auto"), 1)
    expect_length(coef(deeptrafo(y ~ atplag(1:2), data = dat), which = "auto"), 2)
  })

  # Misc --------------------------------------------------------------------

  test_that("model with fixed weight", {
    data("wine", package = "ordinal")
    m <- deeptrafo(response ~ temp,
      data = wine,
      weight_options = weight_control(
        warmstart_weights = list(list(), list(), list("temp" = 0))
      )
    )
    expect_equal(unname(coef(m, which_param = "shifting")$temp[1, 1]), 0)
  })

  # Deep --------------------------------------------------------------------

  test_that("deep conditional model", {
    dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100))

    deep_model <- function(x) {
      x %>%
        layer_dense(units = 32, activation = "relu", use_bias = FALSE) %>%
        layer_dropout(rate = 0.2) %>%
        layer_dense(units = 8, activation = "relu")
    }

    fml <- y | d(x) ~ z + s(z)
    m <- deeptrafo(fml, dat, list_of_deep_models = list(d = deep_model))

    check_methods(m, dat[1:10, ], FALSE, FALSE)
  })

  # Shared -----------------------------------------------------------------

  # test_that("shared model", {
  #
  #   dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100))
  #
  #   deep_model <- function(x) x %>%
  #     layer_dense(units = 32, activation = "relu", use_bias = FALSE) %>%
  #     layer_dropout(rate = 0.2) %>%
  #     layer_dense(units = 8, activation = "relu")
  #
  #   fml <- y | x ~ z + s(z) | d(x)
  #   m <- deeptrafo(fml, dat, list_of_deep_models = list(d = deep_model),
  #                  shared_partition = 7)
  #
  #   check_methods(m, dat[1:10, ], FALSE, FALSE)
  #
  # })
}
