
# devtools::load_all("../deepregression/")
# devtools::load_all(".")
devtools::load_all("../../../deepregression/")

context("Test aliases")

# source("tests/testthat/test-funs.R")
source("test-funs.R")

test_alias <- function(rsp, int = NULL, shi = NULL, FUN = dctm,
                       which = c("ordinal", "count", "survival"), ...) {

  which <- match.arg(which)

  DGP <- switch(which,
    "ordinal" = dgp_ordinal,
    "count" = dgp_count,
    "survival" = dgp_surv
  )

  dat <- DGP()
  m <- FUN(response = rsp, intercept = int, shift = shi, data = dat, ...)

  if (which == "ordinal")
    expect_false(any(is.nan(m$model$loss(t(sapply(dat$y, eval_ord)),
                                         fitted(m))$numpy())))
  hist <- fit(m, epochs = 2L)

  if (which == "ordinal")
    expect_equal(m$init_params$trafo_options$order_bsp, 5L)

  expect_false(any(is.nan(hist$metrics$loss)))

}

# Alias -------------------------------------------------------------------

test_that("simple additive model", {

  dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100),
                    f = factor(sample(0:1, 100, TRUE)))

  # DCTM
  m <- dctm(response = ~ y, intercept = ~ f, shift = ~ 0 + z + s(z), data = dat)
  check_methods(m, newdata = dat, test_plots = FALSE)

  # Tram-like aliases
  m <- BoxCoxNN(y | f ~ z + s(z), data = dat)
  check_methods(m, newdata = dat, test_plots = FALSE)
  m <- LehmanNN(y | f ~ z + s(z), data = dat)
  check_methods(m, newdata = dat, test_plots = FALSE)
  m <- ColrNN(y | f ~ z + s(z), data = dat)
  check_methods(m, newdata = dat, test_plots = FALSE)
  expect_error(PolrNN(y | f ~ z + s(z), data = dat))

})

# Ordinal -----------------------------------------------------------------

test_that("unconditional ordinal model", {

  test_alias(~ y)

})

test_that("ordinal model", {

  test_alias(~ y, NULL, ~ x)
  test_alias(~ y, NULL, ~ x, FUN = ontram)

})

test_that("count model with NN component", {

  nn <- keras_model_sequential() %>%
    layer_dense(input_shape = 1L, units = 6L, activation = "relu") %>%
    layer_dense(units = 1L)

  test_alias(~ y, NULL, ~ nn(x), list_of_deep_models = list(nn = nn), which = "count")

})

test_that("survival model with response-varying effects", {

  test_alias(~ y, ~ f, ~ s(z), which = "survival")

})

test_that("autoregressive transformation model", {

  dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100))
  dat$ylag <- lag(dat$y)
  dat$ylag2 <- lag(dat$y, n=2L)
  dat <- na.omit(dat)
  m <- dctm(~ y, ~ s(x), ~ z + s(z), data = dat, lag_formula = ~ ylag + ylag2)
  hist <- fit(m, epochs = 2L)

  expect_false(any(is.nan(hist$metrics$loss)))

})
