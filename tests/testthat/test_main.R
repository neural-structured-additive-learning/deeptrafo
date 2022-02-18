
# devtools::load_all("../deepregression/")
# devtools::load_all(".")
# devtools::load_all("../../../deepregression/")

context("Test deeptrafo")

# FUNs --------------------------------------------------------------------

check_methods <- function(m, newdata, test_plots = TRUE)
{

  # fit
  hist <- m %>% fit(epochs = 2)
  expect_is(hist, "keras_training_history")

  # plot
  if (test_plots) {
    pret1 <- plot(m, which_param = "interacting")
    expect_is(pret1, "list")
    pret2 <- plot(m, which_param = "shifting")
    expect_is(pret2, "list")
  }

  # coef
  ch1 <- coef(m, which = "interacting")
  expect_is(ch1, "list")
  ch2 <- coef(m, which = "shifting")
  expect_is(ch2, "list")

  # fitted
  fitt <- m %>% fitted()
  expect_is(fitt, "matrix")

  # predict
  a <- predict(m, newdata = newdata, type = "trafo")
  this_n <- nrow(newdata)
  expect_equal(dim(a), c(this_n, 1))
  b <- predict(m, newdata = newdata, type = "pdf")
  expect_equal(dim(b), c(this_n, 1))
  expect_true(all(b >= 0))
  c <- predict(m, newdata = newdata, type = "cdf")
  expect_equal(dim(c), c(this_n, 1))
  expect_true(all(c >= 0) & all(c <= 1))
  d <- predict(m, newdata = newdata, type = "interaction")
  expect_equal(dim(d), c(this_n, 1 + ncol(newdata)))
  e <- predict(m, newdata = newdata, type = "shift")
  expect_equal(dim(e), c(this_n, 1 + ncol(newdata)))
  f <- predict(m, newdata = newdata, type = "output")
  expect_equal(nrow(f), this_n, 1)
  expect_gt(ncol(f), 2)
  # g <- predict(m, newdata = newdata[, colnames(newdata) != "y"], type = "trafo")
  # expect_equal(dim(g), c(this_n, this_n))
  # h <- predict(m, newdata = newdata[, colnames(newdata) != "y"], type = "pdf")
  # expect_equal(dim(h), c(this_n, this_n))
  # expect_true(all(h >= 0))
  # i <- predict(m, newdata = newdata[, colnames(newdata) != "y"], type = "cdf")
  # expect_true(all(i >= 0) & all(i <= 1))
  # expect_equal(dim(i), c(this_n, this_n))

  # logLik
  expect_is(logLik(m), "numeric")

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

test_models <- function(fml, which = c("ordinal", "count", "survival"), ...) {

  which <- match.arg(which)

  DGP <- switch(which,
    "ordinal" = dgp_ordinal,
    "count" = dgp_count,
    "survival" = dgp_surv
  )

  dat <- DGP()
  m <- deeptrafo(fml, dat, ...)

  if (which == "ordinal")
    expect_false(any(is.nan(m$model$loss(t(sapply(dat$y, eval_ord)),
                                         fitted(m))$numpy())))
  hist <- fit(m, epochs = 2L)

  if (which == "ordinal")
    expect_equal(m$init_params$trafo_options$order_bsp, 5L)

  expect_false(any(is.nan(hist$metrics$loss)))

  check_methods(m, dat, test_plots = FALSE)

}

# Additive models ---------------------------------------------------------

test_that("simple additive model", {

  dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100),
                    f = factor(sample(0:1, 100, TRUE)))
  fml <- y | f ~ z + s(z)
  m <- deeptrafo(fml, dat)

  check_methods(m, newdata = dat, test_plots = FALSE)

})

test_that("unconditional additive model", {

  dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100))
  fml <- y ~ 1
  m <- deeptrafo(fml, dat)
  hist <- fit(m, epochs = 2L)
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

  # test_models(y | s(x) ~ z)

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
  fit(m, validation_split = NULL, epochs = 10, batch_size = nrow(df))
  coef(m); coef(m, "shifting")

  cf0 <- qlogis((1:4)/5)
  ll0 <- - nrow(df) * log(1/5)

  sp_inv <- function(x) c(x[1], log(exp(diff(x)) - 1), -Inf)

  tmp <- get_weights(m$model)
  tmp[[2]][] <- 0.0
  tmp[[1]][] <- sp_inv(cf0)
  set_weights(m$model, tmp)

  cf <- coef(m, which = "interacting")

  tloss <- nll_ordinal()
  ll <- tloss(response(df$y), fitted(m))$numpy()

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

  # test_models(y | s(x) ~ z, which = "count")

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

  # test_models(y | s(x) ~ z, which = "survival")

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
  dat$ylag <- lag(dat$y)
  dat$ylag2 <- lag(dat$y, n=2L)
  dat <- na.omit(dat)
  fml <- y | s(x) ~ z + s(z)
  m <- deeptrafo(fml, dat, lag_formula = ~ ylag + ylag2)

  check_methods(m, newdata = dat)

})

# Misc --------------------------------------------------------------------

test_that("model with fixed weight", {

  data("wine", package = "ordinal")
  m <- deeptrafo(response ~ temp, data = wine,
                 weight_options = weight_control(
                   warmstart_weights = list(list(), list(), list("temp" = 0))
                   )
                 )
  expect_equal(coef(m, which_param = "shifting")$temp, matrix(0))

})
