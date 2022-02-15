context("Test deeptrafo")

# devtools::load_all("../../../deepregression/")

check_methods <- function(m, newdata, test_plots = TRUE)
{

  # fit
  hist <- m %>% fit(epochs = 2)
  expect_is(hist, "keras_training_history")

  # plot
  if (test_plots) {
    pret1 <- plot(m, which_param = "h1")
    expect_is(pret1, "list")
    pret2 <- plot(m, which_param = "h2")
    expect_is(pret2, "list")
  }

  # coef
  ch1 <- coef(m, which = "h1")
  expect_is(ch1, "list")
  ch2 <- coef(m, which = "h2")
  expect_is(ch2, "list")

  # fitted
  fitt <- m %>% fitted()
  expect_is(fitt, "matrix")

  # predict
  trf_fun <- m %>% predict.deeptrafo(newdata)
  expect_is(trf_fun, "function")
  a <- trf_fun(newdata$y)
  this_n <- nrow(newdata)
  expect_equal(dim(a), c(this_n,1))
  b <- trf_fun(newdata$y, type = "pdf")
  expect_equal(dim(b), c(this_n,1))
  c <- trf_fun(newdata$y, type = "cdf")
  expect_equal(dim(c), c(this_n,1))
  d <- trf_fun(newdata$y, type = "interaction")
  expect_equal(dim(d), c(this_n,1+ncol(newdata)))
  e <- trf_fun(newdata$y, type = "shift")
  expect_equal(dim(e), c(this_n,1+ncol(newdata)))
  f <- trf_fun(newdata$y, type = "output")
  expect_equal(nrow(f), this_n,1)
  expect_gt(ncol(f), 2)
  g <- trf_fun(newdata$y, type = "trafo", grid=TRUE)
  expect_equal(dim(g), c(this_n,this_n))
  h <- trf_fun(newdata$y, type = "pdf", grid=TRUE)
  expect_equal(dim(h), c(this_n,this_n))
  i <- trf_fun(newdata$y, type = "cdf", grid=TRUE)
  expect_equal(dim(i), c(this_n,this_n))

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

test_models <- function(fml, which = c("ordinal", "count"), ...) {

  which <- match.arg(which)

  DGP <- switch(which,
    "ordinal" = dgp_ordinal,
    "count" = dgp_count
  )

  dat <- DGP()
  m <- deeptrafo(fml, dat, ...)

  if (which == "ordered")
    expect_false(any(is.nan(m$model$loss(t(sapply(dat$y, eval_ord)),
                                         fitted(m))$numpy())))
  hist <- fit(m, epochs = 2L)

  if (which == "ordered")
    expect_equal(m$init_params$trafo_options$order_bsp, ncl - 1L)

  expect_false(any(is.nan(hist$metrics$loss)))

  check_methods(m, dat, test_plots = FALSE)

}

# Additive models ---------------------------------------------------------

test_that("simple additive model", {

  dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100))
  fml <- y | s(x) ~ z + s(z)
  m <- deeptrafo(fml, dat)

  check_methods(m, newdata = dat)

})

test_that("unconditional additive model", {

  dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100))
  fml <- y ~ 1
  m <- deeptrafo(fml, dat)
  hist <- fit(m, epochs = 2L)
  expect_false(any(is.nan(hist$metrics$loss)))

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
  fit(m, validation_split = NULL, epochs = 10, batch_size = nrow(df))
  coef(m); coef(m, "h2")

  cf0 <- qlogis((1:4)/5)
  ll0 <- - nrow(df) * log(1/5)

  sp_inv <- function(x) c(x[1], log(exp(diff(x)) - 1), -Inf)

  tmp <- get_weights(m$model)
  tmp[[2]][] <- 0.0
  tmp[[1]][] <- sp_inv(cf0)
  set_weights(m$model, tmp)

  cf <- coef(m, which = "h1")

  tloss <- nll_ordinal()
  ll <- tloss(t(sapply(df$y, eval_ord)), fitted(m))$numpy()

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

# Autoregressive models ---------------------------------------------------

test_that("autoregressive transformation model", {

  dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100))
  dat$ylag <- lag(dat$y)
  dat$ylag2 <- lag(dat$y, n=2L)
  dat <- na.omit(dat)
  fml <- y ~ s(x) | z + s(z)
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
  expect_equal(coef(m, which_param = "h2")$temp, matrix(0))

})
