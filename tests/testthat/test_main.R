context("Test deeptrafo")

devtools::load_all("../../../deepregression/")

check_methods <- function(m, newdata)
{

  # fit
  hist <- m %>% fit(epochs = 2)
  expect_is(hist, "keras_training_history")

  # plot
  pret1 <- plot(m, which_param = "h1")
  expect_is(pret1, "list")
  pret2 <- plot(m, which_param = "h2")
  expect_is(pret2, "list")

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


}

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

test_that("unconditional ordinal model", {

  dat <- data.frame(y = ordered(sample(1:6, 100, replace = TRUE)),
                    x = rnorm(100), z = rnorm(100))
  fml <- y ~ 1
  m <- deeptrafo(fml, dat)
  hist <- fit(m, epochs = 2L)
  expect_equal(m$init_params$trafo_options$order_bsp, 5L)
  expect_false(any(is.nan(hist$metrics$loss)))

})

test_that("ordinal model", {

  ncl <- 6L
  dat <- data.frame(y = ordered(sample(1:ncl, 100, replace = TRUE)),
                    x = rnorm(100), z = rnorm(100))
  fml <- y ~ s(z)
  m <- deeptrafo(fml, dat)
  hist <- fit(m, epochs = 2L)
  expect_equal(m$init_params$trafo_options$order_bsp, ncl - 1L)
  expect_false(any(is.nan(hist$metrics$loss)))

})

test_that("autoregressive transformation model", {

  dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100))
  dat$ylag <- lag(dat$y)
  dat$ylag2 <- lag(dat$y, n=2L)
  dat <- na.omit(dat)
  fml <- y | s(x) ~ z + s(z)
  m <- deeptrafo(fml, dat, lag_formula = ~ ylag + ylag2)

  check_methods(m, newdata = dat)

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

  expect_equal(ll0, ll, tolerance = 1e-5)
  expect_equal(cf0, unname(unlist(cf))[1:4], tol = 1e-4)

})
