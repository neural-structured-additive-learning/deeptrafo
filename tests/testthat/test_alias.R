
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
  dat$ylag2 <- lag(dat$y, n = 2L)
  dat <- na.omit(dat)
  m <- dctm(~ y, ~ s(x), ~ z + s(z) + atplag(ylag) + atplag(ylag2), data = dat)
  hist <- fit(m, epochs = 2L)

  expect_false(any(is.nan(hist$metrics$loss)))

})

# Gompertz ----------------------------------------------------------------

test_that("gompertz base distribution works for ordinal case", {
  library(tram)

  data("wine", package = "ordinal")
  wine$noise <- rnorm(nrow(wine))
  fml <- rating ~ 0 + temp
  optimizer <- optimizer_adam(learning_rate = 0.1, decay = 5e-4)
  m <- deeptrafo(fml, wine, latent_distr = "gompertz", monitor_metric = NULL,
                 optimizer = optimizer)
  m %>% fit(epochs = 300, batch_size = nrow(wine), validation_split = 0)
  llm <- logLik(m)

  lltm <- c(logLik(tm <- Polr(rating ~ temp, data = wine, method = "cloglog")))

  expect_equal(llm, lltm, tolerance = 1e-3)
})

# Coxph -------------------------------------------------------------------

test_that("CoxphNN gives same results as tram::Coxph", {
  set.seed(1)

  library(tram)
  tord <- 3

  d <- dgp_surv()
  tm <- Coxph(y ~ x, data = d, order = tord, support = range(d$y[, 1]))
  m <- CoxphNN(y ~ 0 + x, data = d, order = tord)

  tmp <- get_weights(m$model)
  .to_gamma <- function(thetas) {
    gammas <- c(thetas[1L], log(exp(diff(thetas)) - 1))
    if (any(!is.finite(gammas))) {
      gammas[!is.finite(gammas)] <- 1e-20
    }
    return(gammas)
  }

  tmp[[1]][] <- .to_gamma(coef(tm, with_baseline = TRUE)[1:(tord + 1)])
  tmp[[2]][] <- coef(tm)
  set_weights(m$model, tmp)

  llm <- - logLik(m)
  lltm <- - c(logLik(tm))

  expect_equal(llm, lltm, tol = 1e-3)
})

# Lm ----------------------------------------------------------------------

test_that("LmNN gives same results as tram::Lm", {

  set.seed(1)
  library(tram)
  df <- data.frame(y = 10 + rnorm(50), x = rnorm(50))
  optimizer <- optimizer_adam(learning_rate = 0.1, decay = 4e-4)
  m <- LmNN(y ~ 0 + x, data = df, optimizer = optimizer)
  # fit(m, epochs = 1800L, validation_split = 0)
  mm <- Lm(y ~ x, data = df, support = range(df$y))
  cfb <- coef(mm, with_baseline = TRUE)

  tmp <- get_weights(m$model)
  tmp[[1]][] <- c(cfb[1], log(exp(cfb[2]) - 1))
  tmp[[2]][] <- - cfb[3]
  set_weights(m$model, tmp)

  expect_equal(c(logLik(mm)), logLik(m), tol = 1e-3)

})

# Survreg -----------------------------------------------------------------

test_that("SurvregNN gives same results as tram::Survreg", {

  set.seed(1)
  library(tram)

  d <- dgp_surv()
  tm <- Survreg(y ~ x, data = d, support = range(d$y[, 1]))
  m <- SurvregNN(y ~ 0 + x, data = d)

  cfb <- coef(tm, with_baseline = TRUE)

  tmp <- get_weights(m$model)
  tmp[[1]][] <- c(cfb[1], log(exp(cfb[2]) - 1))
  tmp[[2]][] <- - cfb[3]
  set_weights(m$model, tmp)

  llm <- - logLik(m)
  lltm <- - c(logLik(tm))

  expect_equal(llm, lltm, tol = 1e-3)

})
