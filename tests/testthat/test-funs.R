# Test FUNs

# FUNs --------------------------------------------------------------------

check_methods <- function(m, newdata, test_plots = TRUE, grid = TRUE)
{

  # fit
  hist <- m %>% fit(epochs = 2, verbose = FALSE)
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

  numb_lags <- 0
  if (m$init_params$is_atm) {
    numb_lags <- max(deeptrafo:::fm_to_lag(m$init_params$lag_formula))
  }

  # predict
  a <- predict(m, newdata = newdata, type = "trafo")
  this_n <- nrow(newdata)
  expect_equal(dim(a), c(this_n - numb_lags, 1))
  b <- predict(m, newdata = newdata, type = "pdf")
  expect_equal(dim(b), c(this_n - numb_lags, 1))
  # expect_true(all(b >= 0))
  c <- predict(m, newdata = newdata, type = "cdf")
  expect_equal(dim(c), c(this_n - numb_lags, 1))
  expect_true(all(c >= 0) & all(c <= 1))
  d <- predict(m, newdata = newdata, type = "interaction")
  expect_equal(dim(d), c(this_n - numb_lags, 1))
  e <- predict(m, newdata = newdata, type = "shift")
  expect_equal(dim(e), c(this_n - numb_lags, 1))
  f <- predict(m, newdata = newdata, type = "terms")
  expect_equal(nrow(f), this_n - numb_lags, 1)
  expect_gt(ncol(f), 2)

  # if (m$init_params$response_type == "ordered") {
  # simulate(m)
  # simulate(m, newdata = newdata)
  # simulate(m, newdata = newdata[1:10, ])
  # }

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
    expect_false(any(is.nan(m$model$loss(m$init_params$y, fitted(m))$numpy())))
  hist <- fit(m, epochs = 2L, verbose = FALSE)

  if (which == "ordinal")
    expect_equal(m$init_params$trafo_options$order_bsp, 5L)

  expect_false(any(is.nan(hist$metrics$loss)))

  check_methods(m, dat, test_plots = FALSE)

}
