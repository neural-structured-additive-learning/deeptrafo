context("Test deeptrafo helperfunctions")

if (.Platform$OS.type != "windows" &&
  reticulate::py_available() &&
  reticulate::py_module_available("tensorflow") &&
  reticulate::py_module_available("keras") &&
  reticulate::py_module_available("tensorflow_probability")) {
  # response helpers --------------------------------------------------------

  test_that("response types", {
    expect_equal(get_response_type(ordered(1)), "ordered")
    expect_equal(get_response_type(1L), "count")
    expect_equal(get_response_type(1.1), "continuous")
  })

  test_that("eval response types", {
    expect_equal(eval_response(ordered(1:2), "ordered"), diag(2))
    expect_equal(eval_response(0:2, "count"), cbind(c(1, 0, 0), y = 0:2))
    expect_equal(eval_response(c(1.1, 3.2), "continuous"), c(1.1, 3.2))
  })

  # test_that("make grids of responses", {
  #   tn <- 50
  #   expect_length(make_grid(rnorm(10), tn)$y, tn)
  #   expect_length(make_grid(1:10)$y, 10L)
  #   expect_length(make_grid(ordered(1:5))$y, 5L)
  #   expect_length(make_grid(survival::Surv(rchisq(10, df = 1), rep(0, 10)), n = tn)$y, tn)
  # })

  # Formula helpers ---------------------------------------------------------

  test_that("formula parts work", {
    expect_equal(forms2form(~y, NULL, NULL, NULL), y ~ 1)
    expect_equal(forms2form(~y, ~x, NULL, NULL), y | x ~ 1)
    expect_equal(forms2form(~y, NULL, ~x, NULL), y ~ x)
    expect_equal(forms2form(~y, ~z, ~x, NULL), y | z ~ x)
    expect_equal(
      forms2form(~y, ~ z + dnn(abc), ~ x + f(g) + lasso(d), NULL),
      y | z + dnn(abc) ~ x + f(g) + lasso(d)
    )
    expect_equal(forms2form(~y, ~x, ~x, ~x), y | x ~ x | x)
    expect_error(forms2form(NULL, ~x)) # response can't be NULL
  })

  test_that("formula is printed correctly", {
    df <- data.frame(y = 1:10, x = 1:10, z = 1:10)
    expect_true(any(grepl("y | 1", capture.output(deeptrafo(y ~ x, data = df)))))
    expect_true(any(grepl("y | z", capture.output(deeptrafo(y | z ~ x, data = df)))))
    expect_true(any(grepl("y | s(z)", capture.output(deeptrafo(y | s(z) ~ 0 + x, data = df)))))
  })


  # ATM helpers -------------------------------------------------------------

  # test_that("atm_lag helpers", {
  #
  #   # form_lin <- ~ 1 + lag1 + lag2
  #   # form_nlin <- ~ 1 + s(lag1) + s(lag2)
  #   # form_mixed <- ~ 1 + s(lag1) + lag2
  #   # form_noint <- ~ 0 + lag1 + lag2
  #   # form_onevar <- ~ lag1
  #   # form_onesm <- ~ s(lag1)
  #   # form_te <- ~ te(lag1, lag2)
  #   # form_mixednl <- ~ te(lag1, lag2) + s(lag3)
  #   # form_lass <- ~ lasso(lag1, lag2, lag3)
  #   # form_lass2 <- ~ lasso(lag1 + lag2 + lag3)
  #   #
  #   # res_form_lin <- apply_atm_lags(form_lin)
  #   # res_form_nlin <- apply_atm_lags(form_nlin)
  #   # res_form_mixed <- apply_atm_lags(form_mixed)
  #   # res_form_noint <- apply_atm_lags(form_noint)
  #   # res_form_onevar <- apply_atm_lags(form_onevar)
  #   # res_form_onesm <- apply_atm_lags(form_onesm)
  #   # res_form_te <- apply_atm_lags(form_te)
  #   # res_form_mixednl <- apply_atm_lags(form_mixednl)
  #   # res_form_lass <- apply_atm_lags(form_lass)
  #   # res_form_lass2 <- apply_atm_lags(form_lass2)
  #   #
  #   # expect_equal(res_form_lin[[1]] , c("lag1", "lag2"))
  #   # expect_equal(res_form_nlin[[1]] , c("lag1", "lag2"))
  #   # expect_equal(res_form_mixed[[1]] , c("lag1", "lag2"))
  #   # expect_equal(res_form_noint[[1]] , c("lag1", "lag2"))
  #   # expect_equal(res_form_onevar[[1]] , c("lag1"))
  #   # expect_equal(res_form_onesm[[1]] , c("lag1"))
  #   # expect_equal(res_form_te[[1]] , c("lag1", "lag2"))
  #   # expect_equal(res_form_mixednl[[1]] , c("lag1", "lag2", "lag3"))
  #   # expect_equal(res_form_lass[[1]] , c("lag1", "lag2", "lag3"))
  #   # expect_equal(res_form_lass2[[1]] , c("lag1", "lag2", "lag3"))
  #   #
  #   #
  #
  #
  # })
  #
}
