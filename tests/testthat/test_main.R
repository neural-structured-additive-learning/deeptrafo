context("Test deeptrafo")

test_that("simple additive model", {
  
  n <- 100
  p <- 10

  x <- matrix(runif(n*p), ncol=p)
  
  # univariate
  y <- rnorm(n)
  C <- rep(1, p)
  X <- data.frame(x = x)
  
  data <- c(X, list(C = C))
  
  mod <- deepclasso(
      y = y,
      data = data,
      list_of_formulas = list(loc = ~ classo(x.1, x.2, x.3, x.4, x.5, x.6, x.7, x.8, x.9, x.10, 
                                             C="C", fac=0.2), scale = ~1)
    )

  mod %>% fit(epochs=10)
  expect_is(mod %>% coef(), "list")
  pred <- mod %>% predict(X)
  expect_is(sum(pred), "numeric")
  
})

