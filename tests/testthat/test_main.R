context("Test deeptrafo")

test_that("simple additive model", {
  
  dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100))
  fml <- y ~ s(x) | z + s(z)
  m <- deeptrafo(fml, dat)
  m %>% fit(epochs = 10)
  plot(m, which_param = "h1")
  plot(m, which_param = "h2")
  coef(m, which = "h1")
  coef(m, which = "h2")
  fitt <- m %>% fitted()
  trf_fun <- m %>% predict.deeptrafo(dat)
  
  a <- trf_fun(dat$y)
  b <- trf_fun(dat$y, type = "pdf")
  c <- trf_fun(dat$y, type = "cdf")
  d <- trf_fun(dat$y, type = "interaction")
  e <- trf_fun(dat$y, type = "shift")
  f <- trf_fun(dat$y, type = "output")
  g <- trf_fun(dat$y, type = "trafo", grid=TRUE)
  h <- trf_fun(dat$y, type = "pdf", grid=TRUE)
  i <- trf_fun(dat$y, type = "cdf", grid=TRUE)
  
})

