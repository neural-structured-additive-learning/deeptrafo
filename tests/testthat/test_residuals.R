context("Test deeptrafo residuals")

if (.Platform$OS.type != "windows" &&
  reticulate::py_available() &&
  reticulate::py_module_available("tensorflow") &&
  reticulate::py_module_available("keras") &&
  reticulate::py_module_available("tensorflow_probability")) {
  if (FALSE) {
    test_that("residuals coincide with tram (ordinal)", {
      library("tram")

      data("wine", package = "ordinal")
      tm <- Polr(rating ~ temp + contact, data = wine)

      m <- deeptrafo(rating ~ 0 + temp + contact,
        data = wine,
        optimizer = optimizer_adam(learning_rate = 0.1, decay = 1e-4)
      )
      fit(m,
        epochs = 3e2, validation_split = NULL, batch_size = nrow(wine),
        verbose = FALSE
      )

      expect_equal(residuals(m), residuals(tm), tolerance = 1e-4)
    })

    # test_that("residuals coincide with tram (ordinal, rv)", {
    #   library(tram)
    #   tn <- 1e2
    #   df <- data.frame(y = ordered(sample.int(5, tn, TRUE)),
    #                    x = factor(sample.int(2, tn, TRUE)),
    #                    z = factor(sample.int(2, tn, TRUE)))
    #
    #   tm <- Polr(y | 0 + x ~ z, data = df)
    #   m <- deeptrafo(y | x ~ 0 + z, data = df, optimizer = optimizer_adam(
    #     learning_rate = 0.1, decay = 1e-3))
    #   fit(m, epochs = 3e3, validation_split = NULL, batch_size = tn)
    #
    #   nd <- expand.grid(y = sort(unique(df$y)),
    #                     x = sort(unique(df$x)),
    #                     z = sort(unique(df$z)))
    #   nd$deeptrafo <- residuals(m, newdata = nd)
    #   nd$tram <- residuals(tm, newdata = nd)
    #   nd
    #
    #   expect_equal(residuals(m), residuals(tm), tolerance = 1e-1)
    # })

    test_that("residuals coincide with tram (continuous)", {
      library(tram)

      tm <- Lm(dist ~ speed, data = cars, support = range(cars$dist))
      m <- LmNN(dist ~ 0 + speed, data = cars)

      cfx <- coef(tm, with_baseline = TRUE)
      tmp <- get_weights(m$model)
      tmp[[1]][] <- c(cfx[1], log(exp(cfx[2]) - 1))
      tmp[[2]][] <- -cfx[length(cfx)]
      set_weights(m$model, tmp)

      expect_equal(residuals(m), residuals(tm), tolerance = 1e-4)
    })

    test_that("residuals coincide with tram (continuous, Bernstein)", {
      library(tram)

      tm <- Colr(dist ~ speed, data = cars, support = range(cars$dist))
      m <- ColrNN(dist ~ 0 + speed, data = cars, order = 6)

      cfx <- coef(tm, with_baseline = TRUE)
      tmp <- get_weights(m$model)
      tmp[[1]][] <- c(cfx[1], log(exp(diff(cfx[-length(cfx)])) - 1))
      tmp[[2]][] <- cfx[length(cfx)]
      set_weights(m$model, tmp)

      expect_equal(
        unname(c(unlist(coef(m, which = "interacting")), unlist(coef(m)))),
        unname(cfx),
        tol = 1e-5
      )

      expect_equal(residuals(m), residuals(tm), tolerance = 1e-4)
    })
  }
}
