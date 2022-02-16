# Demo continuous response
# Lucas Kook
# Feb 2022

set.seed(1234)

# Deps --------------------------------------------------------------------

library(tram)
devtools::load_all("../deepregression/")
devtools::load_all(".")

# Data --------------------------------------------------------------------

data("BostonHousing2", package = "mlbench")

# Model -------------------------------------------------------------------

tm <- Colr(cmedv ~ nox + age, data = BostonHousing2, order = 6)

m <- deeptrafo(cmedv ~ 0 + nox + age, data = BostonHousing2, order = 6,
               optimizer = optimizer_adam(learning_rate = 0.01))
fit(m, epochs = 2e3, validation_split = NULL, batch_size = nrow(BostonHousing2))

coef(tm)
unlist(coef(m, which = "h2"))

# Unconditional case ------------------------------------------------------

tm <- Colr(cmedv ~ 1, data = BostonHousing2, order = 6,
           support = range(BostonHousing2$cmedv))
cfb <- coef(tm, with_baseline = TRUE)

m <- deeptrafo(cmedv ~ 1, data = BostonHousing2, order = 6,
               optimizer = optimizer_adam(learning_rate = 0.001))
coef(m)

tmp <- get_weights(m$model)
tmp[[2]][] <- 0
tmp[[1]][] <- c(cfb[1], log(exp(diff(cfb)) - 1 + 1e-6))
set_weights(m$model, tmp)

coef(m); coef(m, "h2")

logLik(tm)
logLik(m)

m$model$loss(m$init_params$y, fitted(m))

# fit(m, epochs = 3e2, validation_split = NULL, batch_size = nrow(wine))

coef(tm, with_baseline = TRUE)
unlist(coef(m, which = "h1"))
