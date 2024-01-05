# Demo continuous response penalized
# Lucas Kook
# Feb 2022

set.seed(1234)

# Deps --------------------------------------------------------------------

library(tramnet)
devtools::load_all(".")

# Data --------------------------------------------------------------------

data("BostonHousing2", package = "mlbench")
BostonHousing2$nox <- c(scale(BostonHousing2$nox))
BostonHousing2$age <- c(scale(BostonHousing2$age))

# Model -------------------------------------------------------------------

m0 <- Colr(cmedv ~ 1, data = BostonHousing2, order = 6,
           support = range(BostonHousing2$cmedv))
X <- model.matrix(~ nox + age, data = BostonHousing2)[, -1]
tm <- tramnet(m0, x = X, lambda = 1, alpha = 0)

m <- deeptrafo(cmedv ~ 0 + lasso(nox, la = 1) + lasso(age, la = 1),
               data = BostonHousing2, order = 6,
               optimizer = optimizer_adam(learning_rate = 0.1))
fit(m, epochs = 2e3, validation_split = 0, batch_size = nrow(BostonHousing2),
    callbacks = list(callback_reduce_lr_on_plateau("loss", patience = 20),
                     callback_early_stopping("loss", patience = 50)))

coef(tm, with_baseline = TRUE)
unlist(c(coef(m, which = "interacting"), coef(m, which = "shifting")))

logLik(tm)
logLik(m)
