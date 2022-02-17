# Demo survival response
# Lucas Kook
# Feb 2022

set.seed(1234)

# Deps --------------------------------------------------------------------

library(tram)
devtools::load_all("../deepregression/")
devtools::load_all(".")

# Data --------------------------------------------------------------------

data("GBSG2", package = "TH.data")
GBSG2$surv <- survival::Surv(GBSG2$time, GBSG2$cens)
GBSG2$time <- as.double(GBSG2$time)

# Model -------------------------------------------------------------------

tm <- Colr(surv ~ horTh + age, data = GBSG2, order = 10)

m <- deeptrafo(surv ~ 0 + horTh + age, data = GBSG2,
               optimizer = optimizer_adam(learning_rate = 0.01))
fit(m, epochs = 2e3, validation_split = NULL, batch_size = nrow(GBSG2))

coef(tm)
unlist(coef(m, which = "shifting"))

# Unconditional case ------------------------------------------------------

tm <- Colr(surv ~ 1, data = GBSG2, order = 6, support = range(GBSG2$surv))
cfb <- coef(tm, with_baseline = TRUE)

m <- deeptrafo(surv ~ 1, data = GBSG2, order = 6,
               optimizer = optimizer_adam(learning_rate = 0.001))
coef(m)

tmp <- get_weights(m$model)
tmp[[2]][] <- 0
tmp[[1]][] <- c(cfb[1], log(exp(diff(cfb)) - 1))
set_weights(m$model, tmp)

coef(m); coef(m, "shifting")

logLik(tm)
logLik(m)

m$model$loss(m$init_params$y, fitted(m))

# fit(m, epochs = 3e2, validation_split = NULL, batch_size = nrow(wine))

coef(tm, with_baseline = TRUE)
unlist(coef(m, which = "interacting"))
